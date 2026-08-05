/*
 * Terminus alignment panels.
 *
 * Drawn on a canvas rather than as SVG. A panel can be thousands of rows by
 * hundreds of columns; as SVG that is close to a million DOM nodes and the tab
 * dies. Canvas draws it as run-length blocks in a few milliseconds, and
 * hit-testing for the tooltip is arithmetic rather than event delegation.
 *
 * Two kinds of missing sequence are drawn differently, because they mean
 * different things: a grey cell is a model position the hit did not match
 * although the genome extends there, and a blank cell is sequence that does not
 * exist because the contig ended.
 */

(function () {
  'use strict';

  var T = window.TIRmite;
  if (!T || !T.data.msa || !T.data.msa.length) return;

  var container = document.getElementById('msa-panels');
  if (!container) return;

  var LABEL_WIDTH = 132;
  var AXIS_HEIGHT = 18;
  var MIN_ROW_HEIGHT = 2;
  var MAX_ROW_HEIGHT = 12;
  var MAX_CANVAS_HEIGHT = 4000;

  /*
   * Nucleotide colours come from CSS custom properties so the panel follows
   * the surface into dark mode. The hues were picked by validating four-hue
   * subsets of the report palette for colour-vision separation; see the
   * comment beside --base-a in report.css.
   */
  function theme() {
    var css = getComputedStyle(document.documentElement);
    function prop(name, fallback) {
      return css.getPropertyValue(name).trim() || fallback;
    }
    return {
      pad: prop('--neutral-pad', '#8a8a86'),
      surface: prop('--surface-sunken', '#f4f4f2'),
      ink: prop('--ink-muted', '#75756e'),
      rule: prop('--rule', '#e3e3df'),
      bases: {
        A: prop('--base-a', '#008300'),
        C: prop('--base-c', '#2a78d6'),
        G: prop('--base-g', '#eda100'),
        T: prop('--base-t', '#e87ba4'),
        U: prop('--base-t', '#e87ba4'),
      },
      other: prop('--base-other', '#55554f'),
    };
  }

  /* Expand the run-length pad list into a per-column lookup for one row. */
  function padLookup(row, nCols) {
    var kinds = null;
    (row.pad || []).forEach(function (run) {
      if (!kinds) kinds = new Array(nCols);
      var stop = Math.min(nCols, run[0] + run[1]);
      for (var col = run[0]; col < stop; col++) kinds[col] = run[2];
    });
    return kinds;
  }

  function Panel(spec) {
    this.spec = spec;
    this.rowHeight = Math.max(
      MIN_ROW_HEIGHT,
      Math.min(MAX_ROW_HEIGHT, Math.floor(MAX_CANVAS_HEIGHT / Math.max(1, spec.rows.length)))
    );

    var el = document.createElement('div');
    el.className = 'msa';

    var head = document.createElement('div');
    head.className = 'track-head';
    var name = document.createElement('span');
    name.className = 'track-name';
    name.textContent = spec.model;
    var stats = document.createElement('span');
    stats.className = 'track-stats';
    stats.textContent =
      T.formatInt(spec.n_rows_shown) +
      (spec.n_rows_shown < spec.n_rows_total
        ? ' of ' + T.formatInt(spec.n_rows_total)
        : '') +
      ' hits · ' +
      T.formatInt(spec.n_cols) +
      ' columns · ' +
      (spec.aligner === 'mafft'
        ? 'aligned with MAFFT'
        : spec.aligner === 'anchor'
          ? 'placed by model position'
          : 'left-aligned');
    head.appendChild(name);
    head.appendChild(stats);
    el.appendChild(head);

    var wrap = document.createElement('div');
    wrap.className = 'msa-canvas-wrap';
    var canvas = document.createElement('canvas');
    canvas.setAttribute(
      'aria-label',
      'Alignment of ' + spec.n_rows_shown + ' ' + spec.model +
        ' hits. Per-hit values are listed in the tables below.'
    );
    wrap.appendChild(canvas);
    el.appendChild(wrap);

    var caption = document.createElement('div');
    caption.className = 'note';
    caption.style.marginTop = '8px';
    caption.textContent = spec.note || '';
    if (spec.n_rows_shown < spec.n_rows_total) {
      caption.textContent +=
        (caption.textContent ? ' ' : '') +
        'Showing the ' + T.formatInt(spec.n_rows_shown) +
        ' most significant of ' + T.formatInt(spec.n_rows_total) + ' hits.';
    }
    if (caption.textContent) el.appendChild(caption);

    this.el = el;
    this.canvas = canvas;
    this.wrap = wrap;
    this.bind();
  }

  Panel.prototype.draw = function () {
    var spec = this.spec;
    var rows = spec.rows;
    var nCols = spec.n_cols || 1;
    var colours = theme();

    var width = Math.max(320, this.wrap.clientWidth || 640);
    var plotWidth = width - LABEL_WIDTH;
    var colWidth = plotWidth / nCols;
    var height = AXIS_HEIGHT + rows.length * this.rowHeight;

    var ratio = window.devicePixelRatio || 1;
    this.canvas.width = Math.round(width * ratio);
    this.canvas.height = Math.round(height * ratio);
    this.canvas.style.width = width + 'px';
    this.canvas.style.height = height + 'px';

    var ctx = this.canvas.getContext('2d');
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    ctx.clearRect(0, 0, width, height);

    // Coverage bar: how many rows have real sequence in each column. It is what
    // makes the panel readable when a row is one pixel tall.
    var coverage = new Array(nCols).fill(0);
    for (var r = 0; r < rows.length; r++) {
      var seq = rows[r].seq;
      for (var c = 0; c < nCols; c++) {
        if (seq.charAt(c) && seq.charAt(c) !== '-') coverage[c]++;
      }
    }
    var maxCoverage = Math.max(1, rows.length);
    ctx.fillStyle = colours.ink;
    for (c = 0; c < nCols; c++) {
      var barHeight = (coverage[c] / maxCoverage) * (AXIS_HEIGHT - 4);
      ctx.globalAlpha = 0.45;
      ctx.fillRect(
        LABEL_WIDTH + c * colWidth,
        AXIS_HEIGHT - 2 - barHeight,
        Math.max(0.6, colWidth),
        barHeight
      );
    }
    ctx.globalAlpha = 1;

    ctx.font = '10px ui-monospace, SFMono-Regular, Menlo, monospace';
    ctx.textBaseline = 'middle';

    for (r = 0; r < rows.length; r++) {
      var row = rows[r];
      var y = AXIS_HEIGHT + r * this.rowHeight;
      var pads = padLookup(row, nCols);

      if (this.rowHeight >= 9) {
        var hit = T.hitByUid(row.uid);
        ctx.fillStyle = colours.ink;
        ctx.fillText(
          hit ? hit.contig.name + ':' + hit.start : 'hit ' + row.uid,
          6,
          y + this.rowHeight / 2,
          LABEL_WIDTH - 12
        );
      }

      // Run-length: aligned sequence repeats heavily, so filling one rect per
      // run rather than per cell cuts the draw call count by an order of
      // magnitude.
      var runStart = 0;
      var runColour = cellColour(row.seq.charAt(0), pads && pads[0], colours);
      for (c = 1; c <= nCols; c++) {
        var colour = c === nCols
          ? null
          : cellColour(row.seq.charAt(c), pads && pads[c], colours);
        if (colour !== runColour || c === nCols) {
          if (runColour) {
            ctx.fillStyle = runColour;
            ctx.fillRect(
              LABEL_WIDTH + runStart * colWidth,
              y,
              Math.max(0.6, (c - runStart) * colWidth),
              Math.max(1, this.rowHeight - (this.rowHeight >= 4 ? 1 : 0))
            );
          }
          runStart = c;
          runColour = colour;
        }
      }

      // Letters wherever the cells are big enough to hold one. Colour alone
      // would otherwise be the sole encoding of base identity.
      if (this.rowHeight >= 8 && colWidth >= 6) {
        ctx.textAlign = 'center';
        for (c = 0; c < nCols; c++) {
          var ch = row.seq.charAt(c);
          if (!ch || ch === '-' || ch === '.') continue;
          if (pads && pads[c] === 'm') continue;
          var fill = cellColour(ch, null, colours);
          // The palette spans a wide lightness range, so a fixed ink colour
          // would be unreadable on half of it.
          ctx.fillStyle = inkOn(fill);
          ctx.fillText(
            ch.toUpperCase(),
            LABEL_WIDTH + (c + 0.5) * colWidth,
            y + this.rowHeight / 2
          );
        }
        ctx.textAlign = 'left';
      }
    }

    ctx.strokeStyle = colours.rule;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(LABEL_WIDTH, AXIS_HEIGHT - 0.5);
    ctx.lineTo(width, AXIS_HEIGHT - 0.5);
    ctx.stroke();

    this.colWidth = colWidth;
  };

  /*
   * Pick black or white for a letter drawn on `background`, whichever has the
   * higher WCAG contrast ratio against it.
   */
  function inkOn(background) {
    if (!background || background.charAt(0) !== '#' || background.length < 7) {
      return '#ffffff';
    }
    var channels = [1, 3, 5].map(function (i) {
      var c = parseInt(background.substr(i, 2), 16) / 255;
      return c <= 0.04045 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4);
    });
    var luminance =
      0.2126 * channels[0] + 0.7152 * channels[1] + 0.0722 * channels[2];
    var onWhite = 1.05 / (luminance + 0.05);
    var onBlack = (luminance + 0.05) / 0.05;
    return onBlack >= onWhite ? '#111111' : '#ffffff';
  }

  function cellColour(char, pad, colours) {
    // A model pad is real sequence the alignment did not claim: grey, not
    // blank. A gap is sequence that does not exist: nothing is drawn.
    if (pad === 'm') return colours.pad;
    if (!char || char === '-' || char === '.') return null;
    return colours.bases[char.toUpperCase()] || colours.other;
  }

  Panel.prototype.bind = function () {
    var self = this;
    var tooltip = document.getElementById('tooltip');

    this.canvas.addEventListener('mousemove', function (event) {
      if (!tooltip) return;
      var box = self.canvas.getBoundingClientRect();
      var x = event.clientX - box.left;
      var y = event.clientY - box.top;
      var rowIndex = Math.floor((y - AXIS_HEIGHT) / self.rowHeight);
      if (rowIndex < 0 || rowIndex >= self.spec.rows.length || x < LABEL_WIDTH) {
        tooltip.dataset.show = '0';
        return;
      }
      var row = self.spec.rows[rowIndex];
      var col = Math.floor((x - LABEL_WIDTH) / self.colWidth);
      showRowTooltip(tooltip, self.spec, row, col, event);
    });

    this.canvas.addEventListener('mouseleave', function () {
      if (tooltip) tooltip.dataset.show = '0';
    });

    this.canvas.addEventListener('click', function (event) {
      var box = self.canvas.getBoundingClientRect();
      var rowIndex = Math.floor((event.clientY - box.top - AXIS_HEIGHT) / self.rowHeight);
      if (rowIndex < 0 || rowIndex >= self.spec.rows.length) return;
      var hit = T.hitByUid(self.spec.rows[rowIndex].uid);
      if (!hit) return;
      // Bring the corresponding annotation into view rather than duplicating
      // the track browser's controls here.
      var track = document.querySelector(
        '.track[data-contig="' + cssEscape(hit.contig.name) + '"]'
      );
      if (!track || !track.__track) return;
      track.__track.built = true;
      track.__track.setView(hit.start - hit.length * 3, hit.end + hit.length * 3);
      track.__track.draw();
      track.scrollIntoView({ block: 'center', behavior: 'smooth' });
    });
  };

  function cssEscape(value) {
    return window.CSS && CSS.escape ? CSS.escape(value) : value.replace(/"/g, '\\"');
  }

  function showRowTooltip(tooltip, spec, row, col, event) {
    var hit = T.hitByUid(row.uid);
    tooltip.innerHTML = '';

    var title = document.createElement('h4');
    var group = row.group_i == null ? null : T.data.groups[row.group_i];
    if (group) {
      var swatch = document.createElement('span');
      swatch.className = 'swatch';
      swatch.style.background = T.groupColour(group);
      title.appendChild(swatch);
    }
    title.appendChild(document.createTextNode(spec.model));
    tooltip.appendChild(title);

    var dl = document.createElement('dl');
    function add(term, value) {
      if (value == null || value === '') return;
      var dt = document.createElement('dt');
      dt.textContent = term;
      var dd = document.createElement('dd');
      dd.textContent = value;
      dl.appendChild(dt);
      dl.appendChild(dd);
    }
    if (hit) {
      add('Location', hit.contig.name + ':' + hit.start.toLocaleString() + '–' +
        hit.end.toLocaleString());
      add('Strand', hit.strand);
      add('E-value', T.formatEvalue(hit.evalue));
    }
    add('Terminus', row.role || 'unpaired');
    add('Column', (col + 1).toLocaleString() + ' of ' + spec.n_cols.toLocaleString());

    var char = row.seq.charAt(col);
    var kind = padKindAt(row, col);
    if (kind === 'm') {
      add('At this column', 'model not matched here — genome continues');
    } else if (!char || char === '-') {
      add('At this column', 'no sequence — contig ends');
    } else {
      add('At this column', char);
    }
    tooltip.appendChild(dl);

    var hint = document.createElement('div');
    hint.className = 'hint';
    hint.textContent = 'Click to show this hit on its contig.';
    tooltip.appendChild(hint);

    tooltip.dataset.show = '1';
    var box = tooltip.getBoundingClientRect();
    var x = event.clientX + 14;
    var y = event.clientY + 14;
    if (x + box.width > window.innerWidth - 8) x = event.clientX - box.width - 14;
    if (y + box.height > window.innerHeight - 8) y = event.clientY - box.height - 14;
    tooltip.style.left = Math.max(8, x + window.scrollX) + 'px';
    tooltip.style.top = Math.max(8, y + window.scrollY) + 'px';
  }

  function padKindAt(row, col) {
    var runs = row.pad || [];
    for (var i = 0; i < runs.length; i++) {
      if (col >= runs[i][0] && col < runs[i][0] + runs[i][1]) return runs[i][2];
    }
    return null;
  }

  var panels = T.data.msa.map(function (spec) {
    var panel = new Panel(spec);
    container.appendChild(panel.el);
    panel.draw();
    return panel;
  });

  var resizeTimer = null;
  window.addEventListener('resize', function () {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(function () {
      panels.forEach(function (panel) {
        panel.draw();
      });
    }, 150);
  });

  if (window.matchMedia) {
    var scheme = window.matchMedia('(prefers-color-scheme: dark)');
    if (scheme.addEventListener) {
      scheme.addEventListener('change', function () {
        panels.forEach(function (panel) {
          panel.draw();
        });
      });
    }
  }
})();

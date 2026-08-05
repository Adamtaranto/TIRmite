/*
 * Genome-browser-style annotation tracks.
 *
 * One SVG per contig, drawn against a base-pair viewport rather than an SVG
 * viewBox. Scaling a viewBox would scale glyph heights, stroke widths and
 * arrowheads along with the axis, so the geometry is recomputed in screen
 * pixels on every redraw instead.
 *
 * Only the hits inside the viewport are drawn. Contigs are sorted so their hit
 * columns form contiguous slices, so finding the visible window is a binary
 * search and a redraw costs O(visible) regardless of how many hits the run
 * produced.
 */

(function () {
  'use strict';

  var T = window.TIRmite;
  if (!T) return;

  var data = T.data;
  var ROW_HEIGHT = 18;
  var GLYPH_HEIGHT = 11;
  var AXIS_HEIGHT = 22;
  var PAD_TOP = 8;
  var MIN_GLYPH_PX = 2;
  var MIN_VIEW_BP = 40;
  /* Long hits that started before the viewport still have to be drawn. */
  var LOOKBACK = 512;

  var container = document.getElementById('tracks');
  if (!container) return;

  var tooltip = document.getElementById('tooltip');
  var modal = document.getElementById('element-modal');

  var tracks = [];
  var showEmpty = false;
  var filterText = '';

  /* ---- geometry ------------------------------------------------------- */

  /*
   * One path per hit carrying every cue at once: direction as an arrowhead,
   * rounded corners, and a jagged edge wherever the model was not fully
   * matched. Drawing these as separate elements would triple the node count
   * for no gain.
   */
  function glyphPath(x0, x1, y, h, opts) {
    var w = Math.max(MIN_GLYPH_PX, x1 - x0);
    x1 = x0 + w;
    var r = Math.min(3, h / 2, w / 4);
    var head = Math.min(7, w / 3);
    var forward = opts.strand !== '-';
    var y1 = y + h;
    var mid = y + h / 2;
    var p = [];

    // Leading edge: the arrow point on the strand's 3' side, otherwise a
    // rounded or jagged blunt end.
    if (forward) {
      p.push('M' + x0 + ',' + y);
      edge(p, x0, y, y1, opts.jagLeft, r, 'left');
      p.push('L' + (x1 - head) + ',' + y1);
      p.push('L' + x1 + ',' + mid);
      p.push('L' + (x1 - head) + ',' + y);
      p.push('Z');
    } else {
      p.push('M' + x1 + ',' + y);
      p.push('L' + (x0 + head) + ',' + y);
      p.push('L' + x0 + ',' + mid);
      p.push('L' + (x0 + head) + ',' + y1);
      p.push('L' + x1 + ',' + y1);
      edge(p, x1, y1, y, opts.jagRight, r, 'right');
      p.push('Z');
    }
    return p.join(' ');
  }

  /*
   * Draw the blunt end of a glyph, from (x, from) to (x, to). A jagged edge
   * means the hit did not match the whole model there; a rounded one means it
   * did, or that we cannot tell.
   */
  function edge(p, x, from, to, jagged, r, side) {
    if (!jagged) {
      var dir = to > from ? 1 : -1;
      var inset = side === 'left' ? r : -r;
      p.push('L' + (x - inset) + ',' + from);
      p.push('Q' + x + ',' + from + ' ' + x + ',' + (from + r * dir));
      p.push('L' + x + ',' + (to - r * dir));
      p.push('Q' + x + ',' + to + ' ' + (x - inset) + ',' + to);
      return;
    }
    // Teeth alternate either side of the edge line so it reads as torn rather
    // than as a nibble taken out of the glyph.
    var amp = 3;
    var steps = 6;
    p.push('L' + x + ',' + from);
    for (var i = 1; i < steps; i++) {
      var ty = from + ((to - from) * i) / steps;
      p.push('L' + (x + (i % 2 ? amp : -amp)) + ',' + ty);
    }
    p.push('L' + x + ',' + to);
  }

  /* ---- one contig track ----------------------------------------------- */

  function Track(contig) {
    this.contig = contig;
    this.view = { start: 1, end: contig.length };
    this.built = false;
    this.frame = null;

    var el = document.createElement('div');
    el.className = 'track';
    el.dataset.contig = contig.name;

    var head = document.createElement('div');
    head.className = 'track-head';
    var name = document.createElement('span');
    name.className = 'track-name';
    name.textContent = contig.name;
    var stats = document.createElement('span');
    stats.className = 'track-stats';
    stats.textContent =
      T.formatBp(contig.length) +
      (contig.length_source === 'inferred' ? ' (estimated)' : '') +
      ' · ' +
      T.formatInt(contig.n_hits) +
      ' hits · ' +
      T.formatInt(contig.n_pairs) +
      ' pairs';
    head.appendChild(name);
    head.appendChild(stats);

    var self = this;
    var bar = document.createElement('div');
    bar.className = 'track-head';
    var viewLabel = document.createElement('span');
    viewLabel.className = 'track-view';
    var buttons = document.createElement('span');
    // Wrapped rather than bound: a bound factor would leave the click event
    // arriving as the zoom anchor, which is how these buttons produced NaN.
    buttons.appendChild(
      this.button('−', 'Zoom out', function () {
        self.zoomBy(2);
      })
    );
    buttons.appendChild(
      this.button('+', 'Zoom in', function () {
        self.zoomBy(0.5);
      })
    );
    buttons.appendChild(
      this.button('Reset', 'Show whole contig', function () {
        self.reset();
      })
    );
    bar.appendChild(viewLabel);
    bar.appendChild(buttons);

    // The glyphs scroll vertically inside a fixed-height body; the axis lives
    // in its own SVG below it so it stays visible on a tall track.
    var body = document.createElement('div');
    body.className = 'track-body';

    var svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('role', 'img');
    svg.setAttribute(
      'aria-label',
      'Terminus hits on ' + contig.name + '. Values are also listed in the tables below.'
    );
    svg.setAttribute('tabindex', '0');
    svg.style.height = contig.n_rows * ROW_HEIGHT + PAD_TOP * 2 + 'px';
    body.appendChild(svg);

    var axisSvg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    axisSvg.setAttribute('class', 'track-axis');
    axisSvg.setAttribute('aria-hidden', 'true');
    axisSvg.style.height = AXIS_HEIGHT + 'px';

    // A real horizontal scrollbar, so panning is discoverable without knowing
    // that the track can be dragged. The spacer's width encodes the zoom: at
    // whole-contig view it matches the viewport and the bar disappears.
    var hscroll = document.createElement('div');
    hscroll.className = 'track-hscroll';
    hscroll.setAttribute('aria-label', 'Scroll along ' + contig.name);
    var spacer = document.createElement('div');
    spacer.className = 'track-hspacer';
    hscroll.appendChild(spacer);

    el.appendChild(head);
    el.appendChild(bar);
    el.appendChild(body);
    el.appendChild(axisSvg);
    el.appendChild(hscroll);

    this.el = el;
    this.body = body;
    this.svg = svg;
    this.axisSvg = axisSvg;
    this.hscroll = hscroll;
    this.spacer = spacer;
    this.viewLabel = viewLabel;
    this.expectedScrollLeft = null;
    this.bind();
  }

  Track.prototype.button = function (label, title, handler) {
    var b = document.createElement('button');
    b.type = 'button';
    b.textContent = label;
    b.title = title;
    b.setAttribute('aria-label', title + ' on this contig');
    b.addEventListener('click', handler);
    return b;
  };

  Track.prototype.bpAt = function (clientX) {
    var box = this.svg.getBoundingClientRect();
    var frac = (clientX - box.left) / Math.max(1, box.width);
    return this.view.start + frac * (this.view.end - this.view.start);
  };

  Track.prototype.bind = function () {
    var self = this;

    this.svg.addEventListener(
      'wheel',
      function (event) {
        // A horizontal swipe means pan; a vertical one means zoom. Reading
        // only deltaY treated a pure horizontal swipe as deltaY === 0 and
        // zoomed in on every event, which is the opposite of what the gesture
        // asks for. Pinch-zoom arrives as ctrl+wheel and stays a zoom.
        var horizontal =
          !event.ctrlKey && Math.abs(event.deltaX) > Math.abs(event.deltaY);

        if (horizontal) {
          if (!event.deltaX) return;
          event.preventDefault();
          var span = self.view.end - self.view.start;
          var box = self.svg.getBoundingClientRect();
          var shift = (event.deltaX / Math.max(1, box.width)) * span;
          self.setView(self.view.start + shift, self.view.end + shift);
          return;
        }

        if (!event.deltaY) return;
        event.preventDefault();
        var anchor = self.bpAt(event.clientX);
        self.zoomBy(event.deltaY > 0 ? 1.18 : 1 / 1.18, anchor);
      },
      { passive: false }
    );

    var dragging = null;
    this.svg.addEventListener('pointerdown', function (event) {
      if (event.button !== 0) return;
      dragging = { x: event.clientX, start: self.view.start, end: self.view.end };
      self.svg.setPointerCapture(event.pointerId);
      self.svg.classList.add('dragging');
    });
    this.svg.addEventListener('pointermove', function (event) {
      if (dragging) {
        var box = self.svg.getBoundingClientRect();
        var span = dragging.end - dragging.start;
        var shift = ((dragging.x - event.clientX) / Math.max(1, box.width)) * span;
        self.setView(dragging.start + shift, dragging.end + shift);
        return;
      }
      self.hover(event);
    });
    function endDrag(event) {
      if (!dragging) return;
      dragging = null;
      self.svg.classList.remove('dragging');
      try {
        self.svg.releasePointerCapture(event.pointerId);
      } catch (err) {
        /* capture may already be gone */
      }
    }
    this.svg.addEventListener('pointerup', endDrag);
    this.svg.addEventListener('pointercancel', endDrag);
    this.svg.addEventListener('pointerleave', function () {
      hideTooltip();
      self.unhighlight();
    });

    // Dragging the track and dragging the scrollbar have to agree, so each
    // one's update must not be read back as a fresh user action. The test is
    // "is this the position we just wrote?" rather than a timed flag: a flag
    // cleared on requestAnimationFrame stays set forever in a background tab,
    // which would kill the scrollbar outright.
    this.hscroll.addEventListener('scroll', function () {
      if (
        self.expectedScrollLeft !== null &&
        Math.abs(self.hscroll.scrollLeft - self.expectedScrollLeft) < 1
      ) {
        return;
      }
      var span = self.view.end - self.view.start;
      var travel = self.spacer.offsetWidth - self.hscroll.clientWidth;
      if (travel <= 0) return;
      var frac = self.hscroll.scrollLeft / travel;
      var start = 1 + frac * Math.max(0, self.contig.length - span);
      self.view = { start: start, end: start + span };
      self.schedule();
    });

    // Clicking opens the feature popup. A hover tooltip cannot carry working
    // buttons -- moving the pointer towards it leaves the glyph that produced
    // it -- so anything the reader needs to act on lives in the popup.
    this.svg.addEventListener('click', function (event) {
      var target = event.target.closest('.glyph');
      if (!target) return;
      var h = T.hitByUid(+target.dataset.uid);
      if (!h) return;
      hideTooltip();
      openFeature(h);
    });

    this.svg.addEventListener('dblclick', function (event) {
      var target = event.target.closest('.glyph');
      if (!target) return;
      var h = T.hitByUid(+target.dataset.uid);
      if (!h) return;
      var pad;
      if (h.element) {
        pad = Math.max(50, Math.round(h.element.length * 0.1));
        self.setView(h.element.start - pad, h.element.end + pad);
      } else {
        pad = Math.max(50, Math.round(h.length * 2));
        self.setView(h.start - pad, h.end + pad);
      }
    });

    /* Keyboard equivalents, so the track is not mouse-only. */
    this.svg.addEventListener('keydown', function (event) {
      var span = self.view.end - self.view.start;
      var step = span * 0.15;
      var handled = true;
      switch (event.key) {
        case 'ArrowLeft':
          self.setView(self.view.start - step, self.view.end - step);
          break;
        case 'ArrowRight':
          self.setView(self.view.start + step, self.view.end + step);
          break;
        case '+':
        case '=':
          self.zoomBy(0.5);
          break;
        case '-':
        case '_':
          self.zoomBy(2);
          break;
        case 'Home':
          self.reset();
          break;
        default:
          handled = false;
      }
      if (handled) event.preventDefault();
    });
  };

  Track.prototype.setView = function (start, end) {
    var length = this.contig.length;
    var span = Math.min(Math.max(end - start, MIN_VIEW_BP), length);
    if (start < 1) start = 1;
    if (start + span > length + 1) start = Math.max(1, length + 1 - span);
    this.view = { start: start, end: start + span };
    this.schedule();
  };

  Track.prototype.zoomBy = function (factor, anchorBp) {
    var span = this.view.end - this.view.start;
    var anchor = anchorBp == null ? this.view.start + span / 2 : anchorBp;
    var frac = (anchor - this.view.start) / span;
    var next = span * factor;
    this.setView(anchor - frac * next, anchor - frac * next + next);
  };

  Track.prototype.reset = function () {
    this.setView(1, this.contig.length);
  };

  /* Redraws are coalesced so a burst of wheel events costs one paint. */
  Track.prototype.schedule = function () {
    var self = this;
    if (this.frame) return;
    this.frame = requestAnimationFrame(function () {
      self.frame = null;
      self.draw();
    });
  };

  Track.prototype.visibleHits = function () {
    var lo = this.contig.hit_lo;
    var hi = this.contig.hit_hi;
    var from = T.firstStartingAt(lo, hi, this.view.start);
    // Walk back a bounded distance: a hit that starts before the viewport can
    // still extend into it, and the columns are ordered by start, not by end.
    from = Math.max(lo, from - LOOKBACK);
    var out = [];
    for (var i = from; i < hi; i++) {
      if (data.hits.start[i] > this.view.end) break;
      if (data.hits.end[i] < this.view.start) continue;
      out.push(i);
    }
    return out;
  };

  Track.prototype.draw = function () {
    var width = this.svg.clientWidth;
    if (!width) return;
    var view = this.view;
    var span = view.end - view.start;

    function x(bp) {
      return ((bp - view.start) / span) * width;
    }

    this.viewLabel.textContent =
      this.contig.name +
      ':' +
      Math.round(view.start).toLocaleString() +
      '–' +
      Math.round(view.end).toLocaleString() +
      '  (' +
      T.formatBp(Math.round(span)) +
      ' shown)';

    while (this.svg.firstChild) this.svg.removeChild(this.svg.firstChild);

    var top = PAD_TOP;
    var indices = this.visibleHits();

    // Links first, so they sit behind the glyphs they connect.
    var linkGroup = svgEl('g', { class: 'links' });
    var glyphGroup = svgEl('g', { class: 'glyphs' });
    var drawnElements = {};

    indices.forEach(function (i) {
      var h = T.hit(i);
      var y = top + h.row * ROW_HEIGHT;
      var group = h.groups[0] || null;
      var colour = T.groupColour(group) || 'var(--ink-muted)';

      if (h.element && !drawnElements[h.element.pair_id]) {
        drawnElements[h.element.pair_id] = true;
        var partnerUid =
          h.element.left_uid === h.uid ? h.element.right_uid : h.element.left_uid;
        var partner = T.hitByUid(partnerUid);
        if (partner) {
          linkGroup.appendChild(link(h, partner, x, top, colour));
        }
      }

      var jagLeft = !!h.truncLeft && !h.clipLeft;
      var jagRight = !!h.truncRight && !h.clipRight;
      var path = svgEl('path', {
        class: 'glyph' + (h.element ? '' : ' unpaired'),
        d: glyphPath(x(h.start), x(h.end + 1), y, GLYPH_HEIGHT, {
          strand: h.strand,
          jagLeft: jagLeft,
          jagRight: jagRight,
        }),
        fill: colour,
        'data-uid': String(h.uid),
      });
      glyphGroup.appendChild(path);

      // A contig-end cap is a different statement from a jagged edge: the
      // model is unmatched because the sequence stopped, not because the hit
      // was partial.
      if (h.clipLeft) glyphGroup.appendChild(cap(x(h.start), y));
      if (h.clipRight) glyphGroup.appendChild(cap(x(h.end + 1), y));

      if (h.overflow) {
        glyphGroup.appendChild(
          svgEl('circle', {
            class: 'cap',
            cx: x(h.start),
            cy: y - 2,
            r: 1.5,
            fill: 'var(--neutral-pad)',
          })
        );
      }
    });

    this.svg.appendChild(linkGroup);
    this.svg.appendChild(glyphGroup);

    if (!indices.length) {
      var note = svgEl('text', {
        class: 'overflow-note',
        x: width / 2,
        y: top + 16,
        'text-anchor': 'middle',
      });
      note.textContent = 'No hits in this window';
      this.svg.appendChild(note);
    }

    // The axis has its own SVG below the scrolling body, so it stays put on a
    // track too tall to show at once.
    while (this.axisSvg.firstChild) {
      this.axisSvg.removeChild(this.axisSvg.firstChild);
    }
    this.axisSvg.appendChild(axis(width, 1, view, span));

    this.updateScrollProxy();
  };

  /*
   * Resize the horizontal scrollbar to reflect the current zoom, and move its
   * thumb to match the viewport.
   */
  Track.prototype.updateScrollProxy = function () {
    var span = this.view.end - this.view.start;
    var viewport = this.hscroll.clientWidth;
    if (!viewport) return;

    var ratio = this.contig.length / Math.max(1, span);
    this.spacer.style.width = Math.round(viewport * ratio) + 'px';

    var travel = this.spacer.offsetWidth - viewport;
    if (travel <= 0) {
      this.expectedScrollLeft = null;
      return;
    }
    var room = Math.max(1, this.contig.length - span);
    var frac = Math.min(1, Math.max(0, (this.view.start - 1) / room));

    // Recorded, not flagged: the scroll handler compares against this to tell
    // its own echo from a real user scroll.
    this.expectedScrollLeft = frac * travel;
    this.hscroll.scrollLeft = this.expectedScrollLeft;
  };

  function cap(x, y) {
    return svgEl('line', {
      class: 'cap',
      x1: x,
      x2: x,
      y1: y - 2,
      y2: y + GLYPH_HEIGHT + 2,
    });
  }

  function link(a, b, x, top, colour) {
    var left = a.start <= b.start ? a : b;
    var right = left === a ? b : a;
    var y1 = top + left.row * ROW_HEIGHT + GLYPH_HEIGHT / 2;
    var y2 = top + right.row * ROW_HEIGHT + GLYPH_HEIGHT / 2;
    var x1 = x(left.end + 1);
    var x2 = x(right.start);
    var d;
    if (left.row === right.row) {
      d = 'M' + x1 + ',' + y1 + ' L' + x2 + ',' + y2;
    } else {
      var mx = (x1 + x2) / 2;
      d = 'M' + x1 + ',' + y1 + ' C' + mx + ',' + y1 + ' ' + mx + ',' + y2 + ' ' + x2 + ',' + y2;
    }
    return svgEl('path', {
      class: 'link',
      d: d,
      stroke: colour,
      'data-pair': left.element ? left.element.pair_id : '',
    });
  }

  /* Ticks at a round interval, chosen so about six fit across the view. */
  function axis(width, y, view, span) {
    var g = svgEl('g', { class: 'axis' });
    g.appendChild(svgEl('line', { x1: 0, x2: width, y1: y, y2: y }));

    var raw = span / 6;
    var magnitude = Math.pow(10, Math.floor(Math.log10(raw)));
    var step = magnitude;
    [1, 2, 5, 10].some(function (m) {
      if (magnitude * m >= raw) {
        step = magnitude * m;
        return true;
      }
      return false;
    });

    for (var bp = Math.ceil(view.start / step) * step; bp <= view.end; bp += step) {
      var px = ((bp - view.start) / span) * width;
      g.appendChild(svgEl('line', { x1: px, x2: px, y1: y, y2: y + 4 }));
      var label = svgEl('text', {
        x: px,
        y: y + 15,
        'text-anchor': px < 24 ? 'start' : px > width - 24 ? 'end' : 'middle',
      });
      label.textContent = tickLabel(bp, step);
      g.appendChild(label);
    }
    return g;
  }

  /*
   * Tick labels take their unit and precision from the tick interval, not from
   * the magnitude of the coordinate. Formatting by magnitude alone makes every
   * tick read "1.36 Mb" once zoomed in past a megabase.
   */
  function tickLabel(bp, step) {
    var unit = step >= 1e6 ? 1e6 : step >= 1e3 ? 1e3 : 1;
    var suffix = unit === 1e6 ? ' Mb' : unit === 1e3 ? ' kb' : '';
    if (unit === 1) return Math.round(bp).toLocaleString();
    var decimals = Math.min(3, Math.max(0, Math.ceil(-Math.log10(step / unit))));
    return (bp / unit).toFixed(decimals) + suffix;
  }

  function svgEl(name, attrs) {
    var el = document.createElementNS('http://www.w3.org/2000/svg', name);
    for (var key in attrs) {
      if (Object.prototype.hasOwnProperty.call(attrs, key)) {
        el.setAttribute(key, attrs[key]);
      }
    }
    return el;
  }

  /* ---- hover ---------------------------------------------------------- */

  Track.prototype.hover = function (event) {
    var target = event.target.closest ? event.target.closest('.glyph') : null;
    if (!target) {
      hideTooltip();
      this.unhighlight();
      return;
    }
    var h = T.hitByUid(+target.dataset.uid);
    if (!h) return;
    this.highlight(h);
    showTooltip(h, event);
  };

  Track.prototype.highlight = function (h) {
    this.unhighlight();
    var uids = [h.uid];
    if (h.element) uids.push(h.element.left_uid, h.element.right_uid);
    uids.forEach(function (uid) {
      var el = this.svg.querySelector('.glyph[data-uid="' + uid + '"]');
      if (el) el.classList.add('highlight');
    }, this);
    if (h.element) {
      var line = this.svg.querySelector(
        '.link[data-pair="' + cssEscape(h.element.pair_id) + '"]'
      );
      if (line) line.classList.add('highlight');
    }
  };

  Track.prototype.unhighlight = function () {
    var marked = this.svg.querySelectorAll('.highlight');
    for (var i = 0; i < marked.length; i++) marked[i].classList.remove('highlight');
  };

  function cssEscape(value) {
    return window.CSS && CSS.escape ? CSS.escape(value) : value.replace(/"/g, '\\"');
  }

  function row(dl, term, value) {
    if (value == null || value === '') return;
    var dt = document.createElement('dt');
    dt.textContent = term;
    var dd = document.createElement('dd');
    dd.textContent = value;
    dl.appendChild(dt);
    dl.appendChild(dd);
  }

  function showTooltip(h, event) {
    if (!tooltip) return;
    tooltip.innerHTML = '';

    var title = document.createElement('h4');
    var group = h.groups[0];
    if (group) {
      var swatch = document.createElement('span');
      swatch.className = 'swatch';
      swatch.style.background = T.groupColour(group);
      title.appendChild(swatch);
    }
    title.appendChild(document.createTextNode(h.model.name));
    tooltip.appendChild(title);

    var dl = document.createElement('dl');
    row(
      dl,
      'Location',
      h.contig.name + ':' + h.start.toLocaleString() + '–' + h.end.toLocaleString()
    );
    row(dl, 'Length', T.formatBp(h.length));
    row(dl, 'Strand', h.strand);
    row(dl, 'E-value', T.formatEvalue(h.evalue));
    row(dl, 'Score', h.score == null ? null : String(h.score));

    if (h.hmmStart != null && h.hmmEnd != null) {
      row(
        dl,
        'Model match',
        h.hmmStart +
          '–' +
          h.hmmEnd +
          (h.model.length ? ' of ' + h.model.length : '') +
          (h.modelCoverage != null ? ' (' + T.formatPercent(h.modelCoverage) + ')' : '')
      );
    } else {
      row(dl, 'Model match', 'coordinates unavailable for this input format');
    }
    if (h.spanCoverage != null) {
      row(dl, 'Span / model', h.spanCoverage.toFixed(2) + '×');
    }
    row(dl, 'Pairing group', h.groups.map(function (g) { return g.label; }).join(', ') || 'none');

    if (h.element) {
      var partnerUid =
        h.element.left_uid === h.uid ? h.element.right_uid : h.element.left_uid;
      var partner = T.hitByUid(partnerUid);
      row(dl, 'Terminus', h.role === 'left' ? 'left' : h.role === 'right' ? 'right' : '–');
      if (partner) {
        row(
          dl,
          'Partner',
          partner.model.name +
            ' at ' +
            partner.start.toLocaleString() +
            '–' +
            partner.end.toLocaleString()
        );
      }
      row(dl, 'Between termini', T.formatBp(h.element.inner_distance));
      row(dl, 'Element', h.element.element_id + ' · ' + T.formatBp(h.element.length));
    } else {
      row(dl, 'Pairing', 'unpaired');
    }
    tooltip.appendChild(dl);

    var notes = [T.truncationNote(h, 'left'), T.truncationNote(h, 'right')].filter(
      Boolean
    );
    if (h.overflow) {
      notes.push('Track row limit reached; this hit shares a row with others.');
    }
    notes.push(
      h.element
        ? 'Click for details and the element sequence.'
        : 'Click for details.'
    );
    if (notes.length) {
      var hint = document.createElement('div');
      hint.className = 'hint';
      hint.textContent = notes.join(' ');
      tooltip.appendChild(hint);
    }

    tooltip.dataset.show = '1';
    position(tooltip, event);
  }

  function position(el, event) {
    var pad = 14;
    var box = el.getBoundingClientRect();
    var x = event.clientX + pad;
    var y = event.clientY + pad;
    if (x + box.width > window.innerWidth - 8) x = event.clientX - box.width - pad;
    if (y + box.height > window.innerHeight - 8) y = event.clientY - box.height - pad;
    el.style.left = Math.max(8, x + window.scrollX) + 'px';
    el.style.top = Math.max(8, y + window.scrollY) + 'px';
  }

  function hideTooltip() {
    if (tooltip) tooltip.dataset.show = '0';
  }
  /* ---- feature popup --------------------------------------------------- */

  /*
   * One popup serves every route into a feature: clicking an annotation,
   * clicking an element name in a table, or following a coordinate link. It is
   * a native <dialog> so it survives the pointer leaving the track -- the
   * reason a hover tooltip could never carry a working copy button.
   *
   * Metadata and sequence live in separate tabs rather than one long panel:
   * the sequence can be tens of kilobases and would bury the details that
   * identify the feature.
   */

  var lastFocus = null;

  function openFeature(h) {
    if (!modal || !h) return;
    var element = h.element;
    var sequence = element
      ? (data.sequences.seq || {})[element.pair_id] || null
      : null;

    lastFocus = document.activeElement;
    modal.innerHTML = '';

    modal.appendChild(featureHead(h, element));

    var body = document.createElement('div');
    body.className = 'modal-body';

    var panels = [{ id: 'details', label: 'Details', build: function () {
      return detailsPanel(h, element);
    } }];

    if (element) {
      panels.push({
        id: 'sequence',
        label: sequence ? 'Element sequence' : 'Element sequence (not embedded)',
        build: function () {
          return sequencePanel(h, element, sequence);
        },
      });
    }

    body.appendChild(tabs(panels));
    modal.appendChild(body);

    if (typeof modal.showModal === 'function') modal.showModal();
    else modal.setAttribute('open', '');
  }

  function featureHead(h, element) {
    var head = document.createElement('div');
    head.className = 'modal-head';

    var titles = document.createElement('div');
    var h3 = document.createElement('h3');
    var group = h.groups[0];
    if (group) {
      var swatch = document.createElement('span');
      swatch.className = 'swatch';
      swatch.style.background = T.groupColour(group);
      swatch.style.marginRight = '7px';
      h3.appendChild(swatch);
    }
    h3.appendChild(
      document.createTextNode(element ? element.element_id : h.model.name)
    );

    var sub = document.createElement('div');
    sub.className = 'sub';
    sub.textContent = element
      ? h.contig.name + ':' + element.start.toLocaleString() + '–' +
        element.end.toLocaleString() + ' · ' + T.formatBp(element.length)
      : h.contig.name + ':' + h.start.toLocaleString() + '–' +
        h.end.toLocaleString() + ' · ' + T.formatBp(h.length) + ' · unpaired';

    titles.appendChild(h3);
    titles.appendChild(sub);

    var actions = document.createElement('div');
    actions.className = 'modal-head-actions';

    var show = document.createElement('button');
    show.type = 'button';
    show.textContent = 'Show on track';
    show.addEventListener('click', function () {
      closeModal();
      if (element) {
        goToLocus(h.contig.name, element.start, element.end);
      } else {
        goToLocus(h.contig.name, h.start, h.end);
      }
    });

    var close = document.createElement('button');
    close.type = 'button';
    close.textContent = 'Close';
    close.addEventListener('click', closeModal);

    actions.appendChild(show);
    actions.appendChild(close);
    head.appendChild(titles);
    head.appendChild(actions);
    return head;
  }

  function closeModal() {
    if (!modal) return;
    if (typeof modal.close === 'function' && modal.open) modal.close();
    else modal.removeAttribute('open');
    if (lastFocus && lastFocus.focus) lastFocus.focus();
    lastFocus = null;
  }

  /* A minimal tab strip: no dependency, and the roles keep it navigable. */
  function tabs(panels) {
    var wrap = document.createElement('div');

    if (panels.length < 2) {
      wrap.appendChild(panels[0].build());
      return wrap;
    }

    var strip = document.createElement('div');
    strip.className = 'tabs';
    strip.setAttribute('role', 'tablist');

    var bodies = [];
    var buttons = [];

    panels.forEach(function (panel, index) {
      var button = document.createElement('button');
      button.type = 'button';
      button.className = 'tab';
      button.textContent = panel.label;
      button.setAttribute('role', 'tab');
      button.setAttribute('aria-selected', index === 0 ? 'true' : 'false');

      var content = document.createElement('div');
      content.className = 'tab-panel';
      content.setAttribute('role', 'tabpanel');
      // Panels are built lazily: rendering a 30 kb sequence into the DOM for a
      // tab the reader may never open is wasted work on every popup.
      if (index === 0) content.appendChild(panel.build());
      else content.hidden = true;

      button.addEventListener('click', function () {
        buttons.forEach(function (other, otherIndex) {
          other.setAttribute('aria-selected', otherIndex === index ? 'true' : 'false');
          bodies[otherIndex].hidden = otherIndex !== index;
        });
        if (!content.firstChild) content.appendChild(panel.build());
      });

      buttons.push(button);
      bodies.push(content);
      strip.appendChild(button);
    });

    wrap.appendChild(strip);
    bodies.forEach(function (content) {
      wrap.appendChild(content);
    });
    return wrap;
  }

  function detailsPanel(h, element) {
    var dl = document.createElement('dl');

    row(dl, 'Model', h.model.name);
    row(
      dl,
      'Hit location',
      h.contig.name + ':' + h.start.toLocaleString() + '–' + h.end.toLocaleString()
    );
    row(dl, 'Hit length', T.formatBp(h.length));
    row(dl, 'Strand', h.strand);
    row(dl, 'E-value', T.formatEvalue(h.evalue));
    row(dl, 'Score', h.score == null ? null : String(h.score));

    if (h.hmmStart != null && h.hmmEnd != null) {
      row(
        dl,
        'Model match',
        h.hmmStart + '–' + h.hmmEnd +
          (h.model.length ? ' of ' + h.model.length : '') +
          (h.modelCoverage != null ? ' (' + T.formatPercent(h.modelCoverage) + ')' : '')
      );
    } else {
      row(dl, 'Model match', 'coordinates unavailable for this input format');
    }
    if (h.spanCoverage != null) {
      row(dl, 'Span / model', h.spanCoverage.toFixed(2) + '×');
    }

    [T.truncationNote(h, 'left'), T.truncationNote(h, 'right')]
      .filter(Boolean)
      .forEach(function (note, index) {
        row(dl, index === 0 ? 'Completeness' : ' ', note);
      });

    row(
      dl,
      'Pairing group',
      h.groups.map(function (g) { return g.label; }).join(', ') || 'none'
    );

    if (element) {
      var left = T.hitByUid(element.left_uid);
      var right = T.hitByUid(element.right_uid);
      row(dl, 'Terminus', h.role || '—');
      if (left) {
        row(
          dl,
          'Left terminus',
          left.model.name + ' ' + left.start.toLocaleString() + '–' +
            left.end.toLocaleString() + ' (' + left.strand + ')'
        );
      }
      if (right) {
        row(
          dl,
          'Right terminus',
          right.model.name + ' ' + right.start.toLocaleString() + '–' +
            right.end.toLocaleString() + ' (' + right.strand + ')'
        );
      }
      row(dl, 'Between termini', T.formatBp(element.inner_distance));
      row(dl, 'Element length', T.formatBp(element.length));
    } else {
      row(dl, 'Pairing', 'unpaired');
    }

    return dl;
  }

  function sequencePanel(h, element, sequence) {
    var panel = document.createElement('div');

    if (!sequence) {
      var missing = document.createElement('div');
      missing.className = 'seq-missing';
      missing.textContent = data.sequences.embedded
        ? 'This element’s sequence was too large to embed in the report. ' +
          'Extract it from the run’s output:'
        : 'Sequences were not embedded in this report. Extract this element ' +
          'from the run’s output:';
      var code = document.createElement('code');
      code.textContent =
        'samtools faidx <genome.fa> ' + h.contig.name + ':' +
        element.start + '-' + element.end;
      missing.appendChild(code);
      panel.appendChild(missing);
      return panel;
    }

    var reversed = false;
    var actions = document.createElement('div');
    actions.className = 'modal-actions';

    var copy = document.createElement('button');
    copy.type = 'button';
    copy.textContent = 'Copy sequence';

    var fasta = document.createElement('button');
    fasta.type = 'button';
    fasta.textContent = 'Copy as FASTA';

    var download = document.createElement('button');
    download.type = 'button';
    download.textContent = 'Download FASTA';

    var flip = document.createElement('button');
    flip.type = 'button';
    flip.textContent = 'Reverse complement';

    var status = document.createElement('span');
    status.className = 'status';
    status.setAttribute('role', 'status');
    // Extraction always takes the plus strand, whichever strand the termini
    // were found on, so the panel says which strand is on screen.
    status.textContent = 'Plus strand, as extracted · ' + sequence.length + ' bp';

    actions.appendChild(copy);
    actions.appendChild(fasta);
    actions.appendChild(download);
    actions.appendChild(flip);
    actions.appendChild(status);
    panel.appendChild(actions);

    var pre = document.createElement('pre');
    pre.className = 'seq';
    pre.textContent = T.wrapSequence(sequence, 60);
    panel.appendChild(pre);

    function current() {
      return reversed ? T.reverseComplement(sequence) : sequence;
    }
    function header() {
      return (
        '>' + element.element_id + ' ' + h.contig.name + ':' +
        element.start + '-' + element.end +
        (reversed ? ' (reverse complement)' : '')
      );
    }

    flip.addEventListener('click', function () {
      reversed = !reversed;
      pre.textContent = T.wrapSequence(current(), 60);
      status.textContent =
        (reversed ? 'Reverse complement' : 'Plus strand, as extracted') +
        ' · ' + sequence.length + ' bp';
    });
    copy.addEventListener('click', function () {
      report(status, T.copyText(current()), 'Sequence copied.');
    });
    fasta.addEventListener('click', function () {
      report(
        status,
        T.copyText(header() + '\n' + T.wrapSequence(current(), 60) + '\n'),
        'FASTA copied.'
      );
    });
    download.addEventListener('click', function () {
      T.downloadText(
        header() + '\n' + T.wrapSequence(current(), 60) + '\n',
        element.element_id + '.fasta'
      );
    });

    return panel;
  }

  function report(status, promise, message) {
    var previous = status.textContent;
    promise.then(
      function () {
        status.textContent = message;
        setTimeout(function () {
          status.textContent = previous;
        }, 2200);
      },
      function () {
        status.textContent = 'Copy failed — select the text and copy manually.';
      }
    );
  }

  /* ---- navigation ------------------------------------------------------ */

  /*
   * Reveal a locus on its track and frame it. Called from the popup and from
   * the coordinate links in the tables, so a reader can move from a row to the
   * picture without hunting for the contig.
   */
  function goToLocus(contigName, start, end) {
    var host = container.querySelector(
      '.track[data-contig="' + cssEscape(contigName) + '"]'
    );

    if (!host) {
      // The track is filtered out. Clear the filters that hide it rather than
      // silently doing nothing.
      var search = document.getElementById('contig-search');
      var toggle = document.getElementById('show-empty');
      if (search && search.value) {
        search.value = '';
        search.dispatchEvent(new Event('input'));
      }
      if (toggle && !toggle.checked) {
        toggle.checked = true;
        toggle.dispatchEvent(new Event('change'));
      }
      host = container.querySelector(
        '.track[data-contig="' + cssEscape(contigName) + '"]'
      );
      if (!host) return;
    }

    var track = host.__track;
    if (!track) return;

    var pad = Math.max(200, Math.round((end - start) * 0.25));
    track.built = true;
    track.setView(start - pad, end + pad);
    track.draw();
    host.scrollIntoView({ block: 'center', behavior: 'smooth' });

    host.classList.add('track-flash');
    setTimeout(function () {
      host.classList.remove('track-flash');
    }, 1200);
  }

  // Exposed so the tables and the alignment panels can drive the same popup
  // and the same navigation rather than reimplementing either.
  T.showFeature = function (uid) {
    openFeature(T.hitByUid(uid));
  };
  T.showElement = function (pairId) {
    for (var i = 0; i < data.elements.length; i++) {
      if (data.elements[i].pair_id === pairId) {
        openFeature(T.hitByUid(data.elements[i].left_uid));
        return;
      }
    }
  };
  T.goToLocus = goToLocus;

  /* ---- assembly -------------------------------------------------------- */

  /*
   * A draft assembly can have tens of thousands of scaffolds. Building every
   * SVG up front would stall the page, so a track is drawn the first time it
   * scrolls into view.
   */
  var observer =
    'IntersectionObserver' in window
      ? new IntersectionObserver(
          function (entries) {
            entries.forEach(function (entry) {
              if (!entry.isIntersecting) return;
              var track = entry.target.__track;
              if (track && !track.built) {
                track.built = true;
                track.draw();
              }
            });
          },
          { rootMargin: '300px' }
        )
      : null;

  function build() {
    container.innerHTML = '';
    var shown = 0;
    tracks = [];

    data.contigs.forEach(function (contig) {
      if (!showEmpty && contig.n_pairs === 0) return;
      if (filterText && contig.name.toLowerCase().indexOf(filterText) === -1) return;
      var track = new Track(contig);
      track.el.__track = track;
      container.appendChild(track.el);
      tracks.push(track);
      shown++;
      if (observer) observer.observe(track.el);
      else {
        track.built = true;
        track.draw();
      }
    });

    if (!shown) {
      var empty = document.createElement('div');
      empty.className = 'empty-state';
      empty.textContent = filterText
        ? 'No sequences match “' + filterText + '”.'
        : 'No sequences carry a terminus pair. Enable “Show sequences without pairs” to see unpaired hits.';
      container.appendChild(empty);
    }

    var counter = document.getElementById('track-count');
    if (counter) {
      counter.textContent =
        shown.toLocaleString() + ' of ' + data.contigs.length.toLocaleString() + ' sequences';
    }
  }

  var search = document.getElementById('contig-search');
  if (search) {
    search.addEventListener('input', function () {
      filterText = search.value.trim().toLowerCase();
      build();
    });
  }

  var toggle = document.getElementById('show-empty');
  if (toggle) {
    showEmpty = toggle.checked;
    toggle.addEventListener('change', function () {
      showEmpty = toggle.checked;
      build();
    });
  }

  var resizeTimer = null;
  window.addEventListener('resize', function () {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(function () {
      tracks.forEach(function (track) {
        if (track.built) track.schedule();
      });
    }, 150);
  });

  if (window.matchMedia) {
    var scheme = window.matchMedia('(prefers-color-scheme: dark)');
    var onScheme = function () {
      tracks.forEach(function (track) {
        if (track.built) track.schedule();
      });
    };
    if (scheme.addEventListener) scheme.addEventListener('change', onScheme);
  }

  build();
})();

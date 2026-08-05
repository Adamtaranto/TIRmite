/*
 * Tables of every hit and every predicted element.
 *
 * These are built in the browser rather than in the HTML because a run can
 * retain hundreds of thousands of hits, and writing that many <tr> elements
 * into the document would make the file enormous and slow to open. The payload
 * already holds the data; the table renders a filtered window of it.
 *
 * Only a capped number of rows is put in the DOM at once. That cap is a
 * rendering limit, not a data limit: the filter searches every row, the count
 * always states the true total, and the download button writes the complete
 * table as TSV.
 */

(function () {
  'use strict';

  var T = window.TIRmite;
  if (!T) return;

  var data = T.data;
  var ROW_CAP = 500;

  function text(value) {
    return value == null ? '' : String(value);
  }

  function coords(contig, start, end) {
    return contig + ':' + start.toLocaleString() + '–' + end.toLocaleString();
  }

  /* Describes why a hit is shorter than its model, in one word per end. */
  function truncationSummary(h) {
    var parts = [];
    if (h.truncLeft) parts.push(h.clipLeft ? 'contig start' : 'partial left');
    if (h.truncRight) parts.push(h.clipRight ? 'contig end' : 'partial right');
    return parts.join(', ') || 'full length';
  }

  /*
   * A cell that navigates to a locus on its track. Rendered as a button rather
   * than an anchor because it moves the viewport rather than following a URL,
   * and a fake href would break opening the report in a new tab.
   */
  function locusCell(contigName, start, end, label) {
    var button = document.createElement('button');
    button.type = 'button';
    button.className = 'locus-link';
    button.textContent = label;
    button.title = 'Show ' + contigName + ':' + start + '–' + end + ' on its track';
    button.addEventListener('click', function () {
      if (T.goToLocus) T.goToLocus(contigName, start, end);
    });
    return button;
  }

  function nameCell(label, onOpen) {
    var button = document.createElement('button');
    button.type = 'button';
    button.className = 'name-link';
    button.textContent = label;
    button.title = 'Open details for ' + label;
    button.addEventListener('click', onOpen);
    return button;
  }

  var HIT_COLUMNS = [
    { label: 'Hit', num: true, get: function (h) { return h.uid; } },
    {
      label: 'Model',
      get: function (h) { return h.model.name; },
      render: function (h) {
        return nameCell(h.model.name, function () {
          T.showFeature(h.uid);
        });
      },
    },
    { label: 'Sequence', get: function (h) { return h.contig.name; } },
    {
      label: 'Start',
      num: true,
      get: function (h) { return h.start; },
      render: function (h) {
        return locusCell(h.contig.name, h.start, h.end, h.start.toLocaleString());
      },
    },
    {
      label: 'End',
      num: true,
      get: function (h) { return h.end; },
      render: function (h) {
        return locusCell(h.contig.name, h.start, h.end, h.end.toLocaleString());
      },
    },
    { label: 'Length', num: true, get: function (h) { return h.length; } },
    { label: 'Strand', get: function (h) { return h.strand; } },
    {
      label: 'E-value',
      num: true,
      get: function (h) { return h.evalue; },
      show: function (h) { return T.formatEvalue(h.evalue); },
    },
    { label: 'Score', num: true, get: function (h) { return h.score; } },
    {
      label: 'Model match',
      num: true,
      get: function (h) { return h.modelCoverage; },
      show: function (h) { return T.formatPercent(h.modelCoverage); },
    },
    { label: 'Completeness', get: truncationSummary },
    {
      label: 'Pairing group',
      get: function (h) {
        return h.groups.map(function (g) { return g.label; }).join(', ') || '—';
      },
    },
    { label: 'Terminus', get: function (h) { return h.role || '—'; } },
    {
      label: 'Element',
      get: function (h) { return h.element ? h.element.element_id : '—'; },
      render: function (h) {
        if (!h.element) return null;
        return nameCell(h.element.element_id, function () {
          T.showElement(h.element.pair_id);
        });
      },
    },
  ];

  var ELEMENT_COLUMNS = [
    {
      label: 'Element',
      get: function (e) { return e.element_id; },
      render: function (e) {
        return nameCell(e.element_id, function () {
          T.showElement(e.pair_id);
        });
      },
    },
    {
      label: 'Pairing group',
      get: function (e) { return data.groups[e.group_i].label; },
    },
    {
      label: 'Sequence',
      get: function (e) { return data.contigs[e.contig_i].name; },
    },
    {
      label: 'Start',
      num: true,
      get: function (e) { return e.start; },
      render: function (e) {
        return locusCell(
          data.contigs[e.contig_i].name, e.start, e.end, e.start.toLocaleString()
        );
      },
    },
    {
      label: 'End',
      num: true,
      get: function (e) { return e.end; },
      render: function (e) {
        return locusCell(
          data.contigs[e.contig_i].name, e.start, e.end, e.end.toLocaleString()
        );
      },
    },
    { label: 'Length', num: true, get: function (e) { return e.length; } },
    {
      label: 'Left terminus',
      get: function (e) {
        var h = T.hitByUid(e.left_uid);
        return h ? h.model.name + ' ' + coords('', h.start, h.end).slice(1) : '';
      },
    },
    {
      label: 'Right terminus',
      get: function (e) {
        var h = T.hitByUid(e.right_uid);
        return h ? h.model.name + ' ' + coords('', h.start, h.end).slice(1) : '';
      },
    },
    {
      label: 'Between termini',
      num: true,
      get: function (e) { return e.inner_distance; },
    },
    {
      label: 'Sequence embedded',
      get: function (e) {
        return (data.sequences.seq || {})[e.pair_id] ? 'yes' : 'no';
      },
    },
  ];

  function cellText(column, row) {
    if (column.show) return column.show(row);
    var value = column.get(row);
    if (value == null) return '–';
    return typeof value === 'number' ? value.toLocaleString() : String(value);
  }

  function DataTable(container, columns, rows, options) {
    this.columns = columns;
    this.rows = rows;
    this.filtered = rows;
    this.sortIndex = null;
    this.ascending = true;
    this.options = options || {};

    var controls = document.createElement('div');
    controls.className = 'controls';

    var label = document.createElement('label');
    label.appendChild(document.createTextNode('Filter '));
    var search = document.createElement('input');
    search.type = 'search';
    search.placeholder = options.placeholder || 'any value';
    search.autocomplete = 'off';
    label.appendChild(search);
    controls.appendChild(label);

    var count = document.createElement('span');
    count.className = 'muted';
    controls.appendChild(count);

    var download = document.createElement('button');
    download.type = 'button';
    download.textContent = 'Download all as TSV';
    controls.appendChild(download);

    var wrap = document.createElement('div');
    wrap.className = 'table-wrap';
    var table = document.createElement('table');
    var thead = document.createElement('thead');
    var headRow = document.createElement('tr');
    var self = this;

    columns.forEach(function (column, index) {
      var th = document.createElement('th');
      th.scope = 'col';
      th.textContent = column.label;
      if (column.num) th.className = 'num';
      th.setAttribute('role', 'button');
      th.setAttribute('tabindex', '0');
      function activate() {
        self.sortBy(index);
      }
      th.addEventListener('click', activate);
      th.addEventListener('keydown', function (event) {
        if (event.key === 'Enter' || event.key === ' ') {
          event.preventDefault();
          activate();
        }
      });
      headRow.appendChild(th);
    });

    thead.appendChild(headRow);
    table.appendChild(thead);
    var tbody = document.createElement('tbody');
    table.appendChild(tbody);
    wrap.appendChild(table);

    container.appendChild(controls);
    container.appendChild(wrap);

    this.headers = headRow.cells;
    this.tbody = tbody;
    this.count = count;

    search.addEventListener('input', function () {
      self.filter(search.value.trim().toLowerCase());
    });
    download.addEventListener('click', function () {
      self.download();
    });

    this.render();
  }

  DataTable.prototype.filter = function (query) {
    var columns = this.columns;
    if (!query) {
      this.filtered = this.rows;
    } else {
      this.filtered = this.rows.filter(function (row) {
        for (var i = 0; i < columns.length; i++) {
          if (cellText(columns[i], row).toLowerCase().indexOf(query) !== -1) {
            return true;
          }
        }
        return false;
      });
    }
    this.render();
  };

  DataTable.prototype.sortBy = function (index) {
    this.ascending = this.sortIndex === index ? !this.ascending : true;
    this.sortIndex = index;

    var column = this.columns[index];
    var direction = this.ascending ? 1 : -1;
    this.filtered = this.filtered.slice().sort(function (a, b) {
      var av = column.get(a);
      var bv = column.get(b);
      // Missing values sort last whichever way the column is pointing.
      if (av == null && bv == null) return 0;
      if (av == null) return 1;
      if (bv == null) return -1;
      if (typeof av === 'string') {
        return av.toLowerCase() < bv.toLowerCase() ? -direction : direction;
      }
      return av === bv ? 0 : (av < bv ? -direction : direction);
    });

    for (var i = 0; i < this.headers.length; i++) {
      this.headers[i].removeAttribute('aria-sort');
    }
    this.headers[index].setAttribute(
      'aria-sort',
      this.ascending ? 'ascending' : 'descending'
    );
    this.render();
  };

  DataTable.prototype.render = function () {
    var columns = this.columns;
    var shown = this.filtered.slice(0, ROW_CAP);

    // Built off-document: appending row by row to a live tbody forces a
    // reflow per row, which is visible at a few hundred rows.
    var fragment = document.createDocumentFragment();
    shown.forEach(function (row) {
      var tr = document.createElement('tr');
      columns.forEach(function (column) {
        var td = document.createElement('td');
        if (column.num) td.className = 'num';
        // A column may render an interactive cell -- a coordinate that jumps
        // to the track, a name that opens the popup -- and falls back to text
        // when there is nothing to link to.
        var node = column.render ? column.render(row) : null;
        if (node) td.appendChild(node);
        else td.textContent = cellText(column, row);
        tr.appendChild(td);
      });
      fragment.appendChild(tr);
    });

    this.tbody.replaceChildren(fragment);

    var total = this.rows.length;
    var matched = this.filtered.length;
    var parts = [];
    if (matched !== total) {
      parts.push(matched.toLocaleString() + ' of ' + total.toLocaleString() + ' match');
    } else {
      parts.push(total.toLocaleString() + ' rows');
    }
    if (matched > ROW_CAP) {
      parts.push('showing the first ' + ROW_CAP.toLocaleString());
    }
    this.count.textContent = parts.join(' · ');
  };

  DataTable.prototype.download = function () {
    var columns = this.columns;
    var lines = [columns.map(function (c) { return c.label; }).join('\t')];
    this.filtered.forEach(function (row) {
      lines.push(
        columns
          .map(function (column) {
            var value = column.get(row);
            return value == null ? 'NA' : text(value);
          })
          .join('\t')
      );
    });

    // A blob URL keeps this working from file://, where the report normally
    // lives; nothing is sent anywhere.
    var blob = new Blob([lines.join('\n') + '\n'], {
      type: 'text/tab-separated-values',
    });
    var url = URL.createObjectURL(blob);
    var link = document.createElement('a');
    link.href = url;
    link.download = this.options.filename || 'tirmite_table.tsv';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    setTimeout(function () {
      URL.revokeObjectURL(url);
    }, 0);
  };

  var hitContainer = document.getElementById('all-hits-table');
  if (hitContainer) {
    var hits = [];
    for (var i = 0; i < data.hits.n; i++) hits.push(T.hit(i));
    new DataTable(hitContainer, HIT_COLUMNS, hits, {
      placeholder: 'model, sequence, group…',
      filename: 'tirmite_hits.tsv',
    });
  }

  var elementContainer = document.getElementById('elements-table');
  if (elementContainer && data.elements.length) {
    new DataTable(elementContainer, ELEMENT_COLUMNS, data.elements.slice(), {
      placeholder: 'element, sequence, group…',
      filename: 'tirmite_elements.tsv',
    });
  }
})();

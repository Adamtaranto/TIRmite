/*
 * Sortable statistics tables.
 *
 * The tables are the report's table view: every number a glyph or tooltip
 * shows is also readable here, so nothing is gated behind hover.
 */

(function () {
  'use strict';

  function cellValue(row, index) {
    var cell = row.cells[index];
    if (!cell) return '';
    // Cells carry data-sort when their text is formatted for reading
    // ("1.2 Mb", "12%") and would sort wrongly as a string.
    var raw = cell.dataset.sort;
    if (raw !== undefined) {
      var n = parseFloat(raw);
      return isNaN(n) ? raw : n;
    }
    var text = cell.textContent.trim();
    var num = parseFloat(text.replace(/,/g, ''));
    return isNaN(num) ? text.toLowerCase() : num;
  }

  function sortTable(table, index, ascending) {
    var body = table.tBodies[0];
    if (!body) return;
    var rows = Array.prototype.slice.call(body.rows);
    rows.sort(function (a, b) {
      var av = cellValue(a, index);
      var bv = cellValue(b, index);
      if (av === bv) return 0;
      if (av === '' || av == null) return 1;
      if (bv === '' || bv == null) return -1;
      return (av < bv ? -1 : 1) * (ascending ? 1 : -1);
    });
    rows.forEach(function (row) {
      body.appendChild(row);
    });
  }

  document.querySelectorAll('table[data-sortable]').forEach(function (table) {
    var headers = table.tHead ? table.tHead.rows[0].cells : [];
    Array.prototype.forEach.call(headers, function (th, index) {
      th.setAttribute('role', 'button');
      th.setAttribute('tabindex', '0');
      th.title = 'Sort by ' + th.textContent.trim();

      function activate() {
        var ascending = th.getAttribute('aria-sort') !== 'ascending';
        Array.prototype.forEach.call(headers, function (other) {
          other.removeAttribute('aria-sort');
        });
        th.setAttribute('aria-sort', ascending ? 'ascending' : 'descending');
        sortTable(table, index, ascending);
      }

      th.addEventListener('click', activate);
      th.addEventListener('keydown', function (event) {
        if (event.key === 'Enter' || event.key === ' ') {
          event.preventDefault();
          activate();
        }
      });
    });
  });
})();

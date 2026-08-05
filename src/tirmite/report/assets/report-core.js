/*
 * Shared runtime for the TIRmite report.
 *
 * Parses the embedded payload, rehydrates the columnar hit arrays into a
 * per-hit view, and exposes the small helpers the other bundles use. No
 * dependencies: the report has to work from a file:// URL with no network.
 */

(function () {
  'use strict';

  var SUPPORTED_SCHEMA = 1;

  var node = document.getElementById('tirmite-report-data');
  var data;
  try {
    data = JSON.parse(node.textContent);
  } catch (err) {
    showFatal('This report could not be read: its embedded data is corrupt.');
    return;
  }

  if (data.schema_version !== SUPPORTED_SCHEMA) {
    // Better a clear message than a blank page.
    showFatal(
      'This report was written by a different version of TIRmite (data format ' +
        data.schema_version +
        ', this viewer understands ' +
        SUPPORTED_SCHEMA +
        '). Regenerate it with the current version.'
    );
    return;
  }

  function showFatal(message) {
    var box = document.createElement('div');
    box.className = 'warnings';
    box.innerHTML = '<h2>Report unavailable</h2>';
    var p = document.createElement('p');
    p.textContent = message;
    box.appendChild(p);
    var wrap = document.querySelector('.wrap') || document.body;
    wrap.insertBefore(box, wrap.firstChild);
  }

  var ROLE_NAMES = { 0: null, 1: 'left', 2: 'right' };

  /*
   * Hits are stored as parallel arrays because that is far smaller and faster
   * to parse than a list of objects. Materialising one object per hit up front
   * would give back everything the columnar form saved, so a hit object is
   * built only when something actually asks for one.
   */
  function hit(i) {
    var cols = data.hits;
    if (i == null || i < 0 || i >= cols.n) return null;
    return {
      i: i,
      uid: cols.uid[i],
      model: data.models[cols.model_i[i]],
      contig: data.contigs[cols.contig_i[i]],
      start: cols.start[i],
      end: cols.end[i],
      length: cols.end[i] - cols.start[i] + 1,
      strand: cols.strand.charAt(i),
      evalue: cols.evalue[i],
      score: cols.score[i],
      hmmStart: cols.hmm_start[i],
      hmmEnd: cols.hmm_end[i],
      modelCoverage: cols.model_cov[i],
      spanCoverage: cols.span_cov[i],
      truncLeft: cols.trunc_l[i],
      truncRight: cols.trunc_r[i],
      clipLeft: cols.clip_l[i] === 1,
      clipRight: cols.clip_r[i] === 1,
      row: cols.row[i],
      overflow: cols.overflow[i] === 1,
      groups: (cols.group_ix[i] || []).map(function (g) {
        return data.groups[g];
      }),
      role: ROLE_NAMES[cols.role[i]] || null,
      element: cols.pair_ix[i] == null ? null : data.elements[cols.pair_ix[i]],
    };
  }

  /* Index from uid to column position, built once and only if needed. */
  var uidIndex = null;
  function hitByUid(uid) {
    if (uidIndex === null) {
      uidIndex = new Map();
      for (var i = 0; i < data.hits.n; i++) uidIndex.set(data.hits.uid[i], i);
    }
    var i = uidIndex.get(uid);
    return i === undefined ? null : hit(i);
  }

  /*
   * First hit in [lo, hi) whose end reaches `bp`. Hits are sorted by
   * (contig, start), which is not the same order as by end, so the search
   * finds the first hit *starting* at or after bp and the caller walks back a
   * bounded distance to catch long hits that started earlier.
   */
  function firstStartingAt(lo, hi, bp) {
    var starts = data.hits.start;
    while (lo < hi) {
      var mid = (lo + hi) >> 1;
      if (starts[mid] < bp) lo = mid + 1;
      else hi = mid;
    }
    return lo;
  }

  function formatBp(n) {
    if (n == null) return '–';
    if (n >= 1e6) return (n / 1e6).toFixed(n >= 1e7 ? 0 : 2) + ' Mb';
    if (n >= 1e4) return (n / 1e3).toFixed(n >= 1e5 ? 0 : 1) + ' kb';
    return n.toLocaleString() + ' bp';
  }

  function formatInt(n) {
    return n == null ? '–' : n.toLocaleString();
  }

  function formatEvalue(v) {
    if (v == null) return '–';
    if (v === 0) return '0';
    return v < 0.001 || v >= 1000 ? v.toExponential(1) : String(+v.toFixed(4));
  }

  function formatPercent(v) {
    return v == null ? '–' : Math.round(v * 100) + '%';
  }

  /*
   * Describes why a hit does not span its whole model. The distinction matters:
   * a partial model match and a hit that ran off the end of the contig are
   * different findings and the report draws them differently.
   */
  function truncationNote(h, side) {
    var bp = side === 'left' ? h.truncLeft : h.truncRight;
    var clipped = side === 'left' ? h.clipLeft : h.clipRight;
    if (bp == null) return null;
    if (bp === 0) return null;
    var where = side === 'left' ? 'lower' : 'higher';
    if (clipped) {
      return (
        formatInt(bp) +
        ' bp of model unmatched at the ' +
        where +
        '-coordinate end — truncated by the contig end.'
      );
    }
    return (
      formatInt(bp) +
      ' bp of model unmatched at the ' +
      where +
      '-coordinate end — incomplete match.'
    );
  }

  function reverseComplement(seq) {
    var map = {
      A: 'T', T: 'A', G: 'C', C: 'G', U: 'A', N: 'N',
      a: 't', t: 'a', g: 'c', c: 'g', u: 'a', n: 'n',
      R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K',
      B: 'V', V: 'B', D: 'H', H: 'D', '-': '-', '.': '.',
    };
    var out = new Array(seq.length);
    for (var i = 0, j = seq.length - 1; j >= 0; i++, j--) {
      var c = seq.charAt(j);
      out[i] = map[c] || c;
    }
    return out.join('');
  }

  function wrapSequence(seq, width) {
    var lines = [];
    for (var i = 0; i < seq.length; i += width) {
      lines.push(seq.slice(i, i + width));
    }
    return lines.join('\n');
  }

  /*
   * Reports are usually opened over file://, where the async Clipboard API is
   * unavailable in several browsers. The execCommand fallback is therefore the
   * common path here rather than an edge case.
   */
  function copyText(text) {
    if (navigator.clipboard && window.isSecureContext) {
      return navigator.clipboard.writeText(text);
    }
    return new Promise(function (resolve, reject) {
      var area = document.createElement('textarea');
      area.value = text;
      area.setAttribute('readonly', '');
      area.style.position = 'fixed';
      area.style.opacity = '0';
      document.body.appendChild(area);
      area.select();
      var ok = false;
      try {
        ok = document.execCommand('copy');
      } catch (err) {
        ok = false;
      }
      document.body.removeChild(area);
      ok ? resolve() : reject(new Error('copy failed'));
    });
  }

  function isDark() {
    return (
      window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches
    );
  }

  /* A group's colour depends on the surface it is drawn against. */
  function groupColour(group) {
    if (!group) return null;
    return (isDark() && group.colour_dark) || group.colour;
  }

  window.TIRmite = {
    data: data,
    hit: hit,
    hitByUid: hitByUid,
    firstStartingAt: firstStartingAt,
    formatBp: formatBp,
    formatInt: formatInt,
    formatEvalue: formatEvalue,
    formatPercent: formatPercent,
    truncationNote: truncationNote,
    reverseComplement: reverseComplement,
    wrapSequence: wrapSequence,
    copyText: copyText,
    groupColour: groupColour,
    isDark: isDark,
  };
})();

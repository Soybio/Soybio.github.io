/* Tools page: forms, tables and plots. Calculations are in seqtools.js (validated against
   Biopython and SciPy by tools/validate_tools.py). */
(function () {
  'use strict';

  var S = window.SoybioSeq;
  var D = window.SoybioSeqData;
  var X = window.SoybioExport;
  var NS = 'http://www.w3.org/2000/svg';
  var FONT = 'Segoe UI, Roboto, Helvetica Neue, Arial, sans-serif';
  var GREEN = '#22573a', BLUE = '#2f6f9f', ORANGE = '#b7791f', RED = '#b83b3b', GREY = '#9aa3a0';
  var MAX_ALIGN_CELLS = 4e6;

  var $ = function (id) { return document.getElementById(id); };
  var last = {};

  // ---------- Small helpers ----------

  function el(tag, attrs, parent) {
    var e = document.createElementNS(NS, tag);
    for (var k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function txt(parent, x, y, s, attrs) {
    var a = { x: x, y: y };
    for (var k in attrs) a[k] = attrs[k];
    var t = el('text', a, parent);
    t.textContent = s;
    return t;
  }
  function fmt(v, d) {
    if (Math.abs(v) < Math.pow(10, -d) / 2) v = 0; // avoid "-0.000"
    return v.toLocaleString('en-US', { minimumFractionDigits: d, maximumFractionDigits: d });
  }
  function int(v) { return v.toLocaleString('en-US'); }
  function num(id) { return parseFloat($(id).value); }
  function show(id, html) { var b = $(id); b.innerHTML = html; b.hidden = false; }
  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]; }); }
  function error(id, msg) { show(id, '<p class="status error">' + esc(msg) + '</p>'); }
  function facts(rows) {
    return '<dl class="facts">' + rows.map(function (r) { return '<dt>' + esc(r[0]) + '</dt><dd>' + r[1] + '</dd>'; }).join('') + '</dl>';
  }
  function table(head, rows, numCols) {
    numCols = numCols || [];
    return '<div class="table-wrap"><table><thead><tr>' + head.map(function (h, i) {
      return '<th scope="col"' + (numCols.indexOf(i) >= 0 ? ' class="num"' : '') + '>' + esc(h) + '</th>';
    }).join('') + '</tr></thead><tbody>' + rows.map(function (r) {
      return '<tr>' + r.map(function (c, i) { return '<td' + (numCols.indexOf(i) >= 0 ? ' class="num"' : '') + '>' + c + '</td>'; }).join('') + '</tr>';
    }).join('') + '</tbody></table></div>';
  }
  function wrapSeq(s, w) {
    var out = [];
    for (var i = 0; i < s.length; i += w) out.push(s.slice(i, i + w));
    return out.join('\n');
  }
  function fasta(name, s) { return '>' + name + '\n' + wrapSeq(s, 60) + '\n'; }
  function saveText(text, name) { X.saveBlob(new Blob([text], { type: 'text/plain;charset=utf-8' }), name); }

  // Reads a sequence field and reports problems. Returns null on error.
  function getSeq(fieldId, type, outId, label) {
    var r = S.readSeq($(fieldId).value, type);
    if (!r.seq) { error(outId, 'Please enter ' + label + '.'); return null; }
    if (r.invalid.length) {
      error(outId, label.charAt(0).toUpperCase() + label.slice(1) + ' contains characters that are not allowed: ' + r.invalid.join(', ') +
        (type === 'dna' ? '. Use A, C, G, T (U is read as T) or N.' : '. Use the 20 standard amino acid letters.'));
      return null;
    }
    return r;
  }

  // ---------- Plots ----------

  function niceStep(range, count) {
    var raw = range / count;
    var mag = Math.pow(10, Math.floor(Math.log10(raw)));
    return [1, 2, 2.5, 5, 10].map(function (m) { return m * mag; }).find(function (s) { return s >= raw; });
  }

  function svgBox(box, W, H, label) {
    box.innerHTML = '';
    box.hidden = false;
    var svg = el('svg', { xmlns: NS, width: W, height: H, viewBox: '0 0 ' + W + ' ' + H, 'font-family': FONT, role: 'img', 'aria-label': label });
    el('rect', { x: 0, y: 0, width: W, height: H, fill: '#ffffff' }, svg);
    box.appendChild(svg);
    return svg;
  }

  // Axis frame with y ticks; returns scale functions.
  function frame(svg, o) {
    var yStep = niceStep(o.ymax - o.ymin || 1, 5);
    var ymin = Math.floor(o.ymin / yStep) * yStep, ymax = Math.ceil(o.ymax / yStep) * yStep;
    if (ymax === ymin) ymax = ymin + yStep;
    var Y = function (v) { return o.top + o.h - (v - ymin) / (ymax - ymin) * o.h; };
    for (var v = ymin; v <= ymax + yStep / 2; v += yStep) {
      var vv = +v.toPrecision(10);
      el('line', { x1: o.left, x2: o.left + o.w, y1: Y(vv), y2: Y(vv), stroke: vv === 0 ? '#8c959f' : '#e3e8e5' }, svg);
      txt(svg, o.left - 6, Y(vv) + 4, String(vv), { 'font-size': 11, fill: '#57606a', 'text-anchor': 'end' });
    }
    txt(svg, 14, o.top + o.h / 2, o.yLabel, { 'font-size': 12, fill: '#1f2328', 'text-anchor': 'middle', 'font-weight': 600,
      transform: 'rotate(-90 14 ' + (o.top + o.h / 2) + ')' });
    if (o.title) txt(svg, o.left, 22, o.title, { 'font-size': 14, 'font-weight': 700, fill: '#1f2328' });
    return Y;
  }

  function xAxis(svg, o, xmin, xmax, label) {
    var step = niceStep(xmax - xmin || 1, 8);
    var X_ = function (v) { return o.left + (v - xmin) / (xmax - xmin || 1) * o.w; };
    for (var v = Math.ceil(xmin / step) * step; v <= xmax; v += step) {
      var vv = +v.toPrecision(10);
      el('line', { x1: X_(vv), x2: X_(vv), y1: o.top + o.h, y2: o.top + o.h + 4, stroke: '#8c959f' }, svg);
      txt(svg, X_(vv), o.top + o.h + 17, String(vv), { 'font-size': 11, fill: '#57606a', 'text-anchor': 'middle' });
    }
    txt(svg, o.left + o.w / 2, o.top + o.h + 36, label, { 'font-size': 12, fill: '#1f2328', 'text-anchor': 'middle', 'font-weight': 600 });
    return X_;
  }

  function linePlot(box, o) {
    var W = Math.max(600, box.clientWidth || 800), H = 300;
    var svg = svgBox(box, W, H, o.title);
    var f = { left: 56, top: 40, w: W - 76, h: H - 90, yLabel: o.yLabel, title: o.title };
    var ys = o.points.map(function (p) { return p.y; });
    var ymin = o.ymin !== undefined ? o.ymin : Math.min.apply(null, ys);
    var ymax = o.ymax !== undefined ? o.ymax : Math.max.apply(null, ys);
    if (o.hline !== undefined) { ymin = Math.min(ymin, o.hline); ymax = Math.max(ymax, o.hline); }
    f.ymin = ymin; f.ymax = ymax;
    var Y = frame(svg, f);
    var Xs = xAxis(svg, f, o.xmin, o.xmax, o.xLabel);
    if (o.hline !== undefined) {
      el('line', { x1: f.left, x2: f.left + f.w, y1: Y(o.hline), y2: Y(o.hline), stroke: RED, 'stroke-dasharray': '5 4' }, svg);
      txt(svg, f.left + f.w, Y(o.hline) - 5, o.hlineLabel, { 'font-size': 11, fill: RED, 'text-anchor': 'end' });
    }
    el('polyline', { points: o.points.map(function (p) { return Xs(p.x) + ',' + Y(p.y); }).join(' '), fill: 'none', stroke: o.color || GREEN, 'stroke-width': 1.6 }, svg);
  }

  function barPlot(box, o) {
    var n = o.labels.length;
    var W = Math.max(600, box.clientWidth || 800, n * (o.minBar || 14) + 100), H = o.height || 300;
    var svg = svgBox(box, W, H, o.title);
    var bottom = o.rotate ? 60 : 50;
    var f = { left: 60, top: 40, w: W - 80, h: H - 40 - bottom, yLabel: o.yLabel, title: o.title };
    var hi = 0;
    o.values.forEach(function (v, i) {
      hi = Math.max(hi, v + (o.errors ? o.errors[i] : 0));
      if (o.points) o.points[i].forEach(function (p) { hi = Math.max(hi, p); });
    });
    f.ymin = 0; f.ymax = hi || 1;
    var Y = frame(svg, f);
    var bw = f.w / n;
    o.values.forEach(function (v, i) {
      var x = f.left + i * bw;
      var c = o.colors ? o.colors[i] : GREEN;
      el('rect', { x: x + bw * 0.15, y: Y(v), width: bw * 0.7, height: Math.max(0, Y(0) - Y(v)), fill: c, 'fill-opacity': 0.85 }, svg);
      if (o.errors && o.errors[i] > 0) {
        var cx = x + bw / 2;
        el('line', { x1: cx, x2: cx, y1: Y(v + o.errors[i]), y2: Y(Math.max(0, v - o.errors[i])), stroke: '#2b2f33' }, svg);
        el('line', { x1: cx - 5, x2: cx + 5, y1: Y(v + o.errors[i]), y2: Y(v + o.errors[i]), stroke: '#2b2f33' }, svg);
      }
      if (o.points) o.points[i].forEach(function (p, k) {
        el('circle', { cx: x + bw / 2 + (k - (o.points[i].length - 1) / 2) * 5, cy: Y(p), r: 2.6, fill: '#1f2328', 'fill-opacity': 0.65 }, svg);
      });
      var lx = x + bw / 2, ly = f.top + f.h + 14;
      txt(svg, lx, ly, o.labels[i], o.rotate
        ? { 'font-size': 10, fill: '#3d454d', 'text-anchor': 'end', transform: 'rotate(-90 ' + (lx + 3) + ' ' + (ly - 8) + ')' }
        : { 'font-size': 11.5, fill: '#3d454d', 'text-anchor': 'middle' });
    });
    if (o.hline !== undefined) {
      el('line', { x1: f.left, x2: f.left + f.w, y1: Y(o.hline), y2: Y(o.hline), stroke: RED, 'stroke-dasharray': '5 4' }, svg);
    }
  }

  function plotButtons(boxId, name) {
    return '<div class="btn-row plot-dl"><button type="button" class="btn btn-outline btn-small" data-svg="' + boxId + '" data-name="' + name +
      '">SVG</button><button type="button" class="btn btn-outline btn-small" data-png="' + boxId + '" data-name="' + name + '">PNG</button></div>';
  }

  // ---------- 1. Protein properties ----------

  function runProtein() {
    var r = getSeq('prot-seq', 'protein', 'prot-out', 'a protein sequence');
    if (!r) return;
    var p = r.seq;
    if (p.length < 2) { error('prot-out', 'Please enter at least 2 amino acids.'); return; }
    var q = S.proteinParams(p);
    var win = parseInt($('prot-win').value, 10);
    var html = (r.name ? '<p class="muted">' + esc(r.name) + '</p>' : '') + facts([
      ['Length', int(q.length) + ' amino acids'],
      ['Molecular weight', fmt(q.mw, 2) + ' Da (' + fmt(q.mw / 1000, 2) + ' kDa)'],
      ['Theoretical pI', fmt(q.pI, 2)],
      ['Net charge at pH 7.0', fmt(q.charge7, 2)],
      ['Negatively charged (Asp + Glu)', int(q.negative)],
      ['Positively charged (Arg + Lys)', int(q.positive)],
      ['Extinction coefficient (280 nm)', int(q.extCystines) + ' M⁻¹cm⁻¹ (all Cys as cystines); ' + int(q.extReduced) + ' M⁻¹cm⁻¹ (all Cys reduced)'],
      ['Abs 0.1% (1 g/L)', fmt(q.abs01Cystines, 3) + ' (cystines); ' + fmt(q.abs01Reduced, 3) + ' (reduced)'],
      ['Instability index', fmt(q.instability, 2) + (q.instability > 40 ? ' (classified as unstable)' : ' (classified as stable)')],
      ['Aliphatic index', fmt(q.aliphatic, 2)],
      ['GRAVY', fmt(q.gravy, 3) + (q.gravy > 0 ? ' (hydrophobic on average)' : ' (hydrophilic on average)')],
      ['Aromaticity', fmt(q.aromaticity, 3)]
    ]);
    if (!q.counts.W && !q.counts.Y) html += '<p class="note small">This protein has no Trp or Tyr, so its absorbance at 280 nm is too low to estimate concentration reliably.</p>';
    var aa = 'ACDEFGHIKLMNPQRSTVWY'.split('');
    html += '<h4>Amino acid composition</h4>' + plotButtons('prot-comp', 'amino_acid_composition') + '<div class="plot-box" id="prot-comp"></div>';
    html += table(['Amino acid', 'Count', '%'], aa.map(function (a) { return [a, int(q.counts[a]), fmt(q.counts[a] / q.length * 100, 1)]; }), [1, 2]);
    html += '<h4>Hydropathy plot (Kyte &amp; Doolittle, window ' + win + ')</h4>' + plotButtons('prot-hyd', 'hydropathy_plot') + '<div class="plot-box" id="prot-hyd"></div>';
    show('prot-out', html);
    barPlot($('prot-comp'), { labels: aa, values: aa.map(function (a) { return q.counts[a] / q.length * 100; }), yLabel: 'Composition (%)', title: 'Amino acid composition' });
    var hp = S.hydropathy(p, win);
    if (hp.length) {
      linePlot($('prot-hyd'), { points: hp, xmin: 1, xmax: p.length, xLabel: 'Position (amino acid)', yLabel: 'Hydropathy score',
        title: 'Kyte & Doolittle hydropathy', hline: win >= 19 ? 1.6 : 0, hlineLabel: win >= 19 ? 'Transmembrane threshold (1.6)' : '' });
    } else {
      $('prot-hyd').innerHTML = '<p class="muted small" style="padding:0.8rem">The sequence is shorter than the window.</p>';
    }
  }

  // ---------- 2. DNA statistics ----------

  function runDna() {
    var r = getSeq('dna-seq', 'dna', 'dna-out', 'a DNA sequence');
    if (!r) return;
    var s = r.seq, c = S.baseCounts(s), L = s.length;
    var win = parseInt($('dna-win').value, 10), step = Math.max(1, Math.round(win / 10));
    var rc = S.revcomp(s);
    last.dnaRc = { name: (r.name || 'sequence') + '_reverse_complement', seq: rc };
    var html = facts([
      ['Length', int(L) + ' bp'],
      ['GC content', fmt(S.gcContent(s), 2) + '%' + (c.N ? ' (N bases not counted)' : '')],
      ['A', int(c.A) + ' (' + fmt(c.A / L * 100, 1) + '%)'],
      ['C', int(c.C) + ' (' + fmt(c.C / L * 100, 1) + '%)'],
      ['G', int(c.G) + ' (' + fmt(c.G / L * 100, 1) + '%)'],
      ['T', int(c.T) + ' (' + fmt(c.T / L * 100, 1) + '%)'],
      ['N', int(c.N)]
    ]);
    html += '<h4>Reverse complement</h4><pre class="seq-out">' + esc(wrapSeq(rc, 60)) + '</pre>' +
      '<p><button type="button" class="btn btn-outline btn-small" id="dna-rc-dl">Download FASTA</button></p>';
    html += '<h4>GC content along the sequence (window ' + win + ' bp, step ' + step + ' bp)</h4>';
    if (L >= win) html += plotButtons('dna-gc', 'gc_content') + '<div class="plot-box" id="dna-gc"></div>';
    else html += '<p class="muted small">The sequence is shorter than the window.</p>';
    show('dna-out', html);
    $('dna-rc-dl').addEventListener('click', function () { saveText(fasta(last.dnaRc.name, last.dnaRc.seq), last.dnaRc.name + '.fasta'); });
    if (L >= win) {
      linePlot($('dna-gc'), { points: S.gcWindows(s, win, step), xmin: 1, xmax: L, ymin: 0, ymax: 100, xLabel: 'Position (bp)',
        yLabel: 'GC (%)', title: 'GC content', hline: S.gcContent(s), hlineLabel: 'Whole sequence' });
    }
  }

  // ---------- 3. ORF finder and translation ----------

  function runOrf() {
    var r = getSeq('orf-seq', 'dna', 'orf-out', 'a DNA sequence');
    if (!r) return;
    var s = r.seq, minAA = parseInt($('orf-min').value, 10);
    var orfs = S.findOrfs(s, minAA);
    last.orfs = { name: r.name || 'sequence', orfs: orfs };
    var html = '<p>' + int(orfs.length) + ' open reading frame' + (orfs.length === 1 ? '' : 's') + ' of at least ' + minAA +
      ' amino acids (ATG to stop codon, standard genetic code) in ' + int(s.length) + ' bp.</p>';
    if (orfs.length) {
      html += plotButtons('orf-map', 'orf_map') + '<div class="plot-box" id="orf-map"></div>';
      html += table(['#', 'Frame', 'Start', 'End', 'Length (bp)', 'Length (aa)', 'Protein'], orfs.slice(0, 100).map(function (o, i) {
        return [String(i + 1), (o.frame > 0 ? '+' : '') + o.frame, int(o.start), int(o.end), int(o.nt), int(o.aa),
          '<button type="button" class="btn btn-outline btn-small" data-orf="' + i + '">Show</button>'];
      }), [0, 2, 3, 4, 5]);
      html += '<p class="small muted">Positions refer to the entered strand; for minus-strand ORFs, Start is the lower coordinate. Lengths in bp include the stop codon.</p>';
      html += '<div id="orf-protein"></div>';
      html += '<p><button type="button" class="btn btn-outline btn-small" id="orf-dl">Download all ORF proteins (FASTA)</button></p>';
    }
    var fr = parseInt($('orf-frame').value, 10);
    var t = fr > 0 ? s.slice(fr - 1) : S.revcomp(s).slice(-fr - 1);
    var prot = S.translate(t, false);
    html += '<h4>Translation of frame ' + (fr > 0 ? '+' : '') + fr + '</h4><pre class="seq-out">' + esc(wrapSeq(prot, 60)) + '</pre>' +
      '<p class="small muted">* = stop codon; X = codon with N.</p>';
    show('orf-out', html);
    if (orfs.length) {
      orfMap($('orf-map'), orfs, s.length);
      Array.prototype.forEach.call(document.querySelectorAll('[data-orf]'), function (b) {
        b.addEventListener('click', function () {
          var o = orfs[+b.getAttribute('data-orf')];
          $('orf-protein').innerHTML = '<h4>ORF ' + (+b.getAttribute('data-orf') + 1) + ' protein (' + o.aa + ' aa)</h4><pre class="seq-out">' + esc(wrapSeq(o.protein, 60)) + '</pre>';
        });
      });
      $('orf-dl').addEventListener('click', function () {
        saveText(orfs.map(function (o, i) {
          return fasta(last.orfs.name + '_ORF' + (i + 1) + ' frame=' + o.frame + ' ' + o.start + '-' + o.end + ' ' + o.aa + 'aa', o.protein);
        }).join(''), 'orf_proteins.fasta');
      });
    }
  }

  function orfMap(box, orfs, L) {
    var W = Math.max(600, box.clientWidth || 800), H = 270;
    var svg = svgBox(box, W, H, 'ORF map');
    var left = 50, w = W - 70, top = 44, rowH = 30;
    txt(svg, left, 22, 'Open reading frames in six frames', { 'font-size': 14, 'font-weight': 700, fill: '#1f2328' });
    var Xs = function (v) { return left + (v - 1) / Math.max(1, L - 1) * w; };
    var frames = [1, 2, 3, -1, -2, -3];
    frames.forEach(function (f, i) {
      var y = top + i * rowH;
      el('line', { x1: left, x2: left + w, y1: y + 10, y2: y + 10, stroke: '#d0d7d3' }, svg);
      txt(svg, left - 8, y + 14, (f > 0 ? '+' : '') + f, { 'font-size': 11.5, fill: '#3d454d', 'text-anchor': 'end' });
    });
    orfs.forEach(function (o, i) {
      var y = top + frames.indexOf(o.frame) * rowH;
      var x1 = Xs(o.start), x2 = Xs(o.end);
      var rect = el('rect', { x: x1, y: y + 3, width: Math.max(2, x2 - x1), height: 14, fill: o.strand > 0 ? GREEN : BLUE, 'fill-opacity': i === 0 ? 0.95 : 0.6 }, svg);
      el('title', {}, rect).textContent = 'ORF ' + (i + 1) + ': ' + o.start + '–' + o.end + ', ' + o.aa + ' aa';
    });
    var axisY = top + 6 * rowH + 4;
    var step = niceStep(L, 8);
    for (var v = 0; v <= L; v += step) {
      var xx = Xs(Math.max(1, v));
      el('line', { x1: xx, x2: xx, y1: axisY, y2: axisY + 4, stroke: '#8c959f' }, svg);
      txt(svg, xx, axisY + 17, String(v), { 'font-size': 11, fill: '#57606a', 'text-anchor': 'middle' });
    }
    txt(svg, left + w / 2, axisY + 36, 'Position (bp)', { 'font-size': 12, 'font-weight': 600, fill: '#1f2328', 'text-anchor': 'middle' });
  }

  // ---------- 4. Codon usage ----------

  function runCodon() {
    var r = getSeq('codon-seq', 'dna', 'codon-out', 'a coding sequence');
    if (!r) return;
    var s = r.seq;
    var notes = [];
    if (s.length % 3) notes.push('The length (' + s.length + ' bp) is not a multiple of 3; the last ' + (s.length % 3) + ' base(s) are ignored.');
    if (s.slice(0, 3) !== 'ATG') notes.push('The sequence does not start with ATG.');
    var prot = S.translate(s, false);
    var inner = prot.slice(0, -1).indexOf('*');
    if (inner >= 0) notes.push('There is an internal stop codon at codon ' + (inner + 1) + '. Check the reading frame.');
    if (prot.slice(-1) !== '*') notes.push('The sequence does not end with a stop codon.');
    var u = S.codonUsage(s);
    last.codon = u;
    var gcPos = [0, 0, 0], nPos = [0, 0, 0];
    for (var i = 0; i + 3 <= s.length; i += 3) for (var k = 0; k < 3; k++) {
      var b = s.charAt(i + k);
      if (b === 'N') continue;
      nPos[k]++;
      if (b === 'G' || b === 'C') gcPos[k]++;
    }
    var pc = function (k) { return nPos[k] ? fmt(gcPos[k] / nPos[k] * 100, 1) + '%' : '–'; };
    var html = (notes.length ? '<div class="note small"><p>' + notes.map(esc).join('<br>') + '</p></div>' : '') + facts([
      ['Codons', int(u.codons) + (u.other ? ' (' + u.other + ' with N, not counted)' : '')],
      ['GC content', fmt(S.gcContent(s), 2) + '%'],
      ['GC at codon positions 1 / 2 / 3', pc(0) + ' / ' + pc(1) + ' / ' + pc(2)]
    ]);
    html += '<h4>Relative synonymous codon usage (RSCU)</h4><p class="small muted">RSCU = observed count ÷ count expected if all synonymous codons were used equally. 1 = no preference; above 1 = preferred.</p>';
    html += plotButtons('codon-plot', 'rscu') + '<div class="plot-box" id="codon-plot"></div>';
    html += '<p><button type="button" class="btn btn-small" id="codon-csv">Download table (CSV)</button></p>';
    html += table(['Amino acid', 'Codon', 'Count', 'Fraction', 'RSCU'], u.rows.map(function (x) {
      return [x.aa === '*' ? 'Stop' : x.aa, x.codon, int(x.count), x.fraction === null ? '–' : fmt(x.fraction, 2), x.rscu === null ? '–' : fmt(x.rscu, 2)];
    }), [2, 3, 4]);
    show('codon-out', html);
    var rows = u.rows.filter(function (x) { return x.aa !== 'M' && x.aa !== 'W'; });
    var palette = [GREEN, BLUE];
    var aaList = [];
    rows.forEach(function (x) { if (aaList.indexOf(x.aa) < 0) aaList.push(x.aa); });
    barPlot($('codon-plot'), { labels: rows.map(function (x) { return (x.aa === '*' ? 'Stop' : x.aa) + ' ' + x.codon; }),
      values: rows.map(function (x) { return x.rscu || 0; }), colors: rows.map(function (x) { return palette[aaList.indexOf(x.aa) % 2]; }),
      yLabel: 'RSCU', title: 'Relative synonymous codon usage', rotate: true, minBar: 13, height: 340, hline: 1 });
    $('codon-csv').addEventListener('click', function () {
      X.saveCsv([['amino_acid', 'codon', 'count', 'fraction', 'RSCU']].concat(last.codon.rows.map(function (x) {
        return [x.aa === '*' ? 'Stop' : x.aa, x.codon, x.count, x.fraction === null ? '' : x.fraction.toFixed(4), x.rscu === null ? '' : x.rscu.toFixed(4)];
      })), 'codon_usage.csv');
    });
  }

  // ---------- 5. Pairwise alignment ----------

  function runAlign() {
    var type = $('aln-type').value;
    var a = getSeq('aln-a', type, 'aln-out', 'sequence 1');
    if (!a) return;
    var b = getSeq('aln-b', type, 'aln-out', 'sequence 2');
    if (!b) return;
    if (a.seq.length * b.seq.length > MAX_ALIGN_CELLS) {
      error('aln-out', 'The sequences are too long for alignment in the browser (length 1 × length 2 must be below 4 million, for example 2,000 × 2,000).');
      return;
    }
    var opt = { type: type, mode: $('aln-mode').value, open: num('aln-open'), extend: num('aln-ext'), match: num('aln-match'), mismatch: num('aln-mis') };
    if (!(opt.open >= 0 && opt.extend >= 0)) { error('aln-out', 'Gap penalties must be 0 or positive numbers.'); return; }
    var r = S.align(a.seq, b.seq, opt);
    var name1 = (a.name || 'Seq1').split(/\s/)[0].slice(0, 12), name2 = (b.name || 'Seq2').split(/\s/)[0].slice(0, 12);
    var lines = [];
    var pa = r.startA, pb = r.startB, pad = Math.max(name1.length, name2.length, 4);
    for (var i = 0; i < r.length; i += 60) {
      var sa = r.a.slice(i, i + 60), sb = r.b.slice(i, i + 60), sm = r.mid.slice(i, i + 60);
      var ea = pa + sa.replace(/-/g, '').length - 1, eb = pb + sb.replace(/-/g, '').length - 1;
      lines.push(name1.padEnd(pad) + ' ' + String(pa).padStart(6) + ' ' + sa + ' ' + ea);
      lines.push(''.padEnd(pad) + '        ' + sm);
      lines.push(name2.padEnd(pad) + ' ' + String(pb).padStart(6) + ' ' + sb + ' ' + eb);
      lines.push('');
      pa = ea + 1; pb = eb + 1;
    }
    var pct = function (v) { return int(v) + '/' + int(r.length) + ' (' + fmt(r.length ? v / r.length * 100 : 0, 1) + '%)'; };
    var rows = [
      ['Mode', opt.mode === 'global' ? 'Global (Needleman–Wunsch), end gaps penalised' : 'Local (Smith–Waterman)'],
      ['Scoring', (type === 'protein' ? 'BLOSUM62' : 'match ' + opt.match + ', mismatch ' + opt.mismatch) + '; gap open ' + opt.open + ', gap extend ' + opt.extend],
      ['Score', fmt(r.score, 1)],
      ['Alignment length', int(r.length)],
      ['Identity', pct(r.identity)]
    ];
    if (type === 'protein') rows.push(['Similarity', pct(r.similarity)]);
    rows.push(['Gaps', pct(r.gaps)]);
    last.aln = lines.join('\n');
    show('aln-out', facts(rows) + '<pre class="seq-out aln">' + esc(last.aln) + '</pre>' +
      '<p class="small muted">| identical; ' + (type === 'protein' ? ': similar (positive BLOSUM62 score); ' : '') + '. different; space = gap.</p>' +
      '<p><button type="button" class="btn btn-outline btn-small" id="aln-dl">Download alignment (text)</button></p>');
    $('aln-dl').addEventListener('click', function () { saveText(last.aln, 'alignment.txt'); });
  }

  function alnDefaults() {
    var prot = $('aln-type').value === 'protein';
    $('aln-dna-opts').hidden = prot;
    $('aln-open').value = 10;
    $('aln-ext').value = 0.5;
  }

  // ---------- 6. Primer check and in-silico PCR ----------

  function primerSeq(id) {
    var s = $(id).value.replace(/\s/g, '').toUpperCase().replace(/U/g, 'T');
    return s;
  }

  function runPrimer() {
    var f = primerSeq('pr-f'), rv = primerSeq('pr-r');
    if (!f && !rv) { error('pr-out', 'Please enter at least one primer.'); return; }
    var bad = (f + rv).replace(/[ACGT]/g, '');
    if (bad) { error('pr-out', 'Primers may only contain A, C, G and T (found: ' + Array.from(new Set(bad.split(''))).join(', ') + ').'); return; }
    var opt = { na: num('pr-na'), mg: num('pr-mg'), dntps: num('pr-dntp'), primerNM: num('pr-conc') };
    if (!(opt.na > 0) || !(opt.primerNM > 0) || !(opt.mg >= 0) || !(opt.dntps >= 0)) { error('pr-out', 'Please enter valid concentrations (Na⁺ and primer above 0).'); return; }
    var list = [];
    if (f) list.push(['Forward', S.primerInfo(f, opt)]);
    if (rv) list.push(['Reverse', S.primerInfo(rv, opt)]);
    var flag = function (ok) { return ok ? '' : ' <span class="badge badge-medium">check</span>'; };
    var rows = list.map(function (x) {
      var p = x[1];
      return [x[0], '<span class="mono">' + esc(p.seq) + '</span>', String(p.length) + flag(p.length >= 18 && p.length <= 25),
        fmt(p.gc, 1) + flag(p.gc >= 40 && p.gc <= 60), fmt(p.tm, 1) + flag(p.tm >= 55 && p.tm <= 65), String(p.tmWallace),
        String(p.gcClamp) + flag(p.gcClamp >= 1 && p.gcClamp <= 3), String(p.run) + flag(p.run <= 4),
        p.self.any + ' / ' + p.self.threeEnd + flag(p.self.threeEnd < 4)];
    });
    var html = table(['Primer', 'Sequence (5′→3′)', 'Length', 'GC %', 'Tm (°C)', 'Tm Wallace', 'G/C in last 5', 'Longest run', 'Self-complementarity (any / 3′)'], rows, [2, 3, 4, 5, 6, 7, 8]);
    if (f && rv) {
      var pf = list[0][1], pr = list[1][1];
      var cross = S.complementarity(f, rv), cross2 = S.complementarity(rv, f);
      var dTm = Math.abs(pf.tm - pr.tm);
      html += facts([
        ['Tm difference', fmt(dTm, 1) + ' °C' + flag(dTm <= 3)],
        ['Primer–primer complementarity', 'longest stretch ' + cross.any + ' bp; at the 3′ end of forward ' + cross.threeEnd + ' bp, of reverse ' + cross2.threeEnd + ' bp' +
          flag(cross.threeEnd < 4 && cross2.threeEnd < 4)]
      ]);
    }
    var tmpl = S.readSeq($('pr-template').value, 'dna');
    if (tmpl.seq) {
      if (tmpl.invalid.length) html += '<p class="status error">The template contains characters that are not allowed: ' + esc(tmpl.invalid.join(', ')) + '.</p>';
      else if (!f || !rv) html += '<p class="muted">Enter both primers to find PCR products.</p>';
      else {
        var prods = S.pcrProducts(tmpl.seq, f, rv);
        html += '<h4>In-silico PCR on the template (' + int(tmpl.seq.length) + ' bp)</h4>';
        if (!prods.length) html += '<p>No product: the primers do not both match the template exactly in a facing orientation.</p>';
        else {
          html += table(['Product', 'Start', 'End', 'Size (bp)', 'GC %'], prods.map(function (p, i) {
            return [String(i + 1), int(p.start), int(p.end), int(p.length), fmt(S.gcContent(p.product), 1)];
          }), [1, 2, 3, 4]);
          if (prods.length > 1) html += '<p class="note small">More than one product is predicted.</p>';
          html += '<pre class="seq-out">' + esc(wrapSeq(prods[0].product, 60)) + '</pre>';
        }
        html += '<p class="small muted">Only exact matches of the full primer sequences are found; products with mismatched primers are not predicted.</p>';
      }
    }
    html += '<p class="small muted">Tm: nearest-neighbour method (Allawi &amp; SantaLucia 1997) with salt correction (SantaLucia 1998) and Mg²⁺ as Na⁺ equivalent (von Ahsen et al. 2001), calculated as Biopython Tm_NN. "check" marks values outside common guidelines: length 18–25, GC 40–60%, Tm 55–65 °C, 1–3 G/C in the last five bases, runs of at most 4 identical bases, 3′ complementarity below 4 bp, Tm difference at most 3 °C.</p>';
    show('pr-out', html);
  }

  // ---------- 7. Restriction sites ----------

  function runRestriction() {
    var r = getSeq('re-seq', 'dna', 're-out', 'a DNA sequence');
    if (!r) return;
    var s = r.seq;
    var chosen = Array.prototype.filter.call(document.querySelectorAll('#re-enzymes input'), function (c) { return c.checked; }).map(function (c) { return c.value; });
    if (!chosen.length) { error('re-out', 'Please select at least one enzyme.'); return; }
    var res = D.enzymes.filter(function (e) { return chosen.indexOf(e.name) >= 0; }).map(function (e) {
      var cuts = S.cutSites(s, e);
      return { e: e, cuts: cuts, frags: cuts.length ? S.fragments(s.length, cuts) : [] };
    });
    var cutters = res.filter(function (x) { return x.cuts.length; }), non = res.filter(function (x) { return !x.cuts.length; });
    last.re = { res: res, L: s.length };
    var html = '<p>' + int(s.length) + ' bp, linear. ' + cutters.length + ' of ' + res.length + ' selected enzymes cut this sequence.</p>';
    if (cutters.length) {
      html += plotButtons('re-map', 'restriction_map') + '<div class="plot-box" id="re-map"></div>';
      html += table(['Enzyme', 'Site', 'Cuts', 'Cut positions', 'Fragments (bp)'], cutters.map(function (x) {
        return [esc(x.e.name), '<span class="mono">' + esc(siteLabel(x.e)) + '</span>', String(x.cuts.length),
          x.cuts.map(int).join(', '), x.frags.map(int).join(', ')];
      }), [2]);
      html += '<p class="small muted">Cut position = first base after the cut on the entered (top) strand. Fragment sizes are for a single-enzyme digest of the linear sequence.</p>';
    }
    html += '<h4>Enzymes that do not cut</h4><p>' + (non.length ? non.map(function (x) { return esc(x.e.name); }).join(', ') : 'None (all selected enzymes cut).') + '</p>';
    html += '<p><button type="button" class="btn btn-small" id="re-csv">Download results (CSV)</button></p>';
    show('re-out', html);
    if (cutters.length) reMap($('re-map'), cutters, s.length);
    $('re-csv').addEventListener('click', function () {
      X.saveCsv([['enzyme', 'site', 'cuts', 'cut_positions', 'fragments_bp']].concat(last.re.res.map(function (x) {
        return [x.e.name, x.e.site, x.cuts.length, x.cuts.join(' '), x.frags.join(' ')];
      })), 'restriction_sites.csv');
    });
  }

  // Recognition site with the cut marked for enzymes that cut inside the site.
  function siteLabel(e) {
    if (e.fst5 >= 0 && e.fst5 <= e.site.length && e.pal) return e.site.slice(0, e.fst5) + '^' + e.site.slice(e.fst5);
    return e.site + ' (' + e.fst5 + '/' + (e.fst5 - e.ovhg) + ')';
  }

  function reMap(box, cutters, L) {
    var shown = cutters.filter(function (x) { return x.cuts.length <= 4; });
    var W = Math.max(600, box.clientWidth || 800), rowH = 18;
    var H = 120 + shown.length * rowH;
    var svg = svgBox(box, W, H, 'Restriction map');
    var left = 30, w = W - 60, y0 = 70;
    txt(svg, left, 22, 'Restriction map (enzymes with 1–4 cuts)', { 'font-size': 14, 'font-weight': 700, fill: '#1f2328' });
    var Xs = function (v) { return left + (v - 1) / Math.max(1, L - 1) * w; };
    el('rect', { x: left, y: y0 - 6, width: w, height: 12, fill: '#e6efe9', stroke: GREEN }, svg);
    txt(svg, left, y0 + 24, '1', { 'font-size': 11, fill: '#57606a' });
    txt(svg, left + w, y0 + 24, int(L), { 'font-size': 11, fill: '#57606a', 'text-anchor': 'end' });
    shown.forEach(function (x, i) {
      var ly = y0 + 44 + i * rowH;
      x.cuts.forEach(function (c) {
        el('line', { x1: Xs(c), x2: Xs(c), y1: y0 - 10, y2: ly - 4, stroke: '#8c959f', 'stroke-width': 0.8 }, svg);
      });
      var lx = Xs(x.cuts[0]);
      var anchor = lx > left + w - 160 ? 'end' : 'start';
      txt(svg, lx + (anchor === 'end' ? -3 : 3), ly + 4, x.e.name + ' (' + x.cuts.map(int).join(', ') + ')', { 'font-size': 11, fill: '#1f2328', 'text-anchor': anchor });
    });
    if (shown.length < cutters.length) txt(svg, left, H - 10, 'Enzymes with more than 4 cuts are listed in the table only.', { 'font-size': 11, fill: '#57606a' });
  }

  // ---------- 8. qPCR relative expression ----------

  function parseTable(text) {
    return text.split(/\r?\n/).map(function (l) { return l.trim(); }).filter(Boolean).map(function (l) { return l.split(/\t|,|;|\s+/).filter(function (x) { return x !== ''; }); });
  }

  function qpcrRows(outId) {
    var lines = parseTable($('qp-data').value);
    if (lines.length && isNaN(parseFloat(lines[0][2]))) lines.shift(); // header row
    if (!lines.length) { error(outId, 'Please enter your Ct values.'); return null; }
    var nRef = lines[0].length - 3;
    if (nRef < 1) { error(outId, 'Each row needs at least: sample, group, target Ct, reference Ct.'); return null; }
    var rows = [];
    for (var i = 0; i < lines.length; i++) {
      var l = lines[i];
      if (l.length !== nRef + 3) { error(outId, 'Row ' + (i + 1) + ' has ' + l.length + ' columns; expected ' + (nRef + 3) + ' like the first row.'); return null; }
      var nums = l.slice(2).map(Number);
      if (nums.some(function (v) { return !isFinite(v) || v <= 0 || v > 50; })) { error(outId, 'Row ' + (i + 1) + ' has a Ct value that is not a number between 0 and 50.'); return null; }
      rows.push({ sample: l[0], group: l[1], target: nums[0], refs: nums.slice(1) });
    }
    return rows;
  }

  function qpcrGroups() {
    var rows = qpcrRows('qp-out');
    var sel = $('qp-control');
    var prev = sel.value;
    sel.innerHTML = '';
    if (!rows) return;
    var groups = [];
    rows.forEach(function (r) { if (groups.indexOf(r.group) < 0) groups.push(r.group); });
    groups.forEach(function (g) { var o = document.createElement('option'); o.value = g; o.textContent = g; sel.appendChild(o); });
    if (groups.indexOf(prev) >= 0) sel.value = prev;
    $('qp-out').hidden = true;
  }

  function runQpcr() {
    var rows = qpcrRows('qp-out');
    if (!rows) return;
    var control = $('qp-control').value;
    var et = num('qp-et'), er = num('qp-er');
    if (!(et > 1 && et <= 2.2 && er > 1 && er <= 2.2)) { error('qp-out', 'Amplification factors must be between 1 and 2.2 (2 = 100% efficiency).'); return; }
    if (!rows.some(function (r) { return r.group === control; })) { error('qp-out', 'Please choose the control group.'); return; }
    var res = S.relativeExpression(rows, control, et, er);
    last.qpcr = res;
    var useLog = $('qp-log').checked;
    var method = et === 2 && er === 2 ? '2<sup>−ΔΔCt</sup> (Livak &amp; Schmittgen 2001)' : 'efficiency-corrected (Pfaffl 2001)';
    var html = '<p>Method: ' + method + '. Reference genes: ' + rows[0].refs.length + (rows[0].refs.length > 1 ? ' (combined by geometric mean)' : '') +
      '. Calibrator: mean of the "' + esc(control) + '" group.</p>';
    html += plotButtons('qp-plot', 'relative_expression') + '<div class="plot-box" id="qp-plot"></div>';
    html += table(['Group', 'n', 'Mean fold change', 'SD', 'Mean log2 fold change', 'p (Welch t-test vs control)'], res.summary.map(function (g) {
      return [esc(g.group), String(g.n), fmt(g.meanFold, 3), fmt(g.sdFold, 3), fmt(g.meanLog2, 3), g.p === null ? '–' : (g.p < 0.001 ? g.p.toExponential(1) : fmt(g.p, 4))];
    }), [1, 2, 3, 4, 5]);
    html += '<p class="small muted">The t-test compares log2 fold changes of the biological replicates with the control group; it needs at least 2 samples per group.</p>';
    html += table(['Sample', 'Group', 'ΔCt (target − mean reference)', 'Fold change', 'log2 fold change'], res.samples.map(function (s) {
      return [esc(s.sample), esc(s.group), fmt(s.dCt, 2), fmt(s.fold, 3), fmt(s.log2, 3)];
    }), [2, 3, 4]);
    html += '<p><button type="button" class="btn btn-small" id="qp-csv">Download results (CSV)</button></p>';
    show('qp-out', html);
    var groups = res.summary;
    var samples = res.samples;
    barPlot($('qp-plot'), {
      labels: groups.map(function (g) { return g.group; }),
      values: groups.map(function (g) { return useLog ? Math.max(0, g.meanLog2) : g.meanFold; }),
      errors: useLog ? null : groups.map(function (g) { return g.sdFold; }),
      points: useLog ? null : groups.map(function (g) { return samples.filter(function (s) { return s.group === g.group; }).map(function (s) { return s.fold; }); }),
      colors: groups.map(function (g) { return g.group === control ? GREY : GREEN; }),
      yLabel: useLog ? 'log2 fold change (positive part)' : 'Relative expression (fold change)', title: 'Relative expression', hline: useLog ? undefined : 1
    });
    if (useLog && groups.some(function (g) { return g.meanLog2 < 0; })) {
      $('qp-plot').insertAdjacentHTML('afterend', '<p class="small muted" id="qp-lognote">Negative log2 values (down-regulation) are shown as 0 in this simple plot; see the table for the values.</p>');
    }
    $('qp-csv').addEventListener('click', function () {
      X.saveCsv([['sample', 'group', 'dCt', 'fold_change', 'log2_fold_change']].concat(last.qpcr.samples.map(function (s) {
        return [s.sample, s.group, s.dCt.toFixed(4), s.fold.toFixed(6), s.log2.toFixed(6)];
      })), 'relative_expression.csv');
    });
  }

  // ---------- 9. qPCR standard curve ----------

  function runCurve() {
    var lines = parseTable($('sc-data').value);
    if (lines.length && isNaN(parseFloat(lines[0][0]))) lines.shift();
    var pts = [];
    for (var i = 0; i < lines.length; i++) {
      var q = Number(lines[i][0]), ct = Number(lines[i][1]);
      if (lines[i].length < 2 || !(q > 0) || !isFinite(ct)) { error('sc-out', 'Row ' + (i + 1) + ' needs a quantity above 0 and a Ct value.'); return; }
      pts.push({ qty: q, ct: ct });
    }
    var distinct = new Set(pts.map(function (p) { return p.qty; })).size;
    if (distinct < 3) { error('sc-out', 'Please enter at least 3 different dilutions.'); return; }
    var r = S.standardCurve(pts);
    var eff = r.efficiency * 100;
    var html = facts([
      ['Slope', fmt(r.slope, 4)],
      ['Intercept', fmt(r.intercept, 3)],
      ['R²', fmt(r.r2, 4) + (r.r2 < 0.98 ? ' <span class="badge badge-medium">below 0.98</span>' : '')],
      ['Efficiency', fmt(eff, 1) + '%' + (eff < 90 || eff > 110 ? ' <span class="badge badge-medium">outside 90–110%</span>' : '')],
      ['Amplification factor', fmt(1 + r.efficiency, 3) + ' (use this in the relative expression tool)']
    ]);
    html += plotButtons('sc-plot', 'standard_curve') + '<div class="plot-box" id="sc-plot"></div>';
    html += '<p class="small muted">Linear regression of Ct on log10(quantity). Efficiency = 10<sup>(−1/slope)</sup> − 1. A slope of −3.32 means 100% efficiency.</p>';
    show('sc-out', html);
    var box = $('sc-plot');
    var W = Math.max(600, box.clientWidth || 800), H = 320;
    var svg = svgBox(box, W, H, 'Standard curve');
    var f = { left: 56, top: 40, w: W - 80, h: H - 90, yLabel: 'Ct', title: 'Standard curve', ymin: Math.min.apply(null, r.y) - 0.5, ymax: Math.max.apply(null, r.y) + 0.5 };
    var Y = frame(svg, f);
    var xmin = Math.min.apply(null, r.x), xmax = Math.max.apply(null, r.x), pad = (xmax - xmin) * 0.05;
    var Xs = xAxis(svg, f, xmin - pad, xmax + pad, 'log10(quantity)');
    el('line', { x1: Xs(xmin - pad), x2: Xs(xmax + pad), y1: Y(r.intercept + r.slope * (xmin - pad)), y2: Y(r.intercept + r.slope * (xmax + pad)), stroke: RED, 'stroke-width': 1.5 }, svg);
    r.x.forEach(function (x, i) { el('circle', { cx: Xs(x), cy: Y(r.y[i]), r: 4, fill: GREEN, 'fill-opacity': 0.8 }, svg); });
    txt(svg, f.left + f.w - 4, f.top + 14, 'y = ' + fmt(r.slope, 3) + 'x + ' + fmt(r.intercept, 2) + ',  R² = ' + fmt(r.r2, 4) + ',  E = ' + fmt(eff, 1) + '%',
      { 'font-size': 12, fill: '#1f2328', 'text-anchor': 'end' });
  }

  // ---------- Examples ----------

  var EX = {
    protein: '>example artificial protein\nMASSEKKVLVLGAGGVGKSALTVQFVQGIFVEKYDPTIEDSYRKQVEVDCQQCMLEILDTAGTEQFTAMRDLYMKNGQGFALVYSITAQSTFNDLQDLREQILRVKDTEDVPMILVGNKCDLEDERVVGKEQGQNLARQWNNCAFLESSAKSKINVNEIFYDLVRQINRKTPVPGKARKKSSCQLL',
    dna: '>example artificial sequence\nGGCACGAGGTTCAAGCTTCCATGGCTTCTGCAAACGTTCTGGCAGCACTGGTGATCGTCGCCTTCATGTCAGTTTCGGCAACTGCAGAAGGCGAGTTCGACGGTCACCGCAATGGCGGCAAGCTCTACAACGCCATCTCGTGGGACACCGCCGAGGCAGCCATCGCAGCCGTCTCGCAGGGCTACGTTCCAGCCGGATCCTTCGAGACCGCAGAAGCCCTCGCCGGCTTCATCGCGTTCAAGCAGGCCAAGGGATCGTGACGTCATGGCGGCCGCTCGAGTTTAAACTATATAAGCCTGCAGGCTTGAGAATTCATAATGGCAATGG',
    primerF: 'CTTCTGCAAACGTTCTGGCA',
    primerR: 'CGATCCCTTGGCCTGCTTGA',
    qpcr: 'Sample\tGroup\tTarget\tRef1\tRef2\nC1\tControl\t26.10\t19.85\t21.02\nC2\tControl\t25.84\t19.70\t20.91\nC3\tControl\t26.31\t20.02\t21.20\nS1\tSalt\t23.42\t19.91\t21.05\nS2\tSalt\t23.10\t19.76\t20.88\nS3\tSalt\t23.65\t20.10\t21.31\nD1\tDrought\t24.90\t19.80\t21.10\nD2\tDrought\t25.21\t19.95\t21.00\nD3\tDrought\t24.75\t19.88\t21.15',
    curve: 'Quantity\tCt\n1\t17.12\n1\t17.20\n0.1\t20.45\n0.1\t20.51\n0.01\t23.80\n0.01\t23.73\n0.001\t27.14\n0.001\t27.08\n0.0001\t30.52\n0.0001\t30.44'
  };

  // ---------- Setup ----------

  var TOOLS = ['protein', 'dna', 'orf', 'codon', 'align', 'primer', 'restriction', 'qpcr', 'curve'];
  function select(name) {
    if (TOOLS.indexOf(name) < 0) name = 'protein';
    TOOLS.forEach(function (t) {
      $('tool-' + t).hidden = t !== name;
      var b = document.querySelector('[data-tool="' + t + '"]');
      b.setAttribute('aria-pressed', t === name ? 'true' : 'false');
    });
  }
  Array.prototype.forEach.call(document.querySelectorAll('[data-tool]'), function (b) {
    b.addEventListener('click', function () {
      var t = b.getAttribute('data-tool');
      select(t);
      if (history.replaceState) history.replaceState(null, '', '#' + t);
    });
  });
  select(location.hash.slice(1));
  window.addEventListener('hashchange', function () { select(location.hash.slice(1)); });

  var forms = { 'form-protein': [runProtein, 'prot-out'], 'form-dna': [runDna, 'dna-out'], 'form-orf': [runOrf, 'orf-out'],
    'form-codon': [runCodon, 'codon-out'], 'form-align': [runAlign, 'aln-out'], 'form-primer': [runPrimer, 'pr-out'],
    'form-restriction': [runRestriction, 're-out'], 'form-qpcr': [runQpcr, 'qp-out'], 'form-curve': [runCurve, 'sc-out'] };
  Object.keys(forms).forEach(function (id) {
    $(id).addEventListener('submit', function (e) {
      e.preventDefault();
      try { forms[id][0](); } catch (err) { error(forms[id][1], 'Something went wrong: ' + err.message); }
    });
  });

  // Restriction enzyme checkboxes
  var reBox = $('re-enzymes');
  D.enzymes.slice().sort(function (a, b) { return a.name.localeCompare(b.name); }).forEach(function (e) {
    var l = document.createElement('label');
    l.innerHTML = '<input type="checkbox" value="' + e.name + '" checked> ' + e.name;
    l.title = e.site;
    reBox.appendChild(l);
  });
  $('re-all').addEventListener('click', function () { Array.prototype.forEach.call(reBox.querySelectorAll('input'), function (c) { c.checked = true; }); });
  $('re-none').addEventListener('click', function () { Array.prototype.forEach.call(reBox.querySelectorAll('input'), function (c) { c.checked = false; }); });

  $('aln-type').addEventListener('change', alnDefaults);
  alnDefaults();
  $('qp-data').addEventListener('change', qpcrGroups);
  $('qp-data').addEventListener('blur', qpcrGroups);

  var fill = {
    'ex-protein': function () { $('prot-seq').value = EX.protein; runProtein(); },
    'ex-dna': function () { $('dna-seq').value = EX.dna; runDna(); },
    'ex-orf': function () { $('orf-seq').value = EX.dna; $('orf-min').value = 30; runOrf(); },
    'ex-codon': function () {
      // The longest ORF of the example sequence, as a complete coding sequence.
      var s = S.readSeq(EX.dna, 'dna').seq;
      var o = S.findOrfs(s, 30)[0];
      var cds = o.strand > 0 ? s.slice(o.start - 1, o.end) : S.revcomp(s.slice(o.start - 1, o.end));
      $('codon-seq').value = '>example coding sequence\n' + cds;
      runCodon();
    },
    'ex-align': function () {
      $('aln-type').value = 'protein'; alnDefaults();
      $('aln-a').value = '>protein_1\n' + EX.protein.split('\n')[1];
      $('aln-b').value = '>protein_2\n' + EX.protein.split('\n')[1].replace('KSALTVQ', 'KSAITVQ').replace('NDLQDLREQ', 'NDLHDLREQIL').slice(0, 150);
      runAlign();
    },
    'ex-primer': function () { $('pr-f').value = EX.primerF; $('pr-r').value = EX.primerR; $('pr-template').value = EX.dna; runPrimer(); },
    'ex-restriction': function () { $('re-seq').value = EX.dna; runRestriction(); },
    'ex-qpcr': function () { $('qp-data').value = EX.qpcr; qpcrGroups(); $('qp-control').value = 'Control'; runQpcr(); },
    'ex-curve': function () { $('sc-data').value = EX.curve; runCurve(); }
  };
  Object.keys(fill).forEach(function (id) { $(id).addEventListener('click', fill[id]); });

  // Plot downloads (buttons are created with the results)
  document.addEventListener('click', function (e) {
    var b = e.target.closest ? e.target.closest('[data-svg],[data-png]') : null;
    if (!b) return;
    var box = $(b.getAttribute('data-svg') || b.getAttribute('data-png'));
    var svg = box && box.querySelector('svg');
    if (!svg) return;
    if (b.hasAttribute('data-svg')) X.saveSvg(svg, b.getAttribute('data-name') + '.svg');
    else X.savePng(svg, b.getAttribute('data-name') + '.png');
  });
})();

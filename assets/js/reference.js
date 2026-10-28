/* qPCR reference gene finder.
   Ranked list: data/abiotic-rnaseq/reference_genes.json (built by tools/build_reference.py).
   Checking entered genes: TPM values loaded with genedata.js, scored in the browser. */
(function () {
  'use strict';

  var NS = 'http://www.w3.org/2000/svg';
  var FONT = 'Segoe UI, Roboto, Helvetica Neue, Arial, sans-serif';
  var STRESS_COLORS = { Salt: '#2f6f9f', Drought: '#b7791f', Heat: '#b83b3b' };
  var LINE_COLORS = ['#1b7837', '#2166ac', '#b2182b', '#7b3294', '#e08214', '#35978f', '#8c510a', '#de77ae', '#4d4d4d', '#5aae61'];
  var MAX_CHECK = 20;

  var $ = function (id) { return document.getElementById(id); };
  var meta = null;
  var ref = null;
  var state = { stress: 'All', tissue: 'All', minMean: 10, shown: 25, checked: null };

  function fmt(n) { return n.toLocaleString('en-US'); }
  function mean(a) { return a.reduce(function (s, v) { return s + v; }, 0) / a.length; }
  function sd(a) {
    if (a.length < 2) return 0;
    var m = mean(a);
    return Math.sqrt(a.reduce(function (s, v) { return s + (v - m) * (v - m); }, 0) / (a.length - 1));
  }
  function log2p(v) { return Math.log(v + 1) / Math.LN2; }
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
  function cell(tr, value, cls) {
    var td = document.createElement('td');
    if (value instanceof Node) td.appendChild(value); else td.textContent = value;
    if (cls) td.className = cls;
    tr.appendChild(td);
    return td;
  }
  // Main description only (text before the first ';'); the full text is shown on hover.
  function descCell(tr, desc) {
    var td = cell(tr, desc ? desc.split(';')[0].replace(' [Glycine max]', '') : '–', 'desc-cell');
    if (desc) td.title = desc;
    return td;
  }
  function geneLink(gene) {
    var a = document.createElement('a');
    a.href = '/expression.html?genes=' + encodeURIComponent(gene);
    a.textContent = gene;
    a.className = 'mono';
    a.title = 'Show the expression bar plot';
    return a;
  }
  function setStatus(id, msg, isError) {
    var s = $(id);
    s.textContent = msg;
    s.classList.toggle('error', !!isError);
  }
  function subsetKey() { return state.stress + '|' + state.tissue; }
  function subsetName() {
    var s = state.stress === 'All' ? 'all stresses' : state.stress.toLowerCase() + ' stress';
    var t = state.tissue === 'All' ? 'all tissues' : state.tissue.toLowerCase();
    return s + ', ' + t;
  }
  function inSubset(s) {
    return (state.stress === 'All' || s.stress === state.stress) && (state.tissue === 'All' || s.tissue === state.tissue);
  }

  // ---------- Ranked list ----------

  function renderTop() {
    var sub = ref.subsets[subsetKey()];
    var tbody = $('top-table').tBodies[0];
    tbody.innerHTML = '';
    if (!sub) {
      $('top-summary').textContent = 'There are no samples for this combination. Please choose another stress or tissue.';
      $('top-wrap').hidden = true;
      $('top-more').hidden = true;
      return;
    }
    $('top-wrap').hidden = false;
    var rows = sub.top[String(state.minMean)];
    $('top-summary').textContent = 'Most stable genes in ' + subsetName() + ': ' + fmt(sub.samples) + ' samples, ' +
      fmt(sub.expressed) + ' expressed genes ranked. Showing genes with mean TPM ≥ ' + state.minMean + '.';
    rows.slice(0, state.shown).forEach(function (r, i) {
      var info = ref.genes[r[0]] || ['', ''];
      var tr = document.createElement('tr');
      cell(tr, String(i + 1), 'num');
      cell(tr, geneLink(r[0]), 'nowrap');
      descCell(tr, info[1]);
      cell(tr, r[1].toFixed(1), 'num');
      cell(tr, r[2].toFixed(1), 'num');
      cell(tr, r[3].toFixed(3), 'num');
      cell(tr, r[4].toFixed(2), 'num');
      cell(tr, r[5] + ' / ' + sub.contrasts, 'num');
      tbody.appendChild(tr);
    });
    $('top-more').hidden = rows.length <= state.shown;
    $('top-count').textContent = rows.length ? '' : 'No genes match. Try a lower expression level.';
  }

  function topCsv() {
    var sub = ref.subsets[subsetKey()];
    if (!sub) return;
    var rows = [['rank', 'gene', 'ZH13_gene', 'description', 'mean_TPM', 'CV_%', 'SD_log2_TPM', 'max_fold_change', 'DE_comparisons', 'comparisons_in_subset']];
    sub.top[String(state.minMean)].forEach(function (r, i) {
      var info = ref.genes[r[0]] || ['', ''];
      rows.push([i + 1, r[0], info[0], info[1], r[1], r[2], r[3], r[4], r[5], sub.contrasts]);
    });
    saveCsv(rows, 'reference_genes_' + subsetKey().replace('|', '_') + '.csv');
  }

  // ---------- Checking entered genes ----------

  function score(g, idx, cons) {
    var tpm = idx.map(function (i) { return g.t[i]; });
    var lg = tpm.map(log2p);
    var m = mean(tpm);
    var s = sd(lg);
    var lo = Math.min.apply(null, tpm), hi = Math.max.apply(null, tpm);
    var de = 0;
    cons.forEach(function (ci) {
      var f = g.f[ci], p = g.p[ci];
      if (f !== null && p !== null && p < 0.05 && Math.abs(f) >= 1) de++;
    });
    var q = ref.subsets[subsetKey()].q;
    var pct = 0;
    while (pct < 100 && q[pct] < s) pct++;
    return {
      gene: g.gene, desc: g.desc, zh: g.zh, lg: lg, mean: m, cv: m > 0 ? sd(tpm) / m * 100 : NaN, sd: s,
      fold: (hi + 1) / (lo + 1), de: de, pct: pct, lowMin: lo < ref.min_tpm,
      expressed: m >= ref.min_mean
    };
  }

  // geNorm M: mean standard deviation of log2 ratios with every other gene (lower = more stable).
  function genormM(list) {
    return list.map(function (a, i) {
      var sds = [];
      list.forEach(function (b, j) {
        if (i === j) return;
        sds.push(sd(a.lg.map(function (v, k) { return v - b.lg[k]; })));
      });
      return mean(sds);
    });
  }

  function verdict(r) {
    if (!r.expressed) return { text: 'Low expression', cls: 'badge-none' };
    if (r.lowMin) return { text: 'Variable (near zero in some samples)', cls: 'badge-up' };
    if (r.pct <= 10 && r.de === 0) return { text: 'Very stable', cls: 'badge-high' };
    if (r.pct <= 25) return { text: 'Stable', cls: 'badge-high' };
    if (r.pct <= 50) return { text: 'Moderate', cls: 'badge-medium' };
    return { text: 'Variable', cls: 'badge-up' };
  }

  function runCheck() {
    var tokens = SoybioData.parseTokens($('check-genes').value);
    if (!tokens.length) { setStatus('check-status', 'Please enter at least one gene ID.', true); return; }
    if (tokens.length > MAX_CHECK) { setStatus('check-status', 'Please enter up to ' + MAX_CHECK + ' genes.', true); return; }
    setStatus('check-status', 'Loading expression data…');
    SoybioData.resolve(tokens).then(function (items) {
      var seen = {};
      state.checked = items.filter(function (it) {
        if (!it.data || seen[it.data.gene]) return false;
        seen[it.data.gene] = true;
        return true;
      });
      var missing = items.filter(function (it) { return !it.data; });
      $('check-missing').hidden = !missing.length;
      $('check-missing').textContent = missing.length ? 'Not found: ' + missing.map(function (it) { return it.query + ' (' + it.note + ')'; }).join('; ') + '.' : '';
      if (!state.checked.length) {
        $('check-results').hidden = true;
        setStatus('check-status', 'None of the entered IDs were found in the expression data.', true);
        return;
      }
      renderCheck();
      setStatus('check-status', 'Checked ' + state.checked.length + (state.checked.length === 1 ? ' gene.' : ' genes.'));
    }).catch(function () {
      setStatus('check-status', 'The expression data could not be loaded. Please reload the page.', true);
    });
  }

  function renderCheck() {
    if (!state.checked || !state.checked.length) return;
    var sub = ref.subsets[subsetKey()];
    $('check-results').hidden = false;
    var tbody = $('check-table').tBodies[0];
    tbody.innerHTML = '';
    if (!sub) {
      $('check-summary').textContent = 'There are no samples for this combination. Please choose another stress or tissue.';
      $('check-plot').innerHTML = '';
      return;
    }
    var idx = [];
    meta.samples.forEach(function (s, i) { if (inSubset(s)) idx.push(i); });
    var cons = [];
    meta.contrasts.forEach(function (c, i) {
      if ((state.stress === 'All' || c.stress === state.stress) && (state.tissue === 'All' || c.tissue === state.tissue)) cons.push(i);
    });
    var res = state.checked.map(function (it) { var r = score(it.data, idx, cons); r.query = it.query; r.from = it.from; return r; });
    var M = res.length >= 3 ? genormM(res) : null;
    res.forEach(function (r, i) { r.m = M ? M[i] : null; });
    res.sort(function (a, b) { return a.sd - b.sd; });
    state.checkResults = res;
    $('check-summary').textContent = 'Stability in ' + subsetName() + ' (' + fmt(idx.length) + ' samples). Genes are sorted from most to least stable.';
    res.forEach(function (r) {
      var tr = document.createElement('tr');
      var g = cell(tr, geneLink(r.gene), 'nowrap');
      if (r.from && r.from !== 'Wm82.a4') {
        var sm = document.createElement('div');
        sm.className = 'small muted';
        sm.textContent = 'from ' + r.query;
        g.appendChild(sm);
      }
      descCell(tr, r.desc);
      cell(tr, r.mean.toFixed(1), 'num');
      cell(tr, isNaN(r.cv) ? '–' : r.cv.toFixed(1), 'num');
      cell(tr, r.sd.toFixed(3), 'num');
      cell(tr, r.m === null ? '–' : r.m.toFixed(3), 'num');
      cell(tr, r.de + ' / ' + cons.length, 'num');
      cell(tr, r.expressed ? (100 - r.pct) + '%' : '–', 'num');
      var v = verdict(r);
      var b = document.createElement('span');
      b.className = 'badge ' + v.cls;
      b.textContent = v.text;
      cell(tr, b);
      tbody.appendChild(tr);
    });
    $('genorm-note').hidden = !!M;
    renderLines(res, idx);
  }

  // Line plot: mean log2(TPM + 1) per condition group, one line per gene.
  function renderLines(res, idx) {
    var box = $('check-plot');
    box.innerHTML = '';
    var groups = [];
    meta.groups.forEach(function (g) {
      if ((state.stress === 'All' || g.stress === state.stress) && (state.tissue === 'All' || g.tissue === state.tissue)) groups.push(g);
    });
    var pos = {};
    idx.forEach(function (si, k) { pos[si] = k; });
    var shown = res.slice(0, LINE_COLORS.length);
    var series = shown.map(function (r) {
      return groups.map(function (g) { return mean(g.samples.map(function (si) { return r.lg[pos[si]]; })); });
    });
    var W = Math.max(640, box.clientWidth), H = 330;
    var left = 54, right = 16, top = 56, bottom = 46;
    var pw = W - left - right, ph = H - top - bottom;
    var vmax = 1;
    series.forEach(function (s) { s.forEach(function (v) { vmax = Math.max(vmax, v); }); });
    vmax = Math.ceil(vmax);
    var X = function (i) { return left + (groups.length === 1 ? pw / 2 : i / (groups.length - 1) * pw); };
    var Y = function (v) { return top + ph - v / vmax * ph; };
    var svg = el('svg', { xmlns: NS, width: W, height: H, viewBox: '0 0 ' + W + ' ' + H, 'font-family': FONT, role: 'img',
      'aria-label': 'Expression of the checked genes across conditions' });
    el('rect', { x: 0, y: 0, width: W, height: H, fill: '#ffffff' }, svg);
    var lx = left;
    shown.forEach(function (r, i) {
      el('line', { x1: lx, x2: lx + 16, y1: 16, y2: 16, stroke: LINE_COLORS[i], 'stroke-width': 2.5 }, svg);
      txt(svg, lx + 21, 20, r.gene, { 'font-size': 11.5, fill: '#1f2328' });
      lx += 138;
      if (lx > W - 140 && i < shown.length - 1) lx = W;
    });
    var step = vmax <= 6 ? 1 : 2;
    for (var v = 0; v <= vmax; v += step) {
      el('line', { x1: left, x2: left + pw, y1: Y(v), y2: Y(v), stroke: v === 0 ? '#8c959f' : '#e3e8e5' }, svg);
      txt(svg, left - 6, Y(v) + 4, String(v), { 'font-size': 11, fill: '#57606a', 'text-anchor': 'end' });
    }
    txt(svg, 14, top + ph / 2, 'log2(TPM + 1)', { 'font-size': 11.5, fill: '#1f2328', 'text-anchor': 'middle',
      transform: 'rotate(-90 14 ' + (top + ph / 2) + ')' });
    // Stress band under the axis
    groups.forEach(function (g, i) {
      var x0 = i === 0 ? left : (X(i - 1) + X(i)) / 2;
      var x1 = i === groups.length - 1 ? left + pw : (X(i) + X(i + 1)) / 2;
      el('rect', { x: x0, y: top + ph + 8, width: Math.max(1, x1 - x0), height: 7,
        fill: g.treatment === 'Control' ? '#9aa3a0' : STRESS_COLORS[g.stress] }, svg);
    });
    txt(svg, left, H - 10, 'Conditions (' + groups.length + '), grouped by stress and study. Band: grey = control, blue = salt, orange = drought, red = heat.',
      { 'font-size': 11, fill: '#57606a' });
    series.forEach(function (s, i) {
      el('polyline', { points: s.map(function (v, k) { return X(k) + ',' + Y(v); }).join(' '), fill: 'none',
        stroke: LINE_COLORS[i], 'stroke-width': 1.8, 'stroke-linejoin': 'round' }, svg);
    });
    box.appendChild(svg);
  }

  function checkCsv() {
    var res = state.checkResults;
    if (!res) return;
    var rows = [['gene', 'entered_ID', 'description', 'mean_TPM', 'CV_%', 'SD_log2_TPM', 'geNorm_M', 'max_fold_change', 'DE_comparisons', 'more_stable_than_%_of_genes', 'subset']];
    res.forEach(function (r) {
      rows.push([r.gene, r.query, r.desc, r.mean.toFixed(2), isNaN(r.cv) ? '' : r.cv.toFixed(1), r.sd.toFixed(4),
        r.m === null ? '' : r.m.toFixed(4), r.fold.toFixed(2), r.de, r.expressed ? 100 - r.pct : '', subsetName()]);
    });
    saveCsv(rows, 'reference_gene_check.csv');
  }

  function saveCsv(rows, name) {
    var text = '﻿' + rows.map(function (r) {
      return r.map(function (v) {
        var s = v === null || v === undefined ? '' : String(v);
        return /[",\r\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
      }).join(',');
    }).join('\r\n');
    var url = URL.createObjectURL(new Blob([text], { type: 'text/csv;charset=utf-8' }));
    var a = document.createElement('a');
    a.href = url;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 2000);
  }

  // ---------- Setup ----------

  function refresh() {
    renderTop();
    renderCheck();
  }

  Promise.all([SoybioData.loadMeta(), SoybioData.getJSON(SoybioData.BASE + 'reference_genes.json')]).then(function (r) {
    meta = r[0];
    ref = r[1];
    $('check-run').disabled = false;
    setStatus('check-status', 'Ready.');
    renderTop();
  }).catch(function () {
    $('top-summary').textContent = 'The data could not be loaded. Please reload the page.';
  });

  Array.prototype.forEach.call(document.querySelectorAll('#ref-stress button'), function (b) {
    b.addEventListener('click', function () {
      state.stress = b.getAttribute('data-stress');
      Array.prototype.forEach.call(document.querySelectorAll('#ref-stress button'), function (o) {
        o.setAttribute('aria-pressed', o === b ? 'true' : 'false');
      });
      state.shown = 25;
      refresh();
    });
  });
  $('ref-tissue').addEventListener('change', function () { state.tissue = this.value; state.shown = 25; refresh(); });
  $('ref-min').addEventListener('change', function () { state.minMean = +this.value; state.shown = 25; renderTop(); });
  $('top-more').addEventListener('click', function () { state.shown += 25; renderTop(); });
  $('top-csv').addEventListener('click', topCsv);
  $('check-form').addEventListener('submit', function (e) { e.preventDefault(); runCheck(); });
  $('check-csv').addEventListener('click', checkCsv);
  $('check-example').addEventListener('click', function () {
    var top = ref && ref.subsets['All|All'].top['10'];
    if (!top) return;
    // Three top-ranked genes plus two stress-responsive genes for comparison.
    $('check-genes').value = [top[0][0], top[2][0], top[5][0], 'Glyma.10G145300', 'Glyma.14G100000'].join('\n');
    runCheck();
  });
})();

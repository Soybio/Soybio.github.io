/* Expression browser for the abiotic stress RNA-seq data.
   One gene: bar plot of TPM per condition and a table of DE results.
   Several genes: heatmap of genes x conditions (or samples).
   Data files are built by tools/build_expression.py and loaded with genedata.js. Plots are drawn as SVG. */
(function () {
  'use strict';

  var MAX_GENES = 200;
  var NS = 'http://www.w3.org/2000/svg';
  var FONT = 'Segoe UI, Roboto, Helvetica Neue, Arial, sans-serif';
  var COLORS = { Control: '#9aa3a0', Salt: '#2f6f9f', Drought: '#b7791f', Heat: '#b83b3b' };
  var TISSUE_COLORS = { 'Leaf': '#5f9a4e', 'Root': '#8a6a43', 'Pods': '#d1a92a', 'Leaf, root and stem': '#7d8a8c' };
  // ColorBrewer "Reds" for expression and "RdBu" (reversed) for z-scores.
  var SEQ = ['#fff5f0', '#fee0d2', '#fcbba1', '#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d', '#99000d'];
  var DIV = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#f7f7f7', '#fddbc7', '#f4a582', '#d6604d', '#b2182b'];

  var $ = function (id) { return document.getElementById(id); };
  var meta = null;
  var state = { stress: 'All', tissue: 'All', genes: [], single: null };

  // ---------- Small helpers ----------

  function el(tag, attrs, parent) {
    var e = document.createElementNS(NS, tag);
    for (var k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }

  function txt(parent, x, y, str, attrs) {
    var a = { x: x, y: y };
    for (var k in attrs) a[k] = attrs[k];
    var t = el('text', a, parent);
    t.textContent = str;
    return t;
  }

  function fmt(n) { return n.toLocaleString('en-US'); }
  function mean(a) { return a.length ? a.reduce(function (s, v) { return s + v; }, 0) / a.length : 0; }
  function sd(a) {
    if (a.length < 2) return 0;
    var m = mean(a);
    return Math.sqrt(a.reduce(function (s, v) { return s + (v - m) * (v - m); }, 0) / (a.length - 1));
  }
  function log2p(v) { return Math.log(v + 1) / Math.LN2; }
  function radio(name) { var r = document.querySelector('input[name="' + name + '"]:checked'); return r ? r.value : null; }

  var measureCtx = document.createElement('canvas').getContext('2d');
  function fit(str, px, font) {
    measureCtx.font = font;
    if (measureCtx.measureText(str).width <= px) return str;
    while (str.length > 1 && measureCtx.measureText(str + '…').width > px) str = str.slice(0, -1);
    return str + '…';
  }

  function niceTicks(max, count) {
    var raw = max / count;
    var mag = Math.pow(10, Math.floor(Math.log10(raw)));
    var step = [1, 2, 2.5, 5, 10].map(function (m) { return m * mag; }).find(function (s) { return s >= raw; });
    var ticks = [];
    for (var v = 0; v <= max + step * 0.999; v += step) ticks.push(+v.toPrecision(12));
    return ticks;
  }

  function lerpColor(stops, t) {
    t = Math.max(0, Math.min(1, t));
    var p = t * (stops.length - 1);
    var i = Math.min(Math.floor(p), stops.length - 2);
    var f = p - i;
    var a = hex(stops[i]), b = hex(stops[i + 1]);
    return 'rgb(' + [0, 1, 2].map(function (k) { return Math.round(a[k] + (b[k] - a[k]) * f); }).join(',') + ')';
  }
  function hex(h) { return [1, 3, 5].map(function (i) { return parseInt(h.slice(i, i + 2), 16); }); }

  function setStatus(msg, isError) {
    var s = $('status');
    s.textContent = msg;
    s.classList.toggle('error', !!isError);
  }

  // ---------- Tooltip ----------

  var tip = $('plot-tooltip');
  function showTip(evt, html) {
    tip.innerHTML = html;
    tip.hidden = false;
    var x = evt.clientX + 14, y = evt.clientY + 14;
    var w = tip.offsetWidth, h = tip.offsetHeight;
    if (x + w > window.innerWidth - 8) x = evt.clientX - w - 14;
    if (y + h > window.innerHeight - 8) y = evt.clientY - h - 14;
    tip.style.left = x + 'px';
    tip.style.top = y + 'px';
  }
  function hideTip() { tip.hidden = true; }
  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]; }); }

  // Gene data loading and ID lookup are in genedata.js (window.SoybioData).
  var resolve = SoybioData.resolve;

  // ---------- Filters ----------

  function groupVisible(g) {
    return (state.stress === 'All' || g.stress === state.stress) &&
      (state.tissue === 'All' || g.tissue === state.tissue);
  }
  function activeGroups() {
    var out = [];
    meta.groups.forEach(function (g, i) { if (groupVisible(g)) out.push(i); });
    return out;
  }
  function activeSamples() {
    var out = [];
    activeGroups().forEach(function (gi) { out.push.apply(out, meta.groups[gi].samples); });
    return out;
  }
  function studyKey(g) { return g.project + ' · ' + g.tissue; }

  // ---------- Bar plot (one gene) ----------

  function renderBar(g) {
    var box = $('bar-plot');
    box.innerHTML = '';
    var useLog = radio('bar-scale') === 'log';
    var f = useLog ? log2p : function (v) { return v; };
    var gis = activeGroups();
    if (!gis.length) { box.textContent = 'No samples match the selected filters.'; return; }

    var rows = [];
    var lastStress = null, lastStudy = null;
    gis.forEach(function (gi) {
      var grp = meta.groups[gi];
      if (grp.stress !== lastStress) { rows.push({ type: 'stress', grp: grp }); lastStress = grp.stress; lastStudy = null; }
      if (studyKey(grp) !== lastStudy) { rows.push({ type: 'study', grp: grp }); lastStudy = studyKey(grp); }
      var vals = grp.samples.map(function (i) { return f(g.t[i]); });
      rows.push({ type: 'bar', grp: grp, vals: vals, tpm: grp.samples.map(function (i) { return g.t[i]; }), mean: mean(vals), sd: sd(vals) });
    });

    var W = Math.max(660, box.clientWidth);
    var labelW = Math.round(Math.min(330, Math.max(230, W * 0.3)));
    var plotX = labelW + 12, plotW = W - plotX - 28;
    var top = 78, rowH = 18, stressH = 26, studyH = 21;
    var bodyH = rows.reduce(function (s, r) { return s + (r.type === 'bar' ? rowH : r.type === 'stress' ? stressH : studyH); }, 0);
    var H = top + bodyH + 54;

    var maxV = 0;
    rows.forEach(function (r) { if (r.type === 'bar') maxV = Math.max(maxV, r.mean + r.sd, Math.max.apply(null, r.vals)); });
    if (!(maxV > 0)) maxV = 1;
    var ticks = niceTicks(maxV, 5);
    var xmax = ticks[ticks.length - 1];
    var X = function (v) { return plotX + Math.max(0, v) / xmax * plotW; };

    var svg = el('svg', { xmlns: NS, width: W, height: H, viewBox: '0 0 ' + W + ' ' + H, 'font-family': FONT, role: 'img',
      'aria-label': 'Bar plot of ' + g.gene + ' expression' });
    el('rect', { x: 0, y: 0, width: W, height: H, fill: '#ffffff' }, svg);

    // Title and legend
    txt(svg, 8, 22, g.gene, { 'font-size': 16, 'font-weight': 700, fill: '#1f2328' });
    var tw = (function () { measureCtx.font = '700 16px ' + FONT; return measureCtx.measureText(g.gene).width; })();
    if (g.desc) txt(svg, 8 + tw + 10, 22, fit(g.desc, W - tw - 30, '12px ' + FONT), { 'font-size': 12, fill: '#57606a' });
    var lx = 8;
    ['Control', 'Salt', 'Drought', 'Heat'].forEach(function (k) {
      el('rect', { x: lx, y: 36, width: 12, height: 12, rx: 2, fill: COLORS[k] }, svg);
      txt(svg, lx + 17, 46, k === 'Control' ? 'Control' : k + ' stress', { 'font-size': 12, fill: '#1f2328' });
      lx += k === 'Control' ? 80 : 100;
    });
    txt(svg, lx + 8, 46, 'Bars: mean ± SD; dots: samples', { 'font-size': 12, fill: '#57606a' });

    // Axis and grid
    var axisTitle = useLog ? 'Expression, log2(TPM + 1)' : 'Expression (TPM)';
    var yTop = top - 8, yBot = top + bodyH + 4;
    ticks.forEach(function (t) {
      var x = X(t);
      el('line', { x1: x, x2: x, y1: yTop, y2: yBot, stroke: t === 0 ? '#8c959f' : '#e3e8e5', 'stroke-width': 1 }, svg);
      txt(svg, x, yTop - 5, fmtTick(t), { 'font-size': 11, fill: '#57606a', 'text-anchor': 'middle' });
      txt(svg, x, yBot + 15, fmtTick(t), { 'font-size': 11, fill: '#57606a', 'text-anchor': 'middle' });
    });
    txt(svg, plotX + plotW / 2, yBot + 36, axisTitle, { 'font-size': 12, fill: '#1f2328', 'text-anchor': 'middle', 'font-weight': 600 });

    // Rows
    var y = top;
    var jit = [0, -3, 3, -1.5, 1.5, -4, 4, -2.2, 2.2, -0.7, 0.7, -3.5, 3.5, -1, 1];
    rows.forEach(function (r) {
      if (r.type === 'stress') {
        txt(svg, 8, y + 18, r.grp.stress + ' stress', { 'font-size': 13.5, 'font-weight': 700, fill: COLORS[r.grp.stress] });
        y += stressH;
        return;
      }
      if (r.type === 'study') {
        txt(svg, 14, y + 15, fit(studyKey(r.grp) + ' · ' + r.grp.stage, labelW + plotW, '600 11.5px ' + FONT),
          { 'font-size': 11.5, 'font-weight': 600, fill: '#3d454d' });
        y += studyH;
        return;
      }
      var cy = y + rowH / 2;
      var color = r.grp.treatment === 'Control' ? COLORS.Control : COLORS[r.grp.stress];
      txt(svg, labelW, cy + 4, fit(r.grp.label, labelW - 24, '12px ' + FONT), { 'font-size': 12, fill: '#1f2328', 'text-anchor': 'end' });
      el('rect', { x: plotX, y: cy - 5.5, width: Math.max(0, X(r.mean) - plotX), height: 11, fill: color, 'fill-opacity': 0.85 }, svg);
      if (r.sd > 0) {
        var x1 = X(r.mean - r.sd), x2 = X(r.mean + r.sd);
        el('line', { x1: x1, x2: x2, y1: cy, y2: cy, stroke: '#2b2f33', 'stroke-width': 1 }, svg);
        el('line', { x1: x2, x2: x2, y1: cy - 3.5, y2: cy + 3.5, stroke: '#2b2f33', 'stroke-width': 1 }, svg);
        if (r.mean - r.sd > 0) el('line', { x1: x1, x2: x1, y1: cy - 3.5, y2: cy + 3.5, stroke: '#2b2f33', 'stroke-width': 1 }, svg);
      }
      r.vals.forEach(function (v, j) {
        el('circle', { cx: X(v), cy: cy + jit[j % jit.length], r: 2.3, fill: '#1f2328', 'fill-opacity': 0.6 }, svg);
      });
      var hit = el('rect', { x: 0, y: y, width: W, height: rowH, fill: 'transparent' }, svg);
      hit.addEventListener('mousemove', function (e) {
        showTip(e, '<strong>' + esc(studyKey(r.grp)) + '</strong><br>' + esc(r.grp.label) +
          '<br>Mean TPM: ' + mean(r.tpm).toFixed(2) + ' (n = ' + r.tpm.length + ')' +
          '<br>Samples: ' + r.tpm.map(function (v) { return v.toFixed(2); }).join(', '));
      });
      hit.addEventListener('mouseleave', hideTip);
      y += rowH;
    });

    box.appendChild(svg);
  }

  function fmtTick(v) {
    if (v >= 10000) return (v / 1000) + 'k';
    return String(+v.toPrecision(6));
  }

  // ---------- DE table (one gene) ----------

  function fmtP(p) {
    if (p === null) return '–';
    if (p === 0) return '< 1e-300';
    if (p < 0.001) return p.toExponential(1);
    return p.toFixed(3);
  }

  function renderDE(g) {
    var tbody = $('de-table').tBodies[0];
    tbody.innerHTML = '';
    var up = 0, down = 0, tested = 0;
    meta.contrasts.forEach(function (c, i) {
      var lfc = g.f[i], p = g.p[i];
      var call = 'Not tested';
      if (lfc !== null && p !== null) {
        tested++;
        call = (p < 0.05 && Math.abs(lfc) >= 1) ? (lfc > 0 ? 'Up' : 'Down') : 'Not significant';
        if (call === 'Up') up++;
        if (call === 'Down') down++;
      }
      if (!((state.stress === 'All' || c.stress === state.stress) && (state.tissue === 'All' || c.tissue === state.tissue))) return;
      var tr = document.createElement('tr');
      [c.stress, c.project + ' · ' + c.tissue, c.condition, lfc === null ? '–' : lfc.toFixed(2), fmtP(p)].forEach(function (v, k) {
        var td = document.createElement('td');
        td.textContent = v;
        if (k >= 3) td.className = 'num';
        if (k === 1) td.className = 'nowrap';
        tr.appendChild(td);
      });
      var td = document.createElement('td');
      var b = document.createElement('span');
      b.className = 'badge ' + ({ Up: 'badge-up', Down: 'badge-down', 'Not significant': 'badge-none', 'Not tested': 'badge-none' })[call];
      b.textContent = call;
      td.appendChild(b);
      tr.appendChild(td);
      tbody.appendChild(tr);
    });
    $('de-summary').textContent = 'Differentially expressed in ' + (up + down) + ' of ' + meta.contrasts.length +
      ' stress-versus-control comparisons (' + up + ' up, ' + down + ' down; tested in ' + tested + ').';
  }

  function renderGeneInfo(g, item) {
    var box = $('gene-info');
    box.innerHTML = '';
    var h = document.createElement('h3');
    h.textContent = g.gene;
    box.appendChild(h);
    var dl = document.createElement('dl');
    dl.className = 'facts';
    var add = function (k, v) {
      var dt = document.createElement('dt'); dt.textContent = k;
      var dd = document.createElement('dd'); dd.textContent = v;
      dl.appendChild(dt); dl.appendChild(dd);
    };
    if (item && item.from && item.from !== 'Wm82.a4') add('Entered ID', item.query + ' (' + item.from + ')');
    add('ZH13 gene', g.zh || 'No matching ZH13 gene');
    add('Description', g.desc || 'No description');
    var all = g.t;
    add('TPM range', Math.min.apply(null, all).toFixed(2) + ' to ' + Math.max.apply(null, all).toFixed(2) +
      ' across ' + all.length + ' samples');
    box.appendChild(dl);
  }

  function showSingle(g, item) {
    state.single = { g: g, item: item };
    $('single').hidden = false;
    renderGeneInfo(g, item);
    renderBar(g);
    renderDE(g);
  }

  // ---------- Heatmap (several genes) ----------

  function clusterRows(M) {
    var n = M.length;
    var z = M.map(function (row) {
      var m = mean(row);
      var s = Math.sqrt(row.reduce(function (a, v) { return a + (v - m) * (v - m); }, 0));
      return row.map(function (v) { return s > 0 ? (v - m) / s : 0; });
    });
    var D = [];
    for (var i = 0; i < n; i++) {
      D.push(new Float64Array(n));
      for (var j = 0; j < i; j++) {
        var r = 0;
        for (var k = 0; k < z[i].length; k++) r += z[i][k] * z[j][k];
        D[i][j] = D[j][i] = 1 - r;
      }
    }
    var nodes = [];
    for (i = 0; i < n; i++) nodes.push({ leaf: i, size: 1, h: 0 });
    var active = nodes.map(function (_, idx) { return idx; });
    var dist = D.map(function (row) { return Array.prototype.slice.call(row); });
    while (active.length > 1) {
      var best = Infinity, bi = -1, bj = -1;
      for (var a = 0; a < active.length; a++) {
        for (var b = a + 1; b < active.length; b++) {
          var d = dist[active[a]][active[b]];
          if (d < best) { best = d; bi = a; bj = b; }
        }
      }
      var A = active[bi], B = active[bj];
      var na = nodes[A].size, nb = nodes[B].size;
      var id = nodes.length;
      nodes.push({ left: A, right: B, size: na + nb, h: Math.max(best, nodes[A].h, nodes[B].h) });
      dist[id] = [];
      active.forEach(function (c) {
        if (c === A || c === B) return;
        var v = (na * dist[A][c] + nb * dist[B][c]) / (na + nb);
        dist[id][c] = v;
        dist[c][id] = v;
      });
      active.splice(bj, 1);
      active.splice(bi, 1);
      active.push(id);
    }
    var order = [];
    (function walk(k) { var nd = nodes[k]; if (nd.leaf !== undefined) order.push(nd.leaf); else { walk(nd.left); walk(nd.right); } })(nodes.length - 1);
    return { order: order, nodes: nodes, root: nodes.length - 1 };
  }

  function heatmapMatrix(list) {
    var mode = radio('hm-values');
    var bySample = radio('hm-cols') === 'samples';
    var cols = bySample
      ? activeSamples().map(function (i) { var s = meta.samples[i]; return { label: s.id, sub: meta.groups[s.group].label, grp: meta.groups[s.group], idx: [i] }; })
      : activeGroups().map(function (gi) { var g = meta.groups[gi]; return { label: g.project + ' ' + g.tissue + ': ' + g.label, sub: '', grp: g, idx: g.samples }; });
    var raw = list.map(function (g) { return cols.map(function (c) { return mean(c.idx.map(function (i) { return g.t[i]; })); }); });
    var disp = raw.map(function (row) {
      if (mode === 'tpm') return row.slice();
      var lg = row.map(log2p);
      if (mode === 'log') return lg;
      var m = mean(lg), s = sd(lg);
      return lg.map(function (v) { return s > 0 ? (v - m) / s : 0; });
    });
    return { cols: cols, raw: raw, disp: disp, mode: mode };
  }

  function renderHeatmap() {
    var box = $('heatmap');
    box.innerHTML = '';
    var list = state.genes.map(function (it) { return it.data; });
    var hm = heatmapMatrix(list);
    var cols = hm.cols;
    if (!cols.length) { box.textContent = 'No samples match the selected filters.'; return; }

    var n = list.length;
    var clustered = radio('hm-order') === 'cluster' && n >= 3;
    var cl = clustered ? clusterRows(hm.disp) : null;
    var order = cl ? cl.order : list.map(function (_, i) { return i; });

    var avail = Math.max(640, box.clientWidth);
    var dendroW = clustered ? 64 : 0;
    var labelW = 260;
    var gap = 3;
    var studies = 0;
    cols.forEach(function (c, i) { if (i > 0 && studyKey(c.grp) !== studyKey(cols[i - 1].grp)) studies++; });
    var cw = Math.floor((avail - dendroW - labelW - 24 - studies * gap) / cols.length);
    cw = Math.max(4, Math.min(22, cw));
    var ch = n <= 40 ? 16 : n <= 100 ? 12 : 9;
    var rowFont = ch >= 16 ? 12 : ch >= 12 ? 10.5 : 8.5;
    var colX = [];
    var x0 = 10 + dendroW;
    var cx = x0;
    cols.forEach(function (c, i) {
      if (i > 0 && studyKey(c.grp) !== studyKey(cols[i - 1].grp)) cx += gap;
      colX.push(cx);
      cx += cw;
    });
    var gridW = cx - x0;
    var showColLabels = cw >= 9;
    var colFont = 10;
    var colLabelH = 0;
    if (showColLabels) {
      measureCtx.font = colFont + 'px ' + FONT;
      colLabelH = Math.min(260, Math.max.apply(null, cols.map(function (c) { return measureCtx.measureText(c.label).width; }))) + 10;
    }
    var top = 92;
    var stripH = 12;
    var yCells = top + stripH * 2 + 8;
    var H = yCells + n * ch + colLabelH + 26;
    var W = Math.max(avail, x0 + gridW + labelW + 10);

    var svg = el('svg', { xmlns: NS, width: W, height: H, viewBox: '0 0 ' + W + ' ' + H, 'font-family': FONT, role: 'img',
      'aria-label': 'Heatmap of ' + n + ' genes' });
    el('rect', { x: 0, y: 0, width: W, height: H, fill: '#ffffff' }, svg);

    // Color scale
    var vmax = 0;
    hm.disp.forEach(function (row) { row.forEach(function (v) { vmax = Math.max(vmax, Math.abs(v)); }); });
    var lo, hi, stops;
    if (hm.mode === 'z') { hi = Math.min(3, Math.max(1, Math.ceil(vmax * 2) / 2)); lo = -hi; stops = DIV; }
    else { lo = 0; hi = vmax > 0 ? vmax : 1; stops = SEQ; }
    var color = function (v) { return lerpColor(stops, (v - lo) / (hi - lo)); };
    var title = { log: 'log2(TPM + 1)', z: 'Z-score of log2(TPM + 1) per gene', tpm: 'TPM' }[hm.mode];

    txt(svg, 10, 22, n + ' genes × ' + cols.length + (radio('hm-cols') === 'samples' ? ' samples' : ' conditions (mean of replicates)'),
      { 'font-size': 15, 'font-weight': 700, fill: '#1f2328' });
    txt(svg, 10, 40, 'Values: ' + title + (clustered ? '. Genes clustered by expression pattern (average linkage, correlation distance).' : '.'),
      { 'font-size': 12, fill: '#57606a' });

    // Legends
    var gx = 10, gy = 56;
    var defs = el('defs', {}, svg);
    var lg = el('linearGradient', { id: 'hm-grad', x1: 0, x2: 1, y1: 0, y2: 0 }, defs);
    stops.forEach(function (s, i) { el('stop', { offset: (i / (stops.length - 1) * 100) + '%', 'stop-color': s }, lg); });
    el('rect', { x: gx, y: gy, width: 150, height: 10, fill: 'url(#hm-grad)', stroke: '#c4cbc7' }, svg);
    txt(svg, gx, gy + 23, fmtTick(+lo.toFixed(2)), { 'font-size': 10.5, fill: '#57606a' });
    txt(svg, gx + 150, gy + 23, fmtTick(+hi.toFixed(2)), { 'font-size': 10.5, fill: '#57606a', 'text-anchor': 'end' });
    var lx = gx + 180;
    ['Control', 'Salt', 'Drought', 'Heat'].forEach(function (k) {
      el('rect', { x: lx, y: gy, width: 10, height: 10, fill: COLORS[k] }, svg);
      txt(svg, lx + 14, gy + 9, k, { 'font-size': 11, fill: '#1f2328' });
      lx += k === 'Control' ? 66 : k === 'Drought' ? 66 : 50;
    });
    lx += 14;
    Object.keys(TISSUE_COLORS).forEach(function (k) {
      if (!cols.some(function (c) { return c.grp.tissue === k; })) return;
      el('rect', { x: lx, y: gy, width: 10, height: 10, fill: TISSUE_COLORS[k] }, svg);
      txt(svg, lx + 14, gy + 9, k, { 'font-size': 11, fill: '#1f2328' });
      measureCtx.font = '11px ' + FONT;
      lx += 26 + measureCtx.measureText(k).width;
    });

    // Column strips: condition (control or stress) and tissue
    txt(svg, x0 + gridW + 8, top + stripH - 3, 'Condition', { 'font-size': 10, fill: '#57606a' });
    txt(svg, x0 + gridW + 8, top + stripH * 2 - 3, 'Tissue', { 'font-size': 10, fill: '#57606a' });
    cols.forEach(function (c, i) {
      el('rect', { x: colX[i], y: top, width: cw, height: stripH - 1, fill: c.grp.treatment === 'Control' ? COLORS.Control : COLORS[c.grp.stress] }, svg);
      el('rect', { x: colX[i], y: top + stripH, width: cw, height: stripH - 1, fill: TISSUE_COLORS[c.grp.tissue] || '#999' }, svg);
    });

    // Cells
    var cells = el('g', { 'shape-rendering': 'crispEdges' }, svg);
    order.forEach(function (r, ri) {
      var y = yCells + ri * ch;
      hm.disp[r].forEach(function (v, c) {
        el('rect', { x: colX[c], y: y, width: cw, height: ch, fill: color(v) }, cells);
      });
    });
    el('rect', { x: x0 - 0.5, y: yCells - 0.5, width: gridW + 1, height: n * ch + 1, fill: 'none', stroke: '#c4cbc7' }, svg);

    // Row labels (click for the bar plot)
    var lxr = x0 + gridW + 8;
    order.forEach(function (r, ri) {
      var g = list[r];
      var y = yCells + ri * ch + ch / 2 + rowFont * 0.36;
      var label = g.gene + (g.desc && ch >= 12 ? '  ' + g.desc : '');
      var t = txt(svg, lxr, y, fit(label, W - lxr - 8, rowFont + 'px ' + FONT), { 'font-size': rowFont, fill: '#1a5a8a', cursor: 'pointer' });
      t.addEventListener('click', function () { showSingle(g, state.genes[r]); $('single').scrollIntoView({ behavior: 'smooth' }); });
    });

    // Column labels
    if (showColLabels) {
      var yl = yCells + n * ch + 6;
      cols.forEach(function (c, i) {
        var xx = colX[i] + cw / 2 + colFont * 0.35;
        txt(svg, xx, yl, fit(c.label, colLabelH - 10, colFont + 'px ' + FONT),
          { 'font-size': colFont, fill: '#3d454d', 'text-anchor': 'end', transform: 'rotate(-90 ' + xx + ' ' + yl + ')' });
      });
    } else {
      txt(svg, x0, yCells + n * ch + 18, 'Hover over the heatmap to see sample names.', { 'font-size': 11, fill: '#57606a' });
    }

    // Dendrogram
    if (cl) {
      var maxH = cl.nodes[cl.root].h || 1;
      var leafY = {};
      order.forEach(function (r, ri) { leafY[r] = yCells + ri * ch + ch / 2; });
      var dx = function (h) { return 8 + (dendroW - 6) * (1 - h / maxH); };
      var dg = el('g', { stroke: '#57606a', 'stroke-width': 1, fill: 'none' }, svg);
      (function draw(k) {
        var nd = cl.nodes[k];
        if (nd.leaf !== undefined) return { y: leafY[nd.leaf], x: dx(0) };
        var L = draw(nd.left), R = draw(nd.right);
        var x = dx(nd.h);
        el('path', { d: 'M' + L.x + ' ' + L.y + 'H' + x + 'V' + R.y + 'H' + R.x }, dg);
        return { y: (L.y + R.y) / 2, x: x };
      })(cl.root);
    }

    // Hover over cells
    var hit = el('rect', { x: x0, y: yCells, width: gridW, height: n * ch, fill: 'transparent' }, svg);
    hit.addEventListener('mousemove', function (e) {
      var pt = svg.getBoundingClientRect();
      var scale = pt.width / W;
      var mx = (e.clientX - pt.left) / scale, my = (e.clientY - pt.top) / scale;
      var ri = Math.floor((my - yCells) / ch);
      var c = -1;
      for (var i = 0; i < colX.length; i++) if (mx >= colX[i] && mx < colX[i] + cw) { c = i; break; }
      if (ri < 0 || ri >= n || c < 0) { hideTip(); return; }
      var r = order[ri];
      var col = cols[c];
      showTip(e, '<strong>' + esc(list[r].gene) + '</strong><br>' + esc(col.label) + (col.sub ? ' (' + esc(col.sub) + ')' : '') +
        '<br>' + esc(title) + ': ' + hm.disp[r][c].toFixed(2) + '<br>Mean TPM: ' + hm.raw[r][c].toFixed(2));
    });
    hit.addEventListener('mouseleave', hideTip);

    box.appendChild(svg);
    state.heatmap = { hm: hm, order: order, list: list };
  }

  // ---------- Downloads ----------

  function saveBlob(blob, name) {
    var url = URL.createObjectURL(blob);
    var a = document.createElement('a');
    a.href = url;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 2000);
  }

  function svgText(svg) {
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(svg);
  }

  function saveSvg(boxId, name) {
    var svg = $(boxId).querySelector('svg');
    if (svg) saveBlob(new Blob([svgText(svg)], { type: 'image/svg+xml;charset=utf-8' }), name);
  }

  function savePng(boxId, name) {
    var svg = $(boxId).querySelector('svg');
    if (!svg) return;
    var w = +svg.getAttribute('width'), h = +svg.getAttribute('height');
    var img = new Image();
    var url = URL.createObjectURL(new Blob([svgText(svg)], { type: 'image/svg+xml;charset=utf-8' }));
    img.onload = function () {
      var c = document.createElement('canvas');
      c.width = w * 2;
      c.height = h * 2;
      var ctx = c.getContext('2d');
      ctx.scale(2, 2);
      ctx.fillStyle = '#ffffff';
      ctx.fillRect(0, 0, w, h);
      ctx.drawImage(img, 0, 0, w, h);
      URL.revokeObjectURL(url);
      c.toBlob(function (b) { saveBlob(b, name); }, 'image/png');
    };
    img.src = url;
  }

  function csv(rows) {
    return '﻿' + rows.map(function (r) {
      return r.map(function (v) {
        var s = v === null || v === undefined ? '' : String(v);
        return /[",\r\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
      }).join(',');
    }).join('\r\n');
  }

  function saveSingleCsv() {
    var g = state.single && state.single.g;
    if (!g) return;
    var rows = [['sample', 'SRA_run', 'BioProject', 'stress', 'treatment', 'tissue', 'stage', 'duration', 'intensity', 'condition_group', 'TPM']];
    meta.samples.forEach(function (s, i) {
      rows.push([s.id, s.srr, s.project, s.stress, s.treatment, s.tissue, s.stage, s.duration, s.intensity, meta.groups[s.group].label, g.t[i]]);
    });
    saveBlob(new Blob([csv(rows)], { type: 'text/csv;charset=utf-8' }), g.gene + '_TPM.csv');
  }

  function saveHeatmapCsv() {
    var h = state.heatmap;
    if (!h) return;
    var title = { log: 'log2(TPM+1)', z: 'z-score', tpm: 'TPM' }[h.hm.mode];
    var rows = [['gene (' + title + ')', 'description'].concat(h.hm.cols.map(function (c) { return c.label; }))];
    h.order.forEach(function (r) {
      rows.push([h.list[r].gene, h.list[r].desc].concat(h.hm.disp[r].map(function (v) { return +v.toFixed(4); })));
    });
    saveBlob(new Blob([csv(rows)], { type: 'text/csv;charset=utf-8' }), 'heatmap_values.csv');
  }

  // ---------- Main actions ----------

  function plot() {
    var tokens = SoybioData.parseTokens($('genes').value);
    if (!tokens.length) { setStatus('Please enter at least one gene ID.', true); return; }
    if (tokens.length > MAX_GENES) { setStatus('Please enter up to ' + MAX_GENES + ' genes at a time (you entered ' + tokens.length + ').', true); return; }
    setStatus('Loading expression data…');
    $('plot').disabled = true;
    resolve(tokens).then(function (items) {
      $('plot').disabled = false;
      var found = [];
      var seenGene = {};
      items.forEach(function (it) { if (it.data && !seenGene[it.data.gene]) { seenGene[it.data.gene] = true; found.push(it); } });
      var missing = items.filter(function (it) { return !it.data; });
      var note = $('missing');
      if (missing.length) {
        note.hidden = false;
        note.textContent = 'Not shown: ' + missing.map(function (it) { return it.query + ' (' + it.note + ')'; }).join('; ') + '.';
      } else {
        note.hidden = true;
      }
      if (!found.length) {
        setStatus('None of the entered IDs were found in the expression data.', true);
        $('multi').hidden = true;
        $('single').hidden = true;
        return;
      }
      state.genes = found;
      $('results').hidden = false;
      // Put the gene list in the page address so the result can be shared as a link.
      if (window.history && history.replaceState) {
        var q = tokens.length <= 60 ? '?genes=' + encodeURIComponent(tokens.join(',')) : location.pathname;
        history.replaceState(null, '', q);
      }
      if (found.length === 1) {
        $('multi').hidden = true;
        showSingle(found[0].data, found[0]);
        setStatus('Showing 1 gene.');
      } else {
        $('multi').hidden = false;
        $('single').hidden = true;
        state.single = null;
        renderHeatmap();
        setStatus('Showing ' + found.length + ' genes. Click a gene name in the heatmap to see its bar plot.');
      }
      $('results').scrollIntoView({ behavior: 'smooth', block: 'start' });
    }).catch(function () {
      $('plot').disabled = false;
      setStatus('The expression data could not be loaded. Please reload the page and try again.', true);
    });
  }

  function rerender() {
    if (!meta) return;
    if (!$('multi').hidden && state.genes.length > 1) renderHeatmap();
    if (state.single) { renderBar(state.single.g); renderDE(state.single.g); }
  }

  function setupFilters() {
    var tissues = [];
    meta.groups.forEach(function (g) { if (tissues.indexOf(g.tissue) < 0) tissues.push(g.tissue); });
    var sel = $('tissue-filter');
    tissues.forEach(function (t) {
      var o = document.createElement('option');
      o.value = t;
      o.textContent = t;
      sel.appendChild(o);
    });
    sel.addEventListener('change', function () { state.tissue = sel.value; rerender(); });
    Array.prototype.forEach.call(document.querySelectorAll('#stress-filter button'), function (b) {
      b.addEventListener('click', function () {
        state.stress = b.getAttribute('data-stress');
        Array.prototype.forEach.call(document.querySelectorAll('#stress-filter button'), function (o) {
          o.setAttribute('aria-pressed', o === b ? 'true' : 'false');
        });
        rerender();
      });
    });
  }

  SoybioData.loadMeta().then(function (m) {
    meta = m;
    setupFilters();
    $('plot').disabled = false;
    setStatus('Ready. ' + fmt(meta.genes) + ' genes and ' + fmt(meta.samples.length) + ' samples available.');
    var m2 = /[?&]genes=([^&]*)/.exec(location.search);
    if (m2) {
      $('genes').value = decodeURIComponent(m2[1].replace(/\+/g, ' ')).split(/[\s,;]+/).join('\n');
      plot();
    }
  }).catch(function () {
    setStatus('The expression data could not be loaded. Please reload the page.', true);
  });

  $('expr-form').addEventListener('submit', function (e) { e.preventDefault(); plot(); });
  $('genes').addEventListener('keydown', function (e) { if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) plot(); });
  $('ex-single').addEventListener('click', function () { if (!meta) return; $('genes').value = meta.examples.single; plot(); });
  $('ex-multi').addEventListener('click', function () { if (!meta) return; $('genes').value = meta.examples.multiple.join('\n'); plot(); });
  $('clear').addEventListener('click', function () {
    $('genes').value = '';
    $('results').hidden = true;
    state.genes = [];
    state.single = null;
    $('genes').focus();
  });
  ['bar-scale', 'hm-values', 'hm-cols', 'hm-order'].forEach(function (name) {
    Array.prototype.forEach.call(document.querySelectorAll('input[name="' + name + '"]'), function (r) {
      r.addEventListener('change', rerender);
    });
  });
  $('bar-svg').addEventListener('click', function () { saveSvg('bar-plot', (state.single ? state.single.g.gene : 'gene') + '_barplot.svg'); });
  $('bar-png').addEventListener('click', function () { savePng('bar-plot', (state.single ? state.single.g.gene : 'gene') + '_barplot.png'); });
  $('bar-csv').addEventListener('click', saveSingleCsv);
  $('hm-svg').addEventListener('click', function () { saveSvg('heatmap', 'heatmap.svg'); });
  $('hm-png').addEventListener('click', function () { savePng('heatmap', 'heatmap.png'); });
  $('hm-csv').addEventListener('click', saveHeatmapCsv);
  var resizeTimer = null;
  var lastWidth = window.innerWidth;
  window.addEventListener('resize', function () {
    if (window.innerWidth === lastWidth) return;
    lastWidth = window.innerWidth;
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(rerender, 250);
  });
})();

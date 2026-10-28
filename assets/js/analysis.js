/* Gene Analysis page.
   Gene list: enrichment of GO, KEGG and Pfam terms and of stress-responsive genes.
   One gene: annotation, stress response, co-expressed genes and a function prediction
   from the annotation of the co-expressed genes.
   Data: data/annotation/ (tools/build_annotation.py) and data/abiotic-rnaseq/coexp/
   (tools/build_coexpression.py). */
(function () {
  'use strict';

  var ANN = '/data/annotation/';
  var NS = 'http://www.w3.org/2000/svg';
  var FONT = 'Segoe UI, Roboto, Helvetica Neue, Arial, sans-serif';
  var MAX_LIST = 3000;
  var MIN_TERM = 5;      // terms with fewer background genes are not tested
  var DOTS = 15;         // terms shown in the dot plot
  var SOURCES = {
    go_bp: { file: 'go_bp.json', name: 'GO Biological Process' },
    go_mf: { file: 'go_mf.json', name: 'GO Molecular Function' },
    go_cc: { file: 'go_cc.json', name: 'GO Cellular Component' },
    kegg: { file: 'kegg.json', name: 'KEGG pathways' },
    pfam: { file: 'pfam.json', name: 'Pfam domains' }
  };
  var STRESSES = ['Salt', 'Drought', 'Heat'];
  var FDR_COLORS = ['#fde0c5', '#f9a874', '#ec6b4f', '#c2364d', '#7a1d4e'];

  var $ = function (id) { return document.getElementById(id); };
  var cache = {};
  var meta = null;
  var state = { mode: null, list: null, source: 'go_bp', results: {}, single: null };

  // ---------- Helpers ----------

  function load(name) {
    if (!cache[name]) cache[name] = SoybioData.getJSON(ANN + name);
    return cache[name];
  }
  function fmt(n) { return n.toLocaleString('en-US'); }
  function fmtP(p) {
    if (p === 0) return '< 1e-300';
    if (p < 0.001) return p.toExponential(1);
    return p.toFixed(3);
  }
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
  var measureCtx = document.createElement('canvas').getContext('2d');
  function fit(s, px, font) {
    measureCtx.font = font;
    if (measureCtx.measureText(s).width <= px) return s;
    while (s.length > 1 && measureCtx.measureText(s + '…').width > px) s = s.slice(0, -1);
    return s + '…';
  }
  function cell(tr, value, cls) {
    var td = document.createElement('td');
    if (value instanceof Node) td.appendChild(value); else td.textContent = value;
    if (cls) td.className = cls;
    tr.appendChild(td);
    return td;
  }
  function link(href, text, cls) {
    var a = document.createElement('a');
    a.href = href;
    a.textContent = text;
    if (cls) a.className = cls;
    return a;
  }
  function setStatus(msg, isError) {
    var s = $('status');
    s.textContent = msg;
    s.classList.toggle('error', !!isError);
  }
  function shortDesc(d) { return d ? d.split(';')[0].replace(' [Glycine max]', '') : ''; }

  // ---------- Statistics ----------

  var logFact = [0];
  function lf(n) {
    for (var i = logFact.length; i <= n; i++) logFact.push(logFact[i - 1] + Math.log(i));
    return logFact[n];
  }
  function lchoose(n, k) { return lf(n) - lf(k) - lf(n - k); }

  // P(X >= k) for X ~ Hypergeometric(N population, K successes, n draws).
  function hyperUpper(k, K, n, N) {
    var hi = Math.min(K, n);
    if (k > hi) return 0;
    var base = lchoose(N, n);
    var terms = [];
    for (var i = k; i <= hi; i++) terms.push(lchoose(K, i) + lchoose(N - K, n - i) - base);
    var m = Math.max.apply(null, terms);
    var s = terms.reduce(function (a, t) { return a + Math.exp(t - m); }, 0);
    return Math.min(1, Math.exp(m) * s);
  }

  // Benjamini-Hochberg adjusted p-values, in the same order as ps.
  function bh(ps) {
    var order = ps.map(function (p, i) { return i; }).sort(function (a, b) { return ps[a] - ps[b]; });
    var q = new Array(ps.length);
    var prev = 1;
    for (var r = order.length - 1; r >= 0; r--) {
      var i = order[r];
      prev = Math.min(prev, ps[i] * ps.length / (r + 1));
      q[i] = prev;
    }
    return q;
  }

  // Over-representation test of a gene list against an annotation source.
  // db: {terms: [[id, name, global?]], genes: {gene: [term index]}}; bg: Set of background genes or null.
  function enrich(db, list, bg) {
    var universe = Object.keys(db.genes).filter(function (g) { return !bg || bg.has(g); });
    var inUniverse = {};
    universe.forEach(function (g) { inUniverse[g] = true; });
    var hits = list.filter(function (g) { return inUniverse[g]; });
    var N = universe.length, n = hits.length;
    var K = new Int32Array(db.terms.length);
    universe.forEach(function (g) { db.genes[g].forEach(function (t) { K[t]++; }); });
    var k = {}, members = {};
    hits.forEach(function (g) {
      db.genes[g].forEach(function (t) {
        k[t] = (k[t] || 0) + 1;
        (members[t] = members[t] || []).push(g);
      });
    });
    var rows = [];
    Object.keys(k).forEach(function (ts) {
      var t = +ts;
      var term = db.terms[t];
      if (K[t] < MIN_TERM || term[2] === 1 || k[t] < 2) return;
      rows.push({ id: term[0], name: term[1], k: k[t], K: K[t], n: n, N: N,
        fold: (k[t] / n) / (K[t] / N), p: hyperUpper(k[t], K[t], n, N), genes: members[t] });
    });
    var q = bh(rows.map(function (r) { return r.p; }));
    rows.forEach(function (r, i) { r.q = q[i]; });
    rows.sort(function (a, b) { return a.p - b.p || b.fold - a.fold; });
    return { rows: rows, n: n, N: N, listSize: list.length };
  }

  // ---------- ID mapping (no expression data needed) ----------

  function mapIds(tokens) {
    return Promise.all([load('genes.json'), SoybioData.loadIdmap()]).then(function (r) {
      var all = {};
      r[0].forEach(function (g) { all[g] = true; });
      var map = r[1];
      var genes = [], seen = {}, missing = [], converted = 0;
      tokens.forEach(function (q) {
        var c = SoybioData.canonical(q);
        var g = null;
        if (c.type === 'wm') g = all[c.id] ? c.id : (map.a2[c.id] || null);
        else if (c.type === 'zh') g = map.zh[c.id] || null;
        if (g && all[g]) {
          if (g !== c.id) converted++;
          if (!seen[g]) { seen[g] = true; genes.push(g); }
        } else {
          missing.push(q);
        }
      });
      return { genes: genes, missing: missing, converted: converted };
    });
  }

  // ---------- Gene list analysis ----------

  function runList(tokens) {
    var bgTokens = SoybioData.parseTokens($('background').value);
    Promise.all([mapIds(tokens), bgTokens.length ? mapIds(bgTokens) : null]).then(function (r) {
      var m = r[0];
      if (m.genes.length < 2) {
        setStatus(m.genes.length ? 'Only one gene was recognised. For a single gene, enter just that gene.' :
          'None of the entered IDs were recognised.', true);
        return;
      }
      state.mode = 'list';
      state.list = m.genes;
      state.bg = r[1] ? new Set(r[1].genes.concat(m.genes)) : null;
      state.results = {};
      $('single').hidden = true;
      $('list').hidden = false;
      $('results').hidden = false;
      var parts = [fmt(tokens.length) + ' IDs entered', fmt(m.genes.length) + ' genes recognised'];
      if (m.converted) parts.push(fmt(m.converted) + ' converted from ZH13 or Wm82.a2.v1 IDs');
      if (m.missing.length) parts.push(fmt(m.missing.length) + ' not found');
      $('list-summary').textContent = parts.join(', ') + '.' + (state.bg ? ' Background: ' + fmt(state.bg.size) + ' genes you entered.' : '');
      $('list-missing').hidden = !m.missing.length;
      $('list-missing').textContent = m.missing.length ? 'Not found: ' + m.missing.slice(0, 50).join(', ') + (m.missing.length > 50 ? ' …' : '') : '';
      var heat = $('list-heatmap');
      heat.hidden = m.genes.length > 200;
      heat.href = '/expression.html?genes=' + encodeURIComponent(m.genes.slice(0, 60).join(','));
      heat.textContent = m.genes.length > 60 ? 'Expression heatmap of the first 60 genes' : 'Expression heatmap of these genes';
      setStatus('Analysing ' + fmt(m.genes.length) + ' genes…');
      return Promise.all([renderStress(), renderSource()]).then(function () {
        setStatus('Done. ' + fmt(m.genes.length) + ' genes analysed.');
        $('results').scrollIntoView({ behavior: 'smooth', block: 'start' });
      });
    }).catch(function () { setStatus('The annotation data could not be loaded. Please reload the page.', true); });
  }

  function renderStress() {
    return Promise.all([load('stress.json'), load('genes.json')]).then(function (r) {
      var st = r[0];
      var universe = state.bg ? r[1].filter(function (g) { return state.bg.has(g); }) : r[1];
      var inU = {};
      universe.forEach(function (g) { inU[g] = true; });
      var list = state.list.filter(function (g) { return inU[g]; });
      var N = universe.length, n = list.length;
      var tests = [];
      STRESSES.forEach(function (s, si) {
        [['responsive (up or down)', null], ['up-regulated', 0], ['down-regulated', 1]].forEach(function (c) {
          var has = function (g) {
            var v = st.genes[g];
            if (!v) return false;
            return c[1] === null ? v[si * 2] + v[si * 2 + 1] > 0 : v[si * 2 + c[1]] > 0;
          };
          var K = universe.filter(has).length;
          var k = list.filter(has);
          tests.push({ stress: s, cat: c[0], k: k.length, K: K, n: n, N: N, genes: k,
            fold: n && K ? (k.length / n) / (K / N) : 0, p: hyperUpper(k.length, K, n, N) });
        });
      });
      var q = bh(tests.map(function (t) { return t.p; }));
      tests.forEach(function (t, i) { t.q = q[i]; });
      state.stressTests = tests;
      var tbody = $('stress-table').tBodies[0];
      tbody.innerHTML = '';
      tests.forEach(function (t) {
        var tr = document.createElement('tr');
        cell(tr, t.stress + ' stress');
        cell(tr, t.cat);
        cell(tr, t.k + ' / ' + t.n + ' (' + (t.n ? (100 * t.k / t.n).toFixed(1) : '0') + '%)', 'num');
        cell(tr, (100 * t.K / t.N).toFixed(1) + '%', 'num');
        cell(tr, t.fold.toFixed(2), 'num');
        cell(tr, fmtP(t.q), 'num');
        var b = document.createElement('span');
        var sig = t.q < 0.05 && t.fold > 1;
        b.className = 'badge ' + (sig ? 'badge-high' : 'badge-none');
        b.textContent = sig ? 'Enriched' : 'Not enriched';
        cell(tr, b);
        tbody.appendChild(tr);
      });
    });
  }

  function renderSource() {
    var key = state.source;
    return load(SOURCES[key].file).then(function (db) {
      var res = state.results[key] || (state.results[key] = enrich(db, state.list, state.bg));
      var sig = res.rows.filter(function (r) { return r.q < 0.05; });
      $('enrich-summary').textContent = SOURCES[key].name + ': ' + fmt(res.n) + ' of your genes are annotated (background ' +
        fmt(res.N) + ' annotated genes). ' + fmt(res.rows.length) + ' terms tested, ' + fmt(sig.length) + ' significant (FDR < 0.05).';
      renderDots(res, SOURCES[key].name);
      var tbody = $('enrich-table').tBodies[0];
      tbody.innerHTML = '';
      res.rows.filter(function (r) { return r.p < 0.05; }).slice(0, 200).forEach(function (r) {
        var tr = document.createElement('tr');
        cell(tr, r.id, 'nowrap mono');
        cell(tr, r.name || r.id);
        cell(tr, r.k + ' / ' + r.n, 'num');
        cell(tr, r.K + ' / ' + r.N, 'num');
        cell(tr, r.fold.toFixed(2), 'num');
        cell(tr, fmtP(r.p), 'num');
        cell(tr, fmtP(r.q), 'num');
        var det = document.createElement('details');
        var sum = document.createElement('summary');
        sum.textContent = r.k + ' genes';
        det.appendChild(sum);
        var p = document.createElement('div');
        p.className = 'mono small';
        p.textContent = r.genes.join(', ');
        det.appendChild(p);
        cell(tr, det);
        if (r.q < 0.05) tr.className = 'sig-row';
        tbody.appendChild(tr);
      });
      $('enrich-none').hidden = res.rows.some(function (r) { return r.p < 0.05; });
    });
  }

  function renderDots(res, title) {
    var box = $('dot-plot');
    box.innerHTML = '';
    var rows = res.rows.filter(function (r) { return r.q < 0.05; }).slice(0, DOTS);
    if (!rows.length) { box.hidden = true; return; }
    box.hidden = false;
    var W = Math.max(660, box.clientWidth), rowH = 26, top = 70, labelW = Math.min(380, W * 0.45);
    var plotX = labelW + 16, plotW = W - plotX - 150;
    var H = top + rows.length * rowH + 50;
    // x axis: gene ratio (genes with the term / annotated genes in the list).
    var maxRatio = Math.max.apply(null, rows.map(function (r) { return r.k / r.n; }));
    var step = maxRatio <= 0.25 ? 0.05 : maxRatio <= 0.6 ? 0.1 : 0.2;
    var xmax = Math.min(1, Math.ceil(maxRatio * 1.08 / step) * step);
    var X = function (v) { return plotX + v / xmax * plotW; };
    var maxK = Math.max.apply(null, rows.map(function (r) { return r.k; }));
    var R = function (k) { return 3 + 8 * Math.sqrt(k / maxK); };
    var maxL = Math.max.apply(null, rows.map(function (r) { return -Math.log10(Math.max(r.q, 1e-300)); }));
    var color = function (q) {
      var t = -Math.log10(Math.max(q, 1e-300)) / (maxL || 1);
      return FDR_COLORS[Math.min(FDR_COLORS.length - 1, Math.floor(t * FDR_COLORS.length))];
    };
    var svg = el('svg', { xmlns: NS, width: W, height: H, viewBox: '0 0 ' + W + ' ' + H, 'font-family': FONT, role: 'img',
      'aria-label': title + ' enrichment dot plot' });
    el('rect', { x: 0, y: 0, width: W, height: H, fill: '#ffffff' }, svg);
    txt(svg, 10, 22, title + ': top ' + rows.length + ' enriched terms', { 'font-size': 15, 'font-weight': 700, fill: '#1f2328' });
    txt(svg, 10, 40, fmt(res.n) + ' annotated genes in the list; FDR < 0.05; sorted by p-value.', { 'font-size': 12, fill: '#57606a' });
    for (var i = 0; i * step <= xmax + 1e-9; i++) {
      var v = +(i * step).toFixed(2);
      el('line', { x1: X(v), x2: X(v), y1: top - 8, y2: top + rows.length * rowH, stroke: v === 0 ? '#8c959f' : '#e3e8e5' }, svg);
      txt(svg, X(v), top + rows.length * rowH + 16, String(v), { 'font-size': 11, fill: '#57606a', 'text-anchor': 'middle' });
    }
    txt(svg, plotX + plotW / 2, top + rows.length * rowH + 36, 'Gene ratio (genes with the term / annotated genes in your list)', { 'font-size': 12, 'font-weight': 600, fill: '#1f2328', 'text-anchor': 'middle' });
    rows.forEach(function (r, i) {
      var cy = top + i * rowH + rowH / 2;
      txt(svg, labelW, cy + 4, fit(r.name || r.id, labelW - 10, '12px ' + FONT), { 'font-size': 12, fill: '#1f2328', 'text-anchor': 'end' });
      el('line', { x1: plotX, x2: X(r.k / r.n), y1: cy, y2: cy, stroke: '#d0d7d3' }, svg);
      var dot = el('circle', { cx: X(r.k / r.n), cy: cy, r: R(r.k), fill: color(r.q), stroke: '#5c2a3a', 'stroke-width': 0.6 }, svg);
      el('title', {}, dot).textContent = (r.name || r.id) + ': ' + r.k + ' genes, ' + r.fold.toFixed(1) + '× enriched, FDR ' + fmtP(r.q);
    });
    // Legends
    var lx = W - 130, ly = top;
    txt(svg, lx, ly, 'Genes', { 'font-size': 11.5, 'font-weight': 600, fill: '#1f2328' });
    [maxK, Math.max(1, Math.round(maxK / 2)), Math.max(1, Math.round(maxK / 6))].forEach(function (k, i) {
      var cy = ly + 18 + i * 24;
      el('circle', { cx: lx + 10, cy: cy, r: R(k), fill: '#ffffff', stroke: '#57606a' }, svg);
      txt(svg, lx + 28, cy + 4, String(k), { 'font-size': 11, fill: '#57606a' });
    });
    ly += 100;
    txt(svg, lx, ly, '−log10 FDR', { 'font-size': 11.5, 'font-weight': 600, fill: '#1f2328' });
    FDR_COLORS.forEach(function (c, i) {
      el('rect', { x: lx + i * 18, y: ly + 8, width: 18, height: 10, fill: c }, svg);
    });
    txt(svg, lx, ly + 32, '1.3', { 'font-size': 10.5, fill: '#57606a' });
    txt(svg, lx + 90, ly + 32, maxL.toFixed(1), { 'font-size': 10.5, fill: '#57606a', 'text-anchor': 'end' });
    box.appendChild(svg);
  }

  function enrichCsv() {
    var res = state.results[state.source];
    if (!res) return;
    var rows = [['term_id', 'term_name', 'genes_in_list', 'annotated_list_genes', 'genes_in_background', 'background_genes', 'fold_enrichment', 'p_value', 'FDR', 'genes']];
    res.rows.forEach(function (r) { rows.push([r.id, r.name, r.k, r.n, r.K, r.N, r.fold.toFixed(3), r.p, r.q, r.genes.join(' ')]); });
    SoybioExport.saveCsv(rows, 'enrichment_' + state.source + '.csv');
  }

  function stressCsv() {
    var t = state.stressTests;
    if (!t) return;
    var rows = [['stress', 'category', 'genes_in_list', 'list_size', 'genes_in_background', 'background_size', 'fold_enrichment', 'p_value', 'FDR', 'genes']];
    t.forEach(function (r) { rows.push([r.stress, r.cat, r.k, r.n, r.K, r.N, r.fold.toFixed(3), r.p, r.q, r.genes.join(' ')]); });
    SoybioExport.saveCsv(rows, 'stress_enrichment.csv');
  }

  // ---------- Single gene analysis ----------

  function coexpFile(gene) {
    var m = /^Glyma\.(\d{2}|U)G?(\d+)$/.exec(gene);
    var b = Math.floor(parseInt(m[2], 10) / meta.bucket);
    return SoybioData.BASE + 'coexp/' + m[1] + '_' + (b < 10 ? '0' + b : b) + '.json';
  }

  function runSingle(item) {
    var g = item.data;
    state.mode = 'single';
    $('list').hidden = true;
    $('single').hidden = false;
    $('results').hidden = false;
    setStatus('Analysing ' + g.gene + '…');
    var srcs = ['go_bp', 'go_mf', 'go_cc', 'kegg', 'pfam'];
    Promise.all(srcs.map(function (k) { return load(SOURCES[k].file); }).concat([
      SoybioData.getJSON(coexpFile(g.gene)).catch(function () { return {}; })
    ])).then(function (r) {
      var dbs = {};
      srcs.forEach(function (k, i) { dbs[k] = r[i]; });
      var partners = r[srcs.length][g.gene] || [];
      renderSummary(g, item, dbs);
      renderGo(g, dbs);
      renderResponse(g);
      return renderPartners(g, partners, dbs);
    }).then(function () {
      setStatus('Done.');
      $('results').scrollIntoView({ behavior: 'smooth', block: 'start' });
    }).catch(function () { setStatus('The data could not be loaded. Please reload the page.', true); });
  }

  function termNames(db, gene) {
    return (db.genes[gene] || []).map(function (t) { return db.terms[t]; });
  }

  function renderSummary(g, item, dbs) {
    var box = $('gene-summary');
    box.innerHTML = '';
    var h = document.createElement('h3');
    h.className = 'mono';
    h.textContent = g.gene;
    box.appendChild(h);
    var dl = document.createElement('dl');
    dl.className = 'facts';
    var add = function (k, v) {
      var dt = document.createElement('dt'); dt.textContent = k;
      var dd = document.createElement('dd');
      if (v instanceof Node) dd.appendChild(v); else dd.textContent = v;
      dl.appendChild(dt); dl.appendChild(dd);
    };
    if (item.from && item.from !== 'Wm82.a4') add('Entered ID', item.query + ' (' + item.from + ')');
    add('ZH13 gene', g.zh || 'No matching ZH13 gene');
    add('Description', g.desc || 'No description');
    var pf = termNames(dbs.pfam, g.gene).map(function (t) { return t[0]; });
    add('Pfam domains', pf.length ? pf.join(', ') : 'None annotated');
    var ko = dbs.kegg.ko[g.gene];
    var pws = termNames(dbs.kegg, g.gene).filter(function (t) { return !t[2]; }).map(function (t) { return t[1] + ' (' + t[0] + ')'; });
    add('KEGG orthology', ko ? ko.replace(/,/g, ', ') : 'None annotated');
    add('KEGG pathways', pws.length ? pws.join('; ') : 'None annotated');
    add('Expression', link('/expression.html?genes=' + encodeURIComponent(g.gene), 'Bar plot under salt, drought and heat stress'));
    box.appendChild(dl);
  }

  // GO terms of the gene, most specific first (fewest genes in the genome).
  function renderGo(g, dbs) {
    var box = $('gene-go');
    box.innerHTML = '';
    ['go_bp', 'go_mf', 'go_cc'].forEach(function (k) {
      var db = dbs[k];
      var idx = db.genes[g.gene] || [];
      var h = document.createElement('h4');
      h.textContent = SOURCES[k].name.replace('GO ', '') + ' (' + idx.length + ' terms)';
      box.appendChild(h);
      if (!idx.length) {
        var p = document.createElement('p');
        p.className = 'muted small';
        p.textContent = 'No GO terms annotated.';
        box.appendChild(p);
        return;
      }
      if (!db._size) {
        db._size = new Int32Array(db.terms.length);
        Object.keys(db.genes).forEach(function (gg) { db.genes[gg].forEach(function (t) { db._size[t]++; }); });
      }
      var sorted = idx.slice().sort(function (a, b) { return db._size[a] - db._size[b]; });
      var ul = document.createElement('ul');
      ul.className = 'term-list';
      sorted.forEach(function (t, i) {
        var li = document.createElement('li');
        li.textContent = db.terms[t][1] + ' ';
        var s = document.createElement('span');
        s.className = 'muted mono small';
        s.textContent = db.terms[t][0];
        li.appendChild(s);
        if (i >= 8) li.hidden = true;
        ul.appendChild(li);
      });
      box.appendChild(ul);
      if (sorted.length > 8) {
        var more = document.createElement('button');
        more.type = 'button';
        more.className = 'btn btn-outline btn-small';
        more.textContent = 'Show all ' + sorted.length + ' terms';
        more.addEventListener('click', function () {
          Array.prototype.forEach.call(ul.children, function (li) { li.hidden = false; });
          more.remove();
        });
        box.appendChild(more);
      }
    });
  }

  function renderResponse(g) {
    var tbody = $('response-table').tBodies[0];
    tbody.innerHTML = '';
    STRESSES.forEach(function (s) {
      var up = 0, down = 0, n = 0;
      meta.contrasts.forEach(function (c, i) {
        if (c.stress !== s) return;
        n++;
        var f = g.f[i], p = g.p[i];
        if (f !== null && p !== null && p < 0.05 && Math.abs(f) >= 1) { if (f > 0) up++; else down++; }
      });
      var tr = document.createElement('tr');
      cell(tr, s + ' stress');
      cell(tr, String(n), 'num');
      cell(tr, String(up), 'num');
      cell(tr, String(down), 'num');
      var label = up + down === 0 ? 'No response' : up >= 2 * down ? 'Mainly up' : down >= 2 * up ? 'Mainly down' : 'Mixed';
      var b = document.createElement('span');
      b.className = 'badge ' + ({ 'Mainly up': 'badge-up', 'Mainly down': 'badge-down' }[label] || 'badge-none');
      b.textContent = label;
      cell(tr, b);
      tbody.appendChild(tr);
    });
  }

  function renderPartners(g, partners, dbs) {
    var tbody = $('partner-table').tBodies[0];
    tbody.innerHTML = '';
    $('partners-none').hidden = partners.length > 0;
    $('partners-wrap').hidden = !partners.length;
    $('predict').hidden = !partners.length;
    if (!partners.length) return Promise.resolve();
    var ids = partners.map(function (p) { return p[0]; });
    $('partner-heatmap').href = '/expression.html?genes=' + encodeURIComponent([g.gene].concat(ids).join(','));
    state.single = { gene: g.gene, partners: partners };
    renderPrediction(g, ids, dbs);
    return Promise.all(ids.map(SoybioData.getGene)).then(function (data) {
      partners.forEach(function (p, i) {
        var d = data[i];
        var tr = document.createElement('tr');
        cell(tr, String(i + 1), 'num');
        cell(tr, link('/analysis.html?genes=' + encodeURIComponent(p[0]), p[0], 'mono'), 'nowrap');
        cell(tr, p[1].toFixed(2), 'num');
        var td = cell(tr, d ? shortDesc(d.desc) || '–' : '–', 'desc-cell');
        if (d && d.desc) td.title = d.desc;
        tbody.appendChild(tr);
        p.desc = d ? d.desc : '';
      });
    });
  }

  // Function prediction: terms over-represented among the co-expressed genes.
  function renderPrediction(g, ids, dbs) {
    var box = $('predict-body');
    box.innerHTML = '';
    var any = false;
    ['go_bp', 'go_mf', 'kegg', 'pfam'].forEach(function (k) {
      var res = enrich(dbs[k], ids, null);
      var sig = res.rows.filter(function (r) { return r.q < 0.05; }).slice(0, 5);
      if (!sig.length) return;
      any = true;
      var h = document.createElement('h4');
      h.textContent = SOURCES[k].name;
      box.appendChild(h);
      var ul = document.createElement('ul');
      ul.className = 'term-list';
      sig.forEach(function (r) {
        var li = document.createElement('li');
        li.textContent = (r.name || r.id) + ' ';
        var s = document.createElement('span');
        s.className = 'muted small';
        s.textContent = '(' + r.k + ' of ' + res.n + ' annotated partners, ' + r.fold.toFixed(1) + '× enriched, FDR ' + fmtP(r.q) + ')';
        li.appendChild(s);
        ul.appendChild(li);
      });
      box.appendChild(ul);
    });
    if (!any) {
      var p = document.createElement('p');
      p.className = 'muted';
      p.textContent = 'No GO, KEGG or Pfam term is significantly over-represented among the co-expressed genes (FDR < 0.05), so no function can be predicted from co-expression.';
      box.appendChild(p);
    }
  }

  function partnersCsv() {
    var s = state.single;
    if (!s) return;
    var rows = [['gene', 'rank', 'coexpressed_gene', 'pearson_r', 'description']];
    s.partners.forEach(function (p, i) { rows.push([s.gene, i + 1, p[0], p[1], p.desc || '']); });
    SoybioExport.saveCsv(rows, s.gene + '_coexpressed_genes.csv');
  }

  // ---------- Main ----------

  function run() {
    var tokens = SoybioData.parseTokens($('genes').value);
    if (!tokens.length) { setStatus('Please enter at least one gene ID.', true); return; }
    if (tokens.length > MAX_LIST) { setStatus('Please enter up to ' + fmt(MAX_LIST) + ' genes.', true); return; }
    if (window.history && history.replaceState && tokens.length <= 60) {
      history.replaceState(null, '', '?genes=' + encodeURIComponent(tokens.join(',')));
    }
    setStatus('Loading…');
    if (tokens.length === 1) {
      SoybioData.resolve(tokens).then(function (items) {
        if (!items[0].data) { setStatus(tokens[0] + ': ' + items[0].note + '.', true); return; }
        runSingle(items[0]);
      });
    } else {
      runList(tokens);
    }
  }

  SoybioData.loadMeta().then(function (m) {
    meta = m;
    $('run').disabled = false;
    setStatus('Ready.');
    var q = /[?&]genes=([^&]*)/.exec(location.search);
    if (q) {
      $('genes').value = decodeURIComponent(q[1].replace(/\+/g, ' ')).split(/[\s,;]+/).join('\n');
      run();
    }
  }).catch(function () { setStatus('The data could not be loaded. Please reload the page.', true); });

  $('analysis-form').addEventListener('submit', function (e) { e.preventDefault(); run(); });
  $('ex-single').addEventListener('click', function () { $('genes').value = 'Glyma.14G100000'; run(); });
  $('ex-list').addEventListener('click', function () {
    Promise.all([load('stress.json')]).then(function (r) {
      // Example list: genes up-regulated by heat in both heat comparisons.
      var st = r[0].genes;
      var list = Object.keys(st).filter(function (g) { return st[g][4] >= 2; }).sort();
      $('genes').value = list.join('\n');
      run();
    });
  });
  $('clear').addEventListener('click', function () {
    $('genes').value = '';
    $('results').hidden = true;
    $('genes').focus();
  });
  Array.prototype.forEach.call(document.querySelectorAll('#source-filter button'), function (b) {
    b.addEventListener('click', function () {
      state.source = b.getAttribute('data-source');
      Array.prototype.forEach.call(document.querySelectorAll('#source-filter button'), function (o) {
        o.setAttribute('aria-pressed', o === b ? 'true' : 'false');
      });
      if (state.list) renderSource();
    });
  });
  $('enrich-csv').addEventListener('click', enrichCsv);
  $('stress-csv').addEventListener('click', stressCsv);
  $('dot-svg').addEventListener('click', function () { SoybioExport.saveSvg($('dot-plot').querySelector('svg'), 'enrichment_' + state.source + '.svg'); });
  $('dot-png').addEventListener('click', function () { SoybioExport.savePng($('dot-plot').querySelector('svg'), 'enrichment_' + state.source + '.png'); });
  $('partner-csv').addEventListener('click', partnersCsv);
})();

/* Shared loading of the abiotic stress RNA-seq data (built by tools/build_expression.py).
   Used by the Expression and Reference Genes pages. Exposes window.SoybioData. */
(function () {
  'use strict';

  var BASE = '/data/abiotic-rnaseq/';
  var meta = null;
  var idmap = null;
  var chunkCache = {};

  function getJSON(url) {
    return fetch(url).then(function (r) {
      if (!r.ok) throw new Error(url + ': HTTP ' + r.status);
      return r.json();
    });
  }

  function loadMeta() {
    return getJSON(BASE + 'meta.json').then(function (m) { meta = m; return m; });
  }

  function chunkKey(gene) {
    var m = /^Glyma\.(\d{2}|U)G?(\d+)$/.exec(gene);
    if (!m) return null;
    var b = Math.floor(parseInt(m[2], 10) / meta.bucket);
    return m[1] + '_' + (b < 10 ? '0' + b : b);
  }

  function loadChunk(key) {
    if (!chunkCache[key]) {
      chunkCache[key] = getJSON(BASE + 'chunks/' + key + '.json').catch(function () { return {}; });
    }
    return chunkCache[key];
  }

  // Expression data of one Wm82.a4 gene: {gene, zh, desc, t: TPM per sample, f: log2FC and p: padj per contrast}.
  function getGene(gene) {
    var key = chunkKey(gene);
    if (!key) return Promise.resolve(null);
    return loadChunk(key).then(function (c) {
      var g = c[gene];
      return g ? { gene: gene, zh: g.z, desc: g.d, t: g.t, f: g.f, p: g.p } : null;
    });
  }

  function loadIdmap() {
    if (!idmap) idmap = getJSON(BASE + 'idmap.json');
    return idmap;
  }

  // Normalize an entered ID: remove prefixes and transcript/protein suffixes.
  function canonical(tok) {
    var s = tok.replace(/^["']+|["']+$/g, '').replace(/^glyma\.Wm82\.gnm\d+\.ann\d+\./i, '');
    var up = s.toUpperCase()
      .replace(/\.WM82\.A\d+\.V\d+$/, '')
      .replace(/\.P$/, '')
      .replace(/\.(M|T)?\d+$/, '');
    if (/^SOYZH13_/.test(up)) return { type: 'zh', id: 'SoyZH13_' + up.slice(8) };
    if (/^GLYMA\./.test(up)) return { type: 'wm', id: 'Glyma.' + up.slice(6) };
    return { type: 'unknown', id: s };
  }

  // Look up entered IDs (Wm82.a4, ZH13 or Wm82.a2.v1).
  // Returns [{query, data, from, note}] in input order; data is null when not found.
  function resolve(tokens) {
    var items = tokens.map(function (q) { return { query: q, c: canonical(q), data: null, from: '', note: '' }; });
    var first = items.map(function (it) {
      if (it.c.type !== 'wm') return null;
      return getGene(it.c.id).then(function (d) { if (d) { it.data = d; it.from = 'Wm82.a4'; } });
    });
    return Promise.all(first).then(function () {
      var needMap = items.some(function (it) { return !it.data && it.c.type !== 'unknown'; });
      return needMap ? loadIdmap() : null;
    }).then(function (map) {
      return Promise.all(items.map(function (it) {
        if (it.data) return null;
        if (it.c.type === 'unknown') { it.note = 'not a ZH13 or Glyma gene ID'; return null; }
        var target = it.c.type === 'zh' ? map.zh[it.c.id] : map.a2[it.c.id];
        if (!target) {
          it.note = it.c.type === 'zh' ? 'no matching Wm82.a4 gene' : 'not found';
          return null;
        }
        return getGene(target).then(function (d) {
          if (d) { it.data = d; it.from = it.c.type === 'zh' ? 'ZH13' : 'Wm82.a2.v1'; }
          else it.note = 'matching gene ' + target + ' has no expression data';
        });
      }));
    }).then(function () { return items; });
  }

  // Split pasted text into unique IDs (case-insensitive), keeping input order.
  function parseTokens(text) {
    var seen = {};
    var out = [];
    text.split(/[\s,;]+/).forEach(function (t) {
      t = t.trim();
      if (!t) return;
      var k = t.toUpperCase();
      if (seen[k]) return;
      seen[k] = true;
      out.push(t);
    });
    return out;
  }

  window.SoybioData = {
    BASE: BASE,
    getJSON: getJSON,
    loadMeta: loadMeta,
    getGene: getGene,
    resolve: resolve,
    parseTokens: parseTokens
  };
})();

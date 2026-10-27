/* Gene ID Converter: ZH13 <-> Wm82.a4 (and Wm82.a2.v1).
   Loads data/zh13-wm82a4/lookup.json (built by tools/build_data.py) and
   looks up IDs in the browser. */
(function () {
  'use strict';

  var DATA_URL = '/data/zh13-wm82a4/lookup.json';
  var MAX_IDS = 20000;
  var MAX_ROWS_SHOWN = 2000;
  var EXAMPLE = [
    'SoyZH13_01G000024',
    'SoyZH13_01G029800',
    'SoyZH13_01G000001',
    'Glyma.03G110100',
    'Glyma.01G000137',
    'Glyma.01G000322.1',
    'Glyma.06G220000'
  ];

  // Columns of the results table. "basic" columns are shown by default;
  // all columns are shown with "Show all columns" and are always in the CSV.
  var COLUMNS = [
    { key: 'query', label: 'Query', basic: true, cls: 'nowrap mono' },
    { key: 'type', label: 'Input type', cls: 'nowrap' },
    { key: 'zh', label: 'ZH13 gene', basic: true, cls: 'nowrap mono' },
    { key: 'wm', label: 'Wm82.a4 gene', basic: true, cls: 'nowrap mono' },
    { key: 'a2', label: 'Wm82.a2.v1 gene', basic: true, cls: 'nowrap mono' },
    { key: 'conf', label: 'Confidence', basic: true, badge: true },
    { key: 'pid', label: 'Protein identity (%)', basic: true, cls: 'num' },
    { key: 'cds', label: 'CDS lifted (%)', cls: 'num' },
    { key: 'zh_loc', label: 'ZH13 location', cls: 'nowrap' },
    { key: 'wm_loc', label: 'Wm82.a4 location', cls: 'nowrap' },
    { key: 'block', label: 'Synteny block', cls: 'nowrap' },
    { key: 'liftover', label: 'Liftover Wm82.a4 gene', cls: 'nowrap mono' },
    { key: 'evidence', label: 'Evidence', cls: 'wide' },
    { key: 'soyod', label: 'SoyOD converter answer', cls: 'wide' },
    { key: 'agree', label: 'Agreement with SoyOD', cls: 'nowrap' },
    { key: 'desc', label: 'Wm82.a4 description', basic: true, cls: 'desc' }
  ];

  var TYPE_LABEL = { zh: 'ZH13', wm: 'Wm82.a4', a2: 'Wm82.a2.v1', unknown: 'Unknown' };
  var NOT_FOUND = 'Not found';

  var $ = function (id) { return document.getElementById(id); };
  var form = $('converter-form');
  var input = $('ids');
  var statusEl = $('status');
  var convertBtn = $('convert');
  var resultsEl = $('results');
  var allColumnsBox = $('all-columns');

  var db = null;
  var results = [];

  // ---------- Data loading ----------

  // Wraps one table of lookup.json so values can be read by field name.
  function Table(t) {
    var pos = {};
    t.fields.forEach(function (f, i) { pos[f] = i; });
    this.rows = t.rows;
    this.get = function (row, field) {
      var v = row[pos[field]];
      if (v === null || v === undefined) return '';
      if (t.dicts[field]) return t.dicts[field][v];
      // In these columns, 1 means "same as the Wm82.a4 gene of this row".
      if (v === 1 && (field === 'liftover' || field === 'soyod')) return row[pos.wm];
      return v;
    };
  }

  function buildDb(json) {
    var zh = new Table(json.zh);
    var wm = new Table(json.wm);
    var zhIndex = new Map();
    var wmIndex = new Map();
    var a2Index = new Map();
    zh.rows.forEach(function (r, i) { zhIndex.set(zh.get(r, 'zh').toUpperCase(), i); });
    wm.rows.forEach(function (r, i) {
      wmIndex.set(wm.get(r, 'wm').toUpperCase(), i);
      var a2 = wm.get(r, 'a2');
      if (a2) a2Index.set(a2.toUpperCase(), i);
    });
    return { zh: zh, wm: wm, zhIndex: zhIndex, wmIndex: wmIndex, a2Index: a2Index };
  }

  function setStatus(text, isError) {
    statusEl.textContent = text;
    statusEl.classList.toggle('error', !!isError);
  }

  fetch(DATA_URL)
    .then(function (res) {
      if (!res.ok) throw new Error('HTTP ' + res.status);
      return res.json();
    })
    .then(function (json) {
      db = buildDb(json);
      convertBtn.disabled = false;
      setStatus('Ready. ' + fmt(db.zh.rows.length) + ' ZH13 genes and ' +
        fmt(db.wm.rows.length) + ' Wm82.a4 genes loaded.');
    })
    .catch(function () {
      setStatus('The conversion table could not be loaded. Please reload the page, ' +
        'or download the Excel file from the Data page.', true);
    });

  // ---------- Lookup ----------

  function fmt(n) { return n.toLocaleString('en-US'); }

  function parseIds(text) {
    var seen = new Set();
    var ids = [];
    var duplicates = 0;
    text.split(/[\s,;]+/).forEach(function (raw) {
      var id = raw.replace(/^["']+|["']+$/g, '');
      if (!id) return;
      var key = id.toUpperCase();
      if (seen.has(key)) { duplicates++; return; }
      seen.add(key);
      ids.push(id);
    });
    return { ids: ids, duplicates: duplicates };
  }

  // Gene ID forms to try: the ID as entered, then without transcript/protein suffixes.
  function candidates(id) {
    var up = id.toUpperCase();
    var base = up
      .replace(/\.WM82\.A\d+\.V\d+$/, '')
      .replace(/\.P$/, '')
      .replace(/\.(M|T)?\d+$/, '');
    return base !== up ? [up, base] : [up];
  }

  function find(index, id) {
    var keys = candidates(id);
    for (var i = 0; i < keys.length; i++) {
      if (index.has(keys[i])) return index.get(keys[i]);
    }
    return -1;
  }

  function lookup(id, mode) {
    var type = mode;
    if (mode === 'auto') {
      if (/^SoyZH13_/i.test(id)) type = 'zh';
      else if (/^Glyma\./i.test(id)) type = 'wm';
      else type = 'unknown';
    }

    var zi = -1;
    var wi = -1;
    if (type === 'zh') {
      zi = find(db.zhIndex, id);
      if (zi >= 0) {
        var wmId = db.zh.get(db.zh.rows[zi], 'wm');
        if (wmId && db.wmIndex.has(wmId.toUpperCase())) wi = db.wmIndex.get(wmId.toUpperCase());
      }
    } else if (type === 'wm' || type === 'a2') {
      wi = find(type === 'wm' ? db.wmIndex : db.a2Index, id);
      if (wi < 0 && mode === 'auto') {
        wi = find(db.a2Index, id);
        if (wi >= 0) type = 'a2';
      }
      if (wi >= 0) {
        var zhId = db.wm.get(db.wm.rows[wi], 'zh');
        if (zhId && db.zhIndex.has(zhId.toUpperCase())) zi = db.zhIndex.get(zhId.toUpperCase());
      }
    }
    return makeRecord(id, type, zi, wi);
  }

  function makeRecord(query, type, zi, wi) {
    var rec = { query: query, type: TYPE_LABEL[type] };
    COLUMNS.forEach(function (c) { if (!(c.key in rec)) rec[c.key] = ''; });
    if (zi < 0 && wi < 0) {
      rec.conf = NOT_FOUND;
      return rec;
    }
    var Z = db.zh;
    var W = db.wm;
    if (zi >= 0) {
      var z = Z.rows[zi];
      ['zh', 'wm', 'conf', 'evidence', 'liftover', 'pid', 'cds', 'block', 'zh_loc', 'soyod', 'agree']
        .forEach(function (f) { rec[f] = Z.get(z, f); });
    }
    if (wi >= 0) {
      var w = W.rows[wi];
      ['wm', 'zh', 'conf', 'a2', 'wm_loc', 'desc'].forEach(function (f) { rec[f] = W.get(w, f); });
    }
    return rec;
  }

  // ---------- Output ----------

  function badgeClass(conf) {
    if (conf === 'High') return 'badge badge-high';
    if (conf === 'Medium') return 'badge badge-medium';
    if (conf === NOT_FOUND) return 'badge badge-notfound';
    return 'badge badge-none';
  }

  function visibleColumns() {
    return allColumnsBox.checked ? COLUMNS : COLUMNS.filter(function (c) { return c.basic; });
  }

  function render() {
    var cols = visibleColumns();
    var table = $('results-table');
    var thead = table.tHead;
    var tbody = table.tBodies[0];

    var headRow = document.createElement('tr');
    cols.forEach(function (c) {
      var th = document.createElement('th');
      th.scope = 'col';
      th.textContent = c.label;
      if (c.cls && c.cls.indexOf('num') >= 0) th.className = 'num';
      headRow.appendChild(th);
    });
    thead.innerHTML = '';
    thead.appendChild(headRow);

    var frag = document.createDocumentFragment();
    results.slice(0, MAX_ROWS_SHOWN).forEach(function (rec) {
      var tr = document.createElement('tr');
      cols.forEach(function (c) {
        var td = document.createElement('td');
        if (c.cls) td.className = c.cls;
        var value = rec[c.key];
        if (c.badge && value) {
          var span = document.createElement('span');
          span.className = badgeClass(value);
          span.textContent = value;
          td.appendChild(span);
        } else {
          td.textContent = value === '' ? '–' : String(value);
          if (value === '') td.classList.add('muted');
        }
        tr.appendChild(td);
      });
      frag.appendChild(tr);
    });
    tbody.innerHTML = '';
    tbody.appendChild(frag);

    $('results-note').textContent = results.length > MAX_ROWS_SHOWN
      ? 'Showing the first ' + fmt(MAX_ROWS_SHOWN) + ' of ' + fmt(results.length) +
        ' rows. Download the CSV file to get all rows.'
      : '';
  }

  function summarize(duplicates) {
    var paired = 0, none = 0, missing = 0;
    results.forEach(function (r) {
      if (r.conf === NOT_FOUND) missing++;
      else if (r.zh && r.wm) paired++;
      else none++;
    });
    var parts = [
      fmt(paired) + ' with a matching gene',
      fmt(none) + ' with no correspondence',
      fmt(missing) + ' not found in the table'
    ];
    var text = fmt(results.length) + (results.length === 1 ? ' ID: ' : ' IDs: ') + parts.join(', ') + '.';
    if (duplicates) text += ' ' + fmt(duplicates) + ' repeated ' + (duplicates === 1 ? 'ID was' : 'IDs were') + ' skipped.';
    $('summary').textContent = text;
  }

  function convert() {
    if (!db) return;
    var parsed = parseIds(input.value);
    if (!parsed.ids.length) {
      setStatus('Please enter at least one gene ID.', true);
      resultsEl.hidden = true;
      return;
    }
    if (parsed.ids.length > MAX_IDS) {
      setStatus('Too many IDs (' + fmt(parsed.ids.length) + '). Please enter up to ' +
        fmt(MAX_IDS) + ' IDs at a time.', true);
      return;
    }
    var mode = form.elements.idtype.value;
    results = parsed.ids.map(function (id) { return lookup(id, mode); });
    summarize(parsed.duplicates);
    render();
    resultsEl.hidden = false;
    setStatus('Done.');
    resultsEl.scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  function toDelimited(sep, cols) {
    var quote = function (v) {
      var s = String(v);
      if (sep === '\t') return s.replace(/[\t\r\n]+/g, ' ');
      return /[",\r\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
    };
    var lines = [cols.map(function (c) { return quote(c.label); }).join(sep)];
    results.forEach(function (rec) {
      lines.push(cols.map(function (c) { return quote(rec[c.key]); }).join(sep));
    });
    return lines.join('\r\n');
  }

  function downloadCsv() {
    // The byte order mark lets Excel read the file as UTF-8.
    var blob = new Blob(['﻿' + toDelimited(',', COLUMNS)], { type: 'text/csv;charset=utf-8' });
    var url = URL.createObjectURL(blob);
    var a = document.createElement('a');
    a.href = url;
    a.download = 'gene_id_conversion_results.csv';
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
  }

  function copyTable() {
    var text = toDelimited('\t', visibleColumns());
    var btn = $('copy');
    var done = function () {
      btn.textContent = 'Copied';
      setTimeout(function () { btn.textContent = 'Copy table'; }, 1500);
    };
    if (navigator.clipboard && window.isSecureContext) {
      navigator.clipboard.writeText(text).then(done);
    } else {
      var ta = document.createElement('textarea');
      ta.value = text;
      document.body.appendChild(ta);
      ta.select();
      document.execCommand('copy');
      ta.remove();
      done();
    }
  }

  // ---------- Events ----------

  form.addEventListener('submit', function (e) {
    e.preventDefault();
    convert();
  });
  input.addEventListener('keydown', function (e) {
    if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) convert();
  });
  $('example').addEventListener('click', function () {
    input.value = EXAMPLE.join('\n');
    form.elements.idtype.value = 'auto';
    if (db) convert();
  });
  $('clear').addEventListener('click', function () {
    input.value = '';
    results = [];
    resultsEl.hidden = true;
    if (db) setStatus('Ready.');
    input.focus();
  });
  allColumnsBox.addEventListener('change', render);
  $('download').addEventListener('click', downloadCsv);
  $('copy').addEventListener('click', copyTable);
})();

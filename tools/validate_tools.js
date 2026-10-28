// Runs the Tools page calculations (assets/js/seqtools.js) on test cases made by
// validate_tools.py and writes the results for comparison with Biopython/SciPy.
// Usage: node tools/validate_tools.js cases.json results.json
'use strict';
const fs = require('fs');
const path = require('path');

global.window = {};
eval(fs.readFileSync(path.join(__dirname, '..', 'assets', 'js', 'seqdata.js'), 'utf8'));
global.SoybioSeqData = window.SoybioSeqData;
const S = require(path.join(__dirname, '..', 'assets', 'js', 'seqtools.js'));
const D = window.SoybioSeqData;

const cases = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const out = {};

out.protein = cases.protein.map(p => {
  const r = S.proteinParams(p);
  return { mw: r.mw, pI: r.pI, charge7: r.charge7, gravy: r.gravy, instability: r.instability,
    aromaticity: r.aromaticity, ext: [r.extReduced, r.extCystines], hydro9: S.hydropathy(p, 9).map(h => h.y) };
});
out.tm = cases.tm.map(c => S.tmNN(c.seq, c));
out.translate = cases.dna.map(s => S.translate(s, false));
out.revcomp = cases.dna.map(s => S.revcomp(s));
out.gc = cases.dna.map(s => S.gcContent(s));
out.orfs = cases.orf.map(s => S.findOrfs(s, 30).map(o => [o.strand, o.start, o.end, o.protein]));
out.codons = cases.cds.map(s => S.codonUsage(s).rows.map(r => [r.codon, r.count, r.rscu]));
out.cuts = cases.restriction.map(s => {
  const res = {};
  D.enzymes.forEach(e => { res[e.name] = S.cutSites(s, e); });
  return res;
});
out.align = cases.align.map(c => {
  const r = S.align(c.a, c.b, c.opt);
  // Re-score the returned alignment to check the traceback.
  let sc = 0, gapA = 0, gapB = 0;
  for (let i = 0; i < r.a.length; i++) {
    const x = r.a[i], y = r.b[i];
    if (x === '-') { sc -= gapA ? c.opt.extend : c.opt.open; gapA = 1; gapB = 0; }
    else if (y === '-') { sc -= gapB ? c.opt.extend : c.opt.open; gapB = 1; gapA = 0; }
    else {
      gapA = gapB = 0;
      sc += c.opt.type === 'protein' ? D.blosum62[x][y] : (x === y ? c.opt.match : c.opt.mismatch);
    }
  }
  return { score: r.score, rescore: sc, degapA: r.a.replace(/-/g, ''), degapB: r.b.replace(/-/g, ''), startA: r.startA, startB: r.startB };
});
out.comp = cases.comp.map(c => S.complementarity(c[0], c[1]));
out.welch = cases.welch.map(c => { const w = S.welch(c[0], c[1]); return w ? w.p : null; });
out.curve = cases.curve.map(pts => { const r = S.standardCurve(pts); return [r.slope, r.intercept, r.r2, r.efficiency]; });
out.rel = cases.rel.map(c => S.relativeExpression(c.rows, c.control, c.et, c.er).samples.map(s => s.fold));
out.pcr = cases.pcr.map(c => S.pcrProducts(c.t, c.f, c.r).map(p => [p.start, p.end]));

fs.writeFileSync(process.argv[3], JSON.stringify(out));
console.log('results written');

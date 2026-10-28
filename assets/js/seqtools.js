/* Calculations for the Tools page (no page code here, so it can be tested in Node).
   Constants come from seqdata.js (generated from Biopython); formulas follow Biopython
   1.87 and are checked against it by tools/validate_tools.py. Exposes window.SoybioSeq. */
(function (root) {
  'use strict';

  var D = root.SoybioSeqData;
  var AA20 = 'ACDEFGHIKLMNPQRSTVWY';

  // ---------- Input ----------

  // Read pasted text (plain or FASTA; only the first record is used).
  function readSeq(text, type) {
    var lines = String(text || '').split(/\r?\n/);
    var name = '';
    var parts = [];
    var started = false;
    for (var i = 0; i < lines.length; i++) {
      var l = lines[i].trim();
      if (l.charAt(0) === '>') {
        if (started) break;
        name = l.slice(1).trim();
        started = true;
        continue;
      }
      if (l) { parts.push(l); started = true; }
    }
    var seq = parts.join('').replace(/[\s\d]/g, '').toUpperCase();
    var allowed;
    if (type === 'dna') {
      seq = seq.replace(/U/g, 'T');
      allowed = /^[ACGTN]*$/;
    } else {
      seq = seq.replace(/\*$/, '');
      allowed = new RegExp('^[' + AA20 + ']*$');
    }
    var bad = {};
    for (var j = 0; j < seq.length; j++) if (!allowed.test(seq.charAt(j))) bad[seq.charAt(j)] = true;
    return { name: name, seq: seq, invalid: Object.keys(bad) };
  }

  // ---------- DNA ----------

  var COMP = { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N' };
  function complement(s) { return s.replace(/[ACGTN]/g, function (c) { return COMP[c]; }); }
  function revcomp(s) { return complement(s).split('').reverse().join(''); }

  function gcContent(s) {
    var gc = 0, n = 0;
    for (var i = 0; i < s.length; i++) {
      var c = s.charAt(i);
      if (c === 'N') continue;
      n++;
      if (c === 'G' || c === 'C') gc++;
    }
    return n ? gc / n * 100 : 0;
  }

  function baseCounts(s) {
    var c = { A: 0, C: 0, G: 0, T: 0, N: 0 };
    for (var i = 0; i < s.length; i++) c[s.charAt(i)]++;
    return c;
  }

  // GC % in sliding windows; x is the window centre (1-based).
  function gcWindows(s, win, step) {
    var out = [];
    if (s.length < win) return out;
    for (var i = 0; i + win <= s.length; i += step) out.push({ x: i + win / 2 + 0.5, y: gcContent(s.substr(i, win)) });
    return out;
  }

  function translate(s, toStop) {
    var p = '';
    for (var i = 0; i + 3 <= s.length; i += 3) {
      var cod = s.substr(i, 3);
      var aa = D.stops.indexOf(cod) >= 0 ? '*' : (D.codons[cod] || 'X');
      if (aa === '*' && toStop) break;
      p += aa;
    }
    return p;
  }

  // Complete open reading frames (ATG ... stop) in all six frames.
  // Positions are 1-based on the entered (forward) strand; end includes the stop codon.
  function findOrfs(s, minAA) {
    var orfs = [];
    var L = s.length;
    [1, -1].forEach(function (strand) {
      var t = strand === 1 ? s : revcomp(s);
      for (var f = 0; f < 3; f++) {
        var start = -1;
        for (var i = f; i + 3 <= L; i += 3) {
          var cod = t.substr(i, 3);
          if (start < 0 && cod === 'ATG') start = i;
          if (start >= 0 && D.stops.indexOf(cod) >= 0) {
            var aaLen = (i - start) / 3;
            if (aaLen >= minAA) {
              var a = start + 1, b = i + 3;
              orfs.push({
                strand: strand, frame: strand * (f + 1),
                start: strand === 1 ? a : L - b + 1,
                end: strand === 1 ? b : L - a + 1,
                nt: b - a + 1, aa: aaLen,
                protein: translate(t.slice(start, i), false)
              });
            }
            start = -1;
          }
        }
      }
    });
    orfs.sort(function (x, y) { return y.aa - x.aa || x.start - y.start; });
    return orfs;
  }

  // Codon counts and RSCU (relative synonymous codon usage) of a coding sequence.
  function codonUsage(s) {
    var counts = {};
    var all = Object.keys(D.codons).concat(D.stops);
    all.forEach(function (c) { counts[c] = 0; });
    var other = 0;
    for (var i = 0; i + 3 <= s.length; i += 3) {
      var c = s.substr(i, 3);
      if (c in counts) counts[c]++; else other++;
    }
    var byAA = {};
    all.forEach(function (c) {
      var aa = D.stops.indexOf(c) >= 0 ? '*' : D.codons[c];
      (byAA[aa] = byAA[aa] || []).push(c);
    });
    var rows = [];
    Object.keys(byAA).sort().forEach(function (aa) {
      var cods = byAA[aa].sort();
      var total = cods.reduce(function (t, c) { return t + counts[c]; }, 0);
      cods.forEach(function (c) {
        rows.push({ codon: c, aa: aa, count: counts[c], rscu: total ? counts[c] / (total / cods.length) : null,
          fraction: total ? counts[c] / total : null });
      });
    });
    return { rows: rows, codons: Math.floor(s.length / 3), other: other };
  }

  // ---------- Protein (ProtParam) ----------

  function aaCounts(p) {
    var c = {};
    AA20.split('').forEach(function (a) { c[a] = 0; });
    for (var i = 0; i < p.length; i++) c[p.charAt(i)]++;
    return c;
  }

  function molecularWeight(p) {
    var w = 0;
    for (var i = 0; i < p.length; i++) w += D.proteinWeights[p.charAt(i)];
    return w - (p.length - 1) * D.water;
  }

  function chargeAtPH(p, pH) {
    var cnt = aaCounts(p);
    var pos = Object.assign({}, D.posPK), neg = Object.assign({}, D.negPK);
    if (p.charAt(0) in D.pKnterm) pos.Nterm = D.pKnterm[p.charAt(0)];
    if (p.charAt(p.length - 1) in D.pKcterm) neg.Cterm = D.pKcterm[p.charAt(p.length - 1)];
    var n = function (aa) { return aa === 'Nterm' || aa === 'Cterm' ? 1 : (cnt[aa] || 0); };
    var plus = 0, minus = 0;
    Object.keys(pos).forEach(function (aa) { plus += n(aa) / (Math.pow(10, pH - pos[aa]) + 1); });
    Object.keys(neg).forEach(function (aa) { minus += n(aa) / (Math.pow(10, neg[aa] - pH) + 1); });
    return plus - minus;
  }

  // Isoelectric point by bisection, as in Biopython / ExPASy (Bjellqvist pK values).
  function isoelectricPoint(p) {
    var pH = 7.775, lo = 4.05, hi = 12;
    while (hi - lo > 0.0001) {
      if (chargeAtPH(p, pH) > 0) lo = pH; else hi = pH;
      pH = (lo + hi) / 2;
    }
    return pH;
  }

  function proteinParams(p) {
    var cnt = aaCounts(p);
    var L = p.length;
    var pct = function (a) { return cnt[a] / L * 100; };
    var gravy = 0, inst = 0;
    for (var i = 0; i < L; i++) gravy += D.kd[p.charAt(i)];
    for (i = 0; i < L - 1; i++) inst += D.diwv[p.charAt(i)][p.charAt(i + 1)];
    var mw = molecularWeight(p);
    var extRed = cnt.W * 5500 + cnt.Y * 1490;
    var extCys = extRed + Math.floor(cnt.C / 2) * 125;
    return {
      length: L, mw: mw, pI: isoelectricPoint(p), charge7: chargeAtPH(p, 7),
      gravy: gravy / L, instability: 10 / L * inst,
      aromaticity: (cnt.F + cnt.W + cnt.Y) / L,
      aliphatic: pct('A') + 2.9 * pct('V') + 3.9 * (pct('I') + pct('L')),
      extReduced: extRed, extCystines: extCys,
      abs01Reduced: extRed / mw, abs01Cystines: extCys / mw,
      negative: cnt.D + cnt.E, positive: cnt.K + cnt.R,
      counts: cnt
    };
  }

  // Kyte-Doolittle hydropathy, window average plotted at the window centre (1-based).
  function hydropathy(p, win) {
    var out = [];
    for (var i = 0; i + win <= p.length; i++) {
      var s = 0;
      for (var j = i; j < i + win; j++) s += D.kd[p.charAt(j)];
      out.push({ x: i + (win + 1) / 2, y: s / win });
    }
    return out;
  }

  // ---------- Primers ----------

  // Salt correction method 5 of Biopython (SantaLucia 1998), with the Na-equivalent of
  // von Ahsen et al. (2001) for Mg2+ (corrected for dNTPs).
  function saltCorr5(seq, na, mg, dntps) {
    var mon = na;
    if (mg > 0 && dntps < mg) mon += 120 * Math.sqrt(mg - dntps);
    return 0.368 * (seq.length - 1) * Math.log(mon / 1000);
  }

  // Nearest-neighbour Tm of a primer with its perfect complement, as Biopython Tm_NN
  // (DNA_NN3 table, saltcorr=5). Concentrations: na, mg, dntps in mM; primerNM in nM.
  function tmNN(seq, opt) {
    var nn = D.nn;
    var dh = nn.init[0], ds = nn.init[1];
    if (gcContent(seq) === 0) { dh += nn['init_allA/T'][0]; ds += nn['init_allA/T'][1]; }
    else { dh += nn['init_oneG/C'][0]; ds += nn['init_oneG/C'][1]; }
    if (seq.charAt(0) === 'T') { dh += nn['init_5T/A'][0]; ds += nn['init_5T/A'][1]; }
    if (seq.charAt(seq.length - 1) === 'A') { dh += nn['init_5T/A'][0]; ds += nn['init_5T/A'][1]; }
    var ends = seq.charAt(0) + seq.charAt(seq.length - 1);
    var at = (ends.match(/[AT]/g) || []).length, gc = (ends.match(/[GC]/g) || []).length;
    dh += nn['init_A/T'][0] * at + nn['init_G/C'][0] * gc;
    ds += nn['init_A/T'][1] * at + nn['init_G/C'][1] * gc;
    var c = complement(seq);
    for (var i = 0; i < seq.length - 1; i++) {
      var key = seq.substr(i, 2) + '/' + c.substr(i, 2);
      var v = nn[key] || nn[key.split('').reverse().join('')];
      dh += v[0];
      ds += v[1];
    }
    ds += saltCorr5(seq, opt.na, opt.mg, opt.dntps);
    var k = opt.primerNM * 1e-9;
    return (1000 * dh) / (ds + 1.987 * Math.log(k)) - 273.15;
  }

  function tmWallace(seq) {
    var c = baseCounts(seq);
    return 2 * (c.A + c.T) + 4 * (c.G + c.C);
  }

  function longestRun(seq) {
    var best = 1, run = 1;
    for (var i = 1; i < seq.length; i++) {
      run = seq.charAt(i) === seq.charAt(i - 1) ? run + 1 : 1;
      best = Math.max(best, run);
    }
    return seq.length ? best : 0;
  }

  // Longest stretch of consecutive Watson-Crick pairs when a and b anneal antiparallel.
  // threeEnd: the longest stretch that includes the 3' terminal base of a (primer-dimer risk).
  function complementarity(a, b) {
    var br = b.split('').reverse().join('');
    var best = 0, three = 0;
    for (var off = -(br.length - 1); off < a.length; off++) {
      var run = 0;
      for (var i = Math.max(0, off); i < Math.min(a.length, off + br.length); i++) {
        var j = i - off;
        if (COMP[a.charAt(i)] === br.charAt(j) && a.charAt(i) !== 'N') {
          run++;
          best = Math.max(best, run);
          if (i === a.length - 1) three = Math.max(three, run);
        } else {
          run = 0;
        }
      }
    }
    return { any: best, threeEnd: three };
  }

  function primerInfo(seq, opt) {
    var last5 = seq.slice(-5);
    return {
      seq: seq, length: seq.length, gc: gcContent(seq),
      tm: seq.length >= 2 ? tmNN(seq, opt) : NaN, tmWallace: tmWallace(seq),
      gcClamp: (last5.match(/[GC]/g) || []).length,
      run: longestRun(seq), self: complementarity(seq, seq)
    };
  }

  // Exact-match PCR products of a forward and reverse primer on a template (both strands).
  function pcrProducts(template, fwd, rev) {
    var out = [];
    var find = function (s, sub) {
      var pos = [], i = s.indexOf(sub);
      while (i >= 0) { pos.push(i); i = s.indexOf(sub, i + 1); }
      return pos;
    };
    var rc = revcomp(rev);
    find(template, fwd).forEach(function (f) {
      find(template, rc).forEach(function (r) {
        if (r + rc.length > f + fwd.length - 1 && r >= f) {
          out.push({ start: f + 1, end: r + rc.length, length: r + rc.length - f, strand: '+',
            product: template.slice(f, r + rc.length) });
        }
      });
    });
    // Primers may be entered the other way round relative to the template.
    var rcf = revcomp(fwd);
    find(template, rev).forEach(function (r) {
      find(template, rcf).forEach(function (f) {
        if (f >= r && f + rcf.length > r + rev.length - 1) {
          out.push({ start: r + 1, end: f + rcf.length, length: f + rcf.length - r, strand: '−',
            product: template.slice(r, f + rcf.length) });
        }
      });
    });
    return out.sort(function (a, b) { return a.length - b.length; });
  }

  // ---------- Restriction sites ----------

  // Cut positions on the top strand of a linear sequence, as Biopython Restriction.search:
  // each value is the 1-based position of the first base after the cut.
  function cutSites(s, enz) {
    var site = enz.site, size = site.length, rc = revcomp(site);
    var cuts = [];
    for (var i = 0; i + size <= s.length; i++) {
      var w = s.substr(i, size);
      if (w === site) cuts.push(i + 1 + enz.fst5);
      else if (!enz.pal && w === rc) cuts.push(i + 1 - enz.fst3);
    }
    var L = s.length;
    cuts = cuts.filter(function (c) { var crick = c - enz.ovhg; return c > 1 && c <= L && crick > 1 && crick <= L; });
    return cuts.sort(function (a, b) { return a - b; });
  }

  function fragments(L, cuts) {
    var f = [], prev = 1;
    cuts.forEach(function (c) { f.push(c - prev); prev = c; });
    f.push(L - prev + 1);
    return f;
  }

  // ---------- Pairwise alignment ----------

  // Global (Needleman-Wunsch) or local (Smith-Waterman) alignment with affine gaps (Gotoh).
  // A gap of length k costs open + (k - 1) * extend (EMBOSS and Biopython convention).
  function align(a, b, opt) {
    var n = a.length, m = b.length;
    var local = opt.mode === 'local';
    var score = opt.type === 'protein'
      ? function (x, y) { return D.blosum62[x][y]; }
      : function (x, y) { return x === y && x !== 'N' ? opt.match : opt.mismatch; };
    var go = opt.open, ge = opt.extend;
    var NEG = -1e30;
    var W = m + 1;
    // Float32 keeps memory low; with integer or half-integer scores (as offered on the page) results are exact.
    var M = new Float32Array((n + 1) * W), X = new Float32Array((n + 1) * W), Y = new Float32Array((n + 1) * W);
    var tM = new Uint8Array((n + 1) * W), tX = new Uint8Array((n + 1) * W), tY = new Uint8Array((n + 1) * W);
    // States: 0 = M (match), 1 = X (gap in b), 2 = Y (gap in a), 3 = start (local)
    M[0] = 0; X[0] = NEG; Y[0] = NEG;
    for (var i = 1; i <= n; i++) {
      M[i * W] = NEG; Y[i * W] = NEG;
      X[i * W] = local ? NEG : -(go + (i - 1) * ge);
      tX[i * W] = i === 1 ? 0 : 1;
    }
    for (var j = 1; j <= m; j++) {
      M[j] = NEG; X[j] = NEG;
      Y[j] = local ? NEG : -(go + (j - 1) * ge);
      tY[j] = j === 1 ? 0 : 2;
    }
    var best = local ? 0 : NEG, bi = n, bj = m, bs = 0;
    for (i = 1; i <= n; i++) {
      var ai = a.charAt(i - 1);
      for (j = 1; j <= m; j++) {
        var k = i * W + j, d = k - W - 1;
        var s = score(ai, b.charAt(j - 1));
        var mm = M[d], mx = X[d], my = Y[d], t = 0, v = mm;
        if (mx > v) { v = mx; t = 1; }
        if (my > v) { v = my; t = 2; }
        if (local && v < 0) { v = 0; t = 3; }
        M[k] = v + s; tM[k] = t;
        var up = k - W;
        var o1 = M[up] - go, e1 = X[up] - ge, o2 = Y[up] - go;
        if (o1 >= e1 && o1 >= o2) { X[k] = o1; tX[k] = 0; } else if (e1 >= o2) { X[k] = e1; tX[k] = 1; } else { X[k] = o2; tX[k] = 2; }
        var lf = k - 1;
        var p1 = M[lf] - go, f1 = Y[lf] - ge, p2 = X[lf] - go;
        if (p1 >= f1 && p1 >= p2) { Y[k] = p1; tY[k] = 0; } else if (f1 >= p2) { Y[k] = f1; tY[k] = 2; } else { Y[k] = p2; tY[k] = 1; }
        if (local && M[k] > best) { best = M[k]; bi = i; bj = j; bs = 0; }
      }
    }
    if (!local) {
      var e = n * W + m;
      best = M[e]; bs = 0;
      if (X[e] > best) { best = X[e]; bs = 1; }
      if (Y[e] > best) { best = Y[e]; bs = 2; }
    }
    // Traceback
    var ra = [], rb = [];
    i = bi; j = bj;
    var st = bs;
    while (i > 0 || j > 0) {
      k = i * W + j;
      if (st === 0) {
        if (local && M[k] <= 0) break;
        var prev = tM[k];
        ra.push(a.charAt(i - 1)); rb.push(b.charAt(j - 1));
        i--; j--;
        if (prev === 3) break;
        st = prev;
      } else if (st === 1) {
        prev = tX[k];
        ra.push(a.charAt(i - 1)); rb.push('-');
        i--;
        st = prev;
      } else {
        prev = tY[k];
        ra.push('-'); rb.push(b.charAt(j - 1));
        j--;
        st = prev;
      }
      if (local && i === 0 && j === 0) break;
    }
    ra.reverse(); rb.reverse();
    var A = ra.join(''), B = rb.join('');
    var ident = 0, sim = 0, gaps = 0, mid = '';
    for (var q = 0; q < A.length; q++) {
      var x = A.charAt(q), y = B.charAt(q);
      if (x === '-' || y === '-') { gaps++; mid += ' '; }
      else if (x === y) { ident++; sim++; mid += '|'; }
      else if (opt.type === 'protein' && D.blosum62[x][y] > 0) { sim++; mid += ':'; }
      else mid += '.';
    }
    return { score: best, a: A, b: B, mid: mid, length: A.length, identity: ident, similarity: sim, gaps: gaps,
      startA: i + 1, startB: j + 1 };
  }

  // ---------- qPCR statistics ----------

  function mean(a) { return a.reduce(function (s, v) { return s + v; }, 0) / a.length; }
  function sd(a) {
    if (a.length < 2) return 0;
    var m = mean(a);
    return Math.sqrt(a.reduce(function (s, v) { return s + (v - m) * (v - m); }, 0) / (a.length - 1));
  }

  // Regularized incomplete beta function (continued fraction, Numerical Recipes).
  function betacf(a, b, x) {
    var qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < 1e-30) d = 1e-30;
    d = 1 / d;
    var h = d;
    for (var m = 1; m <= 300; m++) {
      var m2 = 2 * m, aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < 1e-30) d = 1e-30;
      c = 1 + aa / c; if (Math.abs(c) < 1e-30) c = 1e-30;
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < 1e-30) d = 1e-30;
      c = 1 + aa / c; if (Math.abs(c) < 1e-30) c = 1e-30;
      d = 1 / d;
      var del = d * c;
      h *= del;
      if (Math.abs(del - 1) < 3e-16) break;
    }
    return h;
  }
  function lgamma(x) {
    var c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    var y = x, t = x + 5.5;
    t -= (x + 0.5) * Math.log(t);
    var s = 1.000000000190015;
    for (var j = 0; j < 6; j++) s += c[j] / ++y;
    return -t + Math.log(2.5066282746310005 * s / x);
  }
  function ibeta(a, b, x) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    var bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
    return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a : 1 - bt * betacf(b, a, 1 - x) / b;
  }
  // Two-sided p-value of Welch's t-test.
  function welch(x, y) {
    if (x.length < 2 || y.length < 2) return null;
    var vx = Math.pow(sd(x), 2) / x.length, vy = Math.pow(sd(y), 2) / y.length;
    if (vx + vy === 0) return null;
    var t = (mean(x) - mean(y)) / Math.sqrt(vx + vy);
    var df = Math.pow(vx + vy, 2) / (vx * vx / (x.length - 1) + vy * vy / (y.length - 1));
    return { t: t, df: df, p: ibeta(df / 2, 0.5, df / (df + t * t)) };
  }

  // Relative expression (Livak 2^-ddCt, or Pfaffl with efficiencies). rows: {sample, group, target, refs: []}.
  // Reference genes are combined by the geometric mean of their relative quantities.
  function relativeExpression(rows, control, effTarget, effRef) {
    var ctl = rows.filter(function (r) { return r.group === control; });
    var ctlTarget = mean(ctl.map(function (r) { return r.target; }));
    var nRef = rows[0].refs.length;
    var ctlRef = [];
    for (var j = 0; j < nRef; j++) ctlRef.push(mean(ctl.map(function (r) { return r.refs[j]; })));
    var samples = rows.map(function (r) {
      var rqT = Math.pow(effTarget, ctlTarget - r.target);
      var logRef = 0;
      for (var j = 0; j < nRef; j++) logRef += Math.log(Math.pow(effRef, ctlRef[j] - r.refs[j]));
      var rqR = Math.exp(logRef / nRef);
      var refMean = mean(r.refs);
      return { sample: r.sample, group: r.group, dCt: r.target - refMean, fold: rqT / rqR, log2: Math.log(rqT / rqR) / Math.LN2 };
    });
    var groups = [];
    samples.forEach(function (s) { if (groups.indexOf(s.group) < 0) groups.push(s.group); });
    groups.sort(function (a, b) { return a === control ? -1 : b === control ? 1 : 0; });
    var summary = groups.map(function (g) {
      var ss = samples.filter(function (s) { return s.group === g; });
      var test = g === control ? null : welch(ss.map(function (s) { return s.log2; }),
        samples.filter(function (s) { return s.group === control; }).map(function (s) { return s.log2; }));
      var folds = ss.map(function (s) { return s.fold; });
      return { group: g, n: ss.length, meanFold: mean(folds), sdFold: sd(folds),
        meanLog2: mean(ss.map(function (s) { return s.log2; })), p: test ? test.p : null };
    });
    return { samples: samples, summary: summary };
  }

  // Standard curve: Ct against log10(quantity). Efficiency = 10^(-1/slope) - 1.
  function standardCurve(points) {
    var x = points.map(function (p) { return Math.log10(p.qty); });
    var y = points.map(function (p) { return p.ct; });
    var mx = mean(x), my = mean(y), sxx = 0, sxy = 0, syy = 0;
    for (var i = 0; i < x.length; i++) {
      sxx += (x[i] - mx) * (x[i] - mx);
      sxy += (x[i] - mx) * (y[i] - my);
      syy += (y[i] - my) * (y[i] - my);
    }
    var slope = sxy / sxx, intercept = my - slope * mx;
    return { slope: slope, intercept: intercept, r2: sxy * sxy / (sxx * syy),
      efficiency: Math.pow(10, -1 / slope) - 1, x: x, y: y };
  }

  var api = {
    readSeq: readSeq, complement: complement, revcomp: revcomp, gcContent: gcContent, baseCounts: baseCounts,
    gcWindows: gcWindows, translate: translate, findOrfs: findOrfs, codonUsage: codonUsage,
    molecularWeight: molecularWeight, chargeAtPH: chargeAtPH, isoelectricPoint: isoelectricPoint,
    proteinParams: proteinParams, hydropathy: hydropathy, tmNN: tmNN, tmWallace: tmWallace,
    complementarity: complementarity, primerInfo: primerInfo, pcrProducts: pcrProducts,
    cutSites: cutSites, fragments: fragments, align: align, welch: welch,
    relativeExpression: relativeExpression, standardCurve: standardCurve, mean: mean, sd: sd
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.SoybioSeq = api;
})(typeof window !== 'undefined' ? window : globalThis);

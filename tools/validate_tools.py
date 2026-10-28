#!/usr/bin/env python3
"""Check the Tools page calculations against Biopython and SciPy.

  python3 tools/validate_tools.py make  cases.json      # random test cases
  node    tools/validate_tools.js cases.json results.json
  python3 tools/validate_tools.py check cases.json results.json

Every number shown by the Tools page is compared with an independent reference
(Biopython 1.87 for sequence calculations, SciPy for statistics).
"""
import json
import math
import random
import sys
from collections import Counter

import numpy as np
from Bio import Restriction
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd
from scipy import stats

AA = "ACDEFGHIKLMNPQRSTVWY"
ENZ = ["AatII", "AflII", "ApaI", "AscI", "AvrII", "BamHI", "BglII", "BsaI", "BbsI", "BsmBI", "ClaI",
       "DraI", "EcoRI", "EcoRV", "HindIII", "HpaI", "KpnI", "MluI", "NcoI", "NdeI", "NheI", "NotI",
       "NsiI", "PacI", "PstI", "SacI", "SacII", "SalI", "SapI", "ScaI", "SmaI", "SpeI", "SphI",
       "StuI", "XbaI", "XhoI", "XmaI"]


def rdna(n, gc=0.5):
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(random.choices("ACGT", weights=w, k=n))


def make(path):
    random.seed(1)
    sites = [str(getattr(Restriction, e).site) for e in ENZ]
    cases = {
        "protein": ["".join(random.choices(AA, k=random.randint(2, 900))) for _ in range(200)]
                   + ["M", "MK", "DDDD", "RRRR", "CCCC", "WWYY", "GGGG"],
        "tm": [{"seq": rdna(random.randint(12, 40), random.uniform(0.2, 0.8)), "na": random.choice([10, 50, 100]),
                "mg": random.choice([0, 1.5, 3]), "dntps": random.choice([0, 0.2, 0.8]),
                "primerNM": random.choice([50, 200, 250, 500])} for _ in range(300)],
        "dna": [rdna(random.randint(1, 600)) for _ in range(100)] + ["ACGTNNACGT", "NNNN"],
        "orf": [rdna(random.randint(300, 3000), random.uniform(0.3, 0.6)) for _ in range(60)],
        "cds": [rdna(3 * random.randint(10, 400)) for _ in range(30)],
        "restriction": [],
        "align": [],
        "comp": [[rdna(random.randint(5, 30)), rdna(random.randint(5, 30))] for _ in range(200)],
        "welch": [], "curve": [], "rel": [], "pcr": [],
    }
    for _ in range(60):
        s = list(rdna(random.randint(20, 3000)))
        for _ in range(random.randint(0, 12)):  # plant sites, also near the ends
            site = random.choice(sites)
            pos = random.choice([0, max(0, len(s) - len(site)), random.randint(0, max(0, len(s) - len(site)))])
            if random.random() < 0.5:
                site = str(Seq(site).reverse_complement())
            s[pos:pos + len(site)] = site
        cases["restriction"].append("".join(s))
    for _ in range(80):
        if random.random() < 0.5:
            a = "".join(random.choices(AA, k=random.randint(1, 300)))
            b = list(a)
            for _ in range(random.randint(0, 60)):
                op = random.random()
                i = random.randint(0, len(b))
                if op < 0.4 and b: b[min(i, len(b) - 1)] = random.choice(AA)
                elif op < 0.7: b.insert(i, random.choice(AA))
                elif b: del b[min(i, len(b) - 1)]
            opt = {"type": "protein", "open": random.choice([10, 11, 5]), "extend": random.choice([0.5, 1, 2])}
        else:
            a = rdna(random.randint(1, 400))
            b = list(a)
            for _ in range(random.randint(0, 80)):
                op = random.random(); i = random.randint(0, len(b))
                if op < 0.4 and b: b[min(i, len(b) - 1)] = random.choice("ACGT")
                elif op < 0.7: b.insert(i, random.choice("ACGT"))
                elif b: del b[min(i, len(b) - 1)]
            opt = {"type": "dna", "match": 5, "mismatch": -4, "open": random.choice([10, 16]), "extend": random.choice([0.5, 4])}
        if not b: b = ["A"] if opt["type"] == "dna" else ["G"]
        opt["mode"] = random.choice(["global", "local"])
        cases["align"].append({"a": a, "b": "".join(b), "opt": opt})
    for _ in range(100):
        x = list(np.round(np.random.default_rng(len(cases["welch"])).normal(0, 1, random.randint(2, 8)), 4))
        y = list(np.round(np.random.default_rng(1000 + len(cases["welch"])).normal(random.uniform(-2, 2), random.uniform(0.3, 2), random.randint(2, 8)), 4))
        cases["welch"].append([[float(v) for v in x], [float(v) for v in y]])
    for _ in range(50):
        slope = random.uniform(-3.6, -3.1)
        pts = [{"qty": 10 ** k, "ct": round(35 + slope * k + random.gauss(0, 0.2), 3)} for k in range(random.randint(3, 7))]
        cases["curve"].append(pts)
    for _ in range(50):
        nref = random.randint(1, 3)
        rows = []
        for g in ["Control", "Salt", "Drought"]:
            for r in range(3):
                rows.append({"sample": f"{g}{r}", "group": g, "target": round(random.uniform(20, 32), 2),
                             "refs": [round(random.uniform(17, 25), 2) for _ in range(nref)]})
        cases["rel"].append({"rows": rows, "control": "Control", "et": random.choice([2, 1.95, 1.9]), "er": random.choice([2, 1.98])})
    for _ in range(60):
        t = rdna(random.randint(200, 1500))
        i = random.randint(0, len(t) - 120); j = random.randint(i + 60, len(t) - 20)
        f = t[i:i + 20]; r = str(Seq(t[j:j + 20]).reverse_complement())
        if random.random() < 0.3: f, r = r, f
        cases["pcr"].append({"t": t, "f": f, "r": r})
    json.dump(cases, open(path, "w"))
    print("cases written")


def close(a, b, tol):
    return (a is None and b is None) or (a is not None and b is not None and abs(a - b) <= tol)


def ref_orfs(s, min_aa):
    out = []
    L = len(s)
    for strand, t in ((1, s), (-1, str(Seq(s).reverse_complement()))):
        for f in range(3):
            start = None
            for i in range(f, L - 2, 3):
                cod = t[i:i + 3]
                if start is None and cod == "ATG":
                    start = i
                if start is not None and cod in ("TAA", "TAG", "TGA"):
                    if (i - start) // 3 >= min_aa:
                        a, b = start + 1, i + 3
                        st, en = (a, b) if strand == 1 else (L - b + 1, L - a + 1)
                        out.append([strand, st, en, str(Seq(t[start:i]).translate())])
                    start = None
    return sorted(map(tuple, out))


def check(cpath, rpath):
    cases, res = json.load(open(cpath)), json.load(open(rpath))
    fails = Counter()
    n = Counter()

    for p, r in zip(cases["protein"], res["protein"]):
        pa = ProteinAnalysis(p)
        ref = {"mw": pa.molecular_weight(), "pI": pa.isoelectric_point(), "charge7": pa.charge_at_pH(7.0),
               "gravy": pa.gravy(), "aromaticity": pa.aromaticity(),
               "instability": pa.instability_index() if len(p) > 1 else None}
        for k, v in ref.items():
            n[k] += 1
            if v is not None and not close(r[k], v, 1e-6): fails[k] += 1
        n["ext"] += 1
        if list(pa.molar_extinction_coefficient()) != r["ext"]: fails["ext"] += 1
        hp = [sum(kd[c] for c in p[i:i + 9]) / 9 for i in range(len(p) - 8)]
        n["hydropathy"] += 1
        if len(hp) != len(r["hydro9"]) or any(abs(a - b) > 1e-9 for a, b in zip(hp, r["hydro9"])): fails["hydropathy"] += 1

    for c, t in zip(cases["tm"], res["tm"]):
        ref = mt.Tm_NN(c["seq"], Na=c["na"], Mg=c["mg"], dNTPs=c["dntps"], dnac1=c["primerNM"], dnac2=0, saltcorr=5)
        n["Tm_NN"] += 1
        if not close(t, ref, 1e-6): fails["Tm_NN"] += 1

    for s, tr, rc, gc in zip(cases["dna"], res["translate"], res["revcomp"], res["gc"]):
        n["translate"] += 1; n["revcomp"] += 1; n["gc"] += 1
        cut = s[:len(s) // 3 * 3]
        if tr != (str(Seq(cut).translate()) if cut else ""): fails["translate"] += 1
        if rc != str(Seq(s).reverse_complement()): fails["revcomp"] += 1
        if not close(gc, gc_fraction(s, "remove") * 100, 1e-9): fails["gc"] += 1

    for s, r in zip(cases["orf"], res["orfs"]):
        n["orfs"] += 1
        if sorted(map(tuple, r)) != ref_orfs(s, 30): fails["orfs"] += 1

    for s, r in zip(cases["cds"], res["codons"]):
        cnt = Counter(s[i:i + 3] for i in range(0, len(s) - 2, 3))
        n["codon usage"] += 1
        if any(cnt.get(c, 0) != k for c, k, _ in r): fails["codon usage"] += 1

    for s, r in zip(cases["restriction"], res["cuts"]):
        seq = Seq(s)
        for e in ENZ:
            n["restriction"] += 1
            if getattr(Restriction, e).search(seq, linear=True) != r[e]: fails["restriction"] += 1

    blosum = substitution_matrices.load("BLOSUM62")
    for c, r in zip(cases["align"], res["align"]):
        o = c["opt"]
        al = PairwiseAligner(mode=o["mode"], open_gap_score=-o["open"], extend_gap_score=-o["extend"])
        if o["type"] == "protein": al.substitution_matrix = blosum
        else: al.match_score, al.mismatch_score = o["match"], o["mismatch"]
        ref = al.score(c["a"], c["b"])
        n["alignment score"] += 1; n["alignment traceback"] += 1
        if not close(r["score"], ref, 1e-6): fails["alignment score"] += 1
        ok = close(r["rescore"], r["score"], 1e-6)
        if o["mode"] == "global": ok = ok and r["degapA"] == c["a"] and r["degapB"] == c["b"]
        else: ok = ok and c["a"][r["startA"] - 1:].startswith(r["degapA"]) and c["b"][r["startB"] - 1:].startswith(r["degapB"])
        if not ok: fails["alignment traceback"] += 1

    for (a, b), r in zip(cases["comp"], res["comp"]):
        comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
        br = b[::-1]; best = three = 0
        for off in range(-(len(br) - 1), len(a)):
            run = 0
            for i in range(max(0, off), min(len(a), off + len(br))):
                if comp[a[i]] == br[i - off]:
                    run += 1; best = max(best, run)
                    if i == len(a) - 1: three = max(three, run)
                else:
                    run = 0
        n["complementarity"] += 1
        if [best, three] != [r["any"], r["threeEnd"]]: fails["complementarity"] += 1

    for (x, y), p in zip(cases["welch"], res["welch"]):
        n["Welch t-test"] += 1
        if not close(p, float(stats.ttest_ind(x, y, equal_var=False).pvalue), 1e-8): fails["Welch t-test"] += 1

    for pts, r in zip(cases["curve"], res["curve"]):
        lr = stats.linregress([math.log10(p["qty"]) for p in pts], [p["ct"] for p in pts])
        ref = [lr.slope, lr.intercept, lr.rvalue ** 2, 10 ** (-1 / lr.slope) - 1]
        n["standard curve"] += 1
        if any(not close(a, b, 1e-9) for a, b in zip(r, ref)): fails["standard curve"] += 1

    for c, folds in zip(cases["rel"], res["rel"]):
        ctl = [r for r in c["rows"] if r["group"] == c["control"]]
        mt_ = np.mean([r["target"] for r in ctl]); mr = np.mean([r["refs"] for r in ctl], axis=0)
        ref = []
        for r in c["rows"]:
            rqt = c["et"] ** (mt_ - r["target"])
            rqr = np.exp(np.mean([np.log(c["er"] ** (mr[j] - r["refs"][j])) for j in range(len(mr))]))
            ref.append(rqt / rqr)
        n["relative expression"] += 1
        if any(not close(a, b, 1e-9) for a, b in zip(folds, ref)): fails["relative expression"] += 1

    for c, r in zip(cases["pcr"], res["pcr"]):
        t, f, rv = c["t"], c["f"], c["r"]
        exp = []
        rc = str(Seq(rv).reverse_complement()); rcf = str(Seq(f).reverse_complement())
        fs = [i for i in range(len(t)) if t.startswith(f, i)]; rs = [i for i in range(len(t)) if t.startswith(rc, i)]
        exp += [[a + 1, b + len(rc)] for a in fs for b in rs if b >= a and b + len(rc) > a + len(f) - 1]
        rs2 = [i for i in range(len(t)) if t.startswith(rv, i)]; fs2 = [i for i in range(len(t)) if t.startswith(rcf, i)]
        exp += [[a + 1, b + len(rcf)] for a in rs2 for b in fs2 if b >= a and b + len(rcf) > a + len(rv) - 1]
        n["in-silico PCR"] += 1
        if sorted(map(tuple, exp)) != sorted(map(tuple, r)) or not r: fails["in-silico PCR"] += 1

    print(f"{'check':24s} {'tests':>7s} {'failed':>7s}")
    for k in n:
        print(f"{k:24s} {n[k]:7d} {fails[k]:7d}")
    print("ALL PASSED" if not sum(fails.values()) else "FAILURES FOUND")


if __name__ == "__main__":
    make(sys.argv[2]) if sys.argv[1] == "make" else check(sys.argv[2], sys.argv[3])

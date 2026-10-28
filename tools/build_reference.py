#!/usr/bin/env python3
"""Build the data for the qPCR reference gene finder.

Uses only files already on the website:
  data/abiotic-rnaseq/abiotic_stress_TPM_matrix.tsv.gz, meta.json, abiotic_stress_DEGs.tsv.gz
  data/zh13-wm82a4/Wm82a4_to_ZH13.tsv (descriptions and ZH13 IDs)
Writes data/abiotic-rnaseq/reference_genes.json.

For each subset of samples (stress x tissue), genes are ranked by the standard deviation
of log2(TPM + 1) across samples (lower = more stable). Only expressed genes are ranked
(mean TPM >= 10 and TPM >= 2 in every sample of the subset).

Run from the website folder:  python3 tools/build_reference.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "data" / "abiotic-rnaseq"
TOP_N = 100
MIN_MEAN, MIN_TPM = 10, 2
LEVELS = [10, 30, 100]  # ranked lists for these minimum mean TPM values
SUBSETS = [("All", "All"), ("Salt", "All"), ("Drought", "All"), ("Heat", "All"),
           ("All", "Leaf"), ("All", "Root"), ("Salt", "Leaf"), ("Salt", "Root"),
           ("Drought", "Leaf"), ("Drought", "Root"), ("Heat", "Leaf")]


def main():
    tpm = pd.read_csv(D / "abiotic_stress_TPM_matrix.tsv.gz", sep="\t", index_col=0)
    meta = json.loads((D / "meta.json").read_text(encoding="utf-8"))
    degs = pd.read_csv(D / "abiotic_stress_DEGs.tsv.gz", sep="\t")
    ids = pd.read_csv(ROOT / "data" / "zh13-wm82a4" / "Wm82a4_to_ZH13.tsv", sep="\t").set_index("Wm82.a4_gene_ID")
    samples = pd.DataFrame(meta["samples"])
    assert list(samples["id"]) == list(tpm.columns)

    genes = tpm.index.to_numpy()
    out = {"fields": ["gene", "mean_tpm", "cv", "sd_log2", "max_fold", "de"],
           "min_mean": MIN_MEAN, "min_tpm": MIN_TPM, "subsets": {}, "genes": {}}
    listed = set()
    for stress, tissue in SUBSETS:
        keep = np.ones(len(samples), bool)
        if stress != "All":
            keep &= (samples["stress"] == stress).to_numpy()
        if tissue != "All":
            keep &= (samples["tissue"] == tissue).to_numpy()
        X = tpm.to_numpy()[:, keep]
        mean = X.mean(axis=1)
        lo, hi = X.min(axis=1), X.max(axis=1)
        sd_log = np.log2(X + 1).std(axis=1, ddof=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            cv = X.std(axis=1, ddof=1) / mean * 100
        fold = (hi + 1) / (lo + 1)
        expressed = (mean >= MIN_MEAN) & (lo >= MIN_TPM)

        cons = [c["id"] for c in meta["contrasts"]
                if (stress == "All" or c["stress"] == stress) and (tissue == "All" or c["tissue"] == tissue)]
        de = degs[degs["contrast"].isin(cons)].groupby("gene")["contrast"].nunique()
        de = de.reindex(genes).fillna(0).astype(int).to_numpy()

        top = {}
        for level in LEVELS:
            idx = np.where(expressed & (mean >= level))[0]
            idx = idx[np.argsort(sd_log[idx], kind="stable")][:TOP_N]
            top[str(level)] = [[genes[i], round(float(mean[i]), 1), round(float(cv[i]), 1), round(float(sd_log[i]), 3),
                                round(float(fold[i]), 2), int(de[i])] for i in idx]
            listed.update(genes[idx])
        out["subsets"][f"{stress}|{tissue}"] = {
            "samples": int(keep.sum()),
            "contrasts": len(cons),
            "expressed": int(expressed.sum()),
            # Percentiles 0..100 of sd_log2 among expressed genes, to rank any gene.
            "q": [round(float(v), 4) for v in np.quantile(sd_log[expressed], np.linspace(0, 1, 101))],
            "top": top,
        }
        print(f"{stress:8s} {tissue:5s} samples={keep.sum():3d} expressed={expressed.sum():5d} "
              f"best={top['10'][0][0]} sd={top['10'][0][3]} n100={len(top['100'])}")

    for g in sorted(listed):
        desc = ids["Wm82.a4_description"].get(g, "")
        zh = ids["ZH13_gene_ID"].get(g, "")
        out["genes"][g] = ["" if pd.isna(zh) else zh, "" if pd.isna(desc) else str(desc)]
    path = D / "reference_genes.json"
    path.write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print(f"{path.name}: {path.stat().st_size / 1e3:.0f} KB, {len(listed)} genes listed")


if __name__ == "__main__":
    main()

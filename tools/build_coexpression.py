#!/usr/bin/env python3
"""Find the most co-expressed genes of each gene in the abiotic stress RNA-seq data.

Method:
  1. Keep genes with mean TPM >= 2 across the 192 samples.
  2. Use log2(TPM + 1) and subtract each gene's mean within every study and tissue
     (BioProject + tissue). This removes differences between studies and tissues, so the
     correlation reflects how genes change together in response to stress treatments.
  3. Drop genes that hardly change after this step (SD < 0.15).
  4. Pearson correlation between all remaining genes; keep the 30 most correlated partners.

Output: data/abiotic-rnaseq/coexp/<key>.json with {gene: [[partner, r], ...]}, using the
same file keys as chunks/ (see build_expression.py), and coexp/info.json.

Run from the website folder (after build_expression.py):  python3 tools/build_coexpression.py
"""
import json
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "data" / "abiotic-rnaseq"
OUT = D / "coexp"
MIN_MEAN = 2
MIN_SD = 0.15
TOP = 30
BLOCK = 1000


def chunk_key(gene, bucket):
    m = re.match(r"^Glyma\.(\d{2}|U)G?(\d+)$", gene)
    return f"{m.group(1)}_{int(m.group(2)) // bucket:02d}"


def main():
    meta = json.loads((D / "meta.json").read_text(encoding="utf-8"))
    tpm = pd.read_csv(D / "abiotic_stress_TPM_matrix.tsv.gz", sep="\t", index_col=0)
    tpm = tpm[tpm.mean(axis=1) >= MIN_MEAN]
    x = np.log2(tpm.to_numpy(dtype=np.float64) + 1)

    study = pd.Series([s["project"] + "|" + s["tissue"] for s in meta["samples"]])
    for key in study.unique():
        cols = np.where(study == key)[0]
        x[:, cols] -= x[:, cols].mean(axis=1, keepdims=True)
    sd = x.std(axis=1, ddof=1)
    keep = sd >= MIN_SD
    genes = tpm.index.to_numpy()[keep]
    z = x[keep]
    z = (z - z.mean(axis=1, keepdims=True)) / z.std(axis=1, keepdims=True)
    z = (z / np.sqrt(z.shape[1])).astype(np.float32)
    n = len(genes)
    print(f"genes with mean TPM >= {MIN_MEAN}: {len(tpm)}, used after SD filter: {n}")

    partners = {}
    for start in range(0, n, BLOCK):
        r = z[start:start + BLOCK] @ z.T
        for i in range(r.shape[0]):
            r[i, start + i] = -2  # exclude the gene itself
        top = np.argpartition(-r, TOP, axis=1)[:, :TOP]
        for i in range(r.shape[0]):
            idx = top[i][np.argsort(-r[i, top[i]])]
            partners[genes[start + i]] = [[genes[j], round(float(r[i, j]), 3)] for j in idx]

    shutil.rmtree(OUT, ignore_errors=True)
    OUT.mkdir(parents=True)
    files = {}
    for g, p in partners.items():
        files.setdefault(chunk_key(g, meta["bucket"]), {})[g] = p
    for key, obj in files.items():
        (OUT / f"{key}.json").write_text(json.dumps(obj, separators=(",", ":")), encoding="utf-8")
    best = np.array([p[0][1] for p in partners.values()])
    info = {"genes": n, "top": TOP, "min_mean_tpm": MIN_MEAN, "min_sd": MIN_SD}
    (OUT / "info.json").write_text(json.dumps(info), encoding="utf-8")
    total = sum(f.stat().st_size for f in OUT.iterdir())
    print(f"files: {len(files)}, total {total / 1e6:.1f} MB; best-partner r median {np.median(best):.2f}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build the website files for the abiotic stress RNA-seq expression data.

Input: the extracted folder "newly_processed_abiotic_RNAseq" (from the zip file) with
  matrices/gene_TPM_matrix.tsv, matrices/gene_counts_matrix.tsv, matrices/alignment_stats.tsv,
  SRR_bioproject_list.xlsx (OpenDocument spreadsheet),
  downstream/tables/contrasts.csv, DE_*.tsv, all_DEGs_long.tsv, deg_stress_membership.tsv
Also uses data/zh13-wm82a4/Wm82a4_to_ZH13.tsv for gene descriptions and ZH13 IDs.

Output (data/abiotic-rnaseq/):
  meta.json          samples, condition groups, contrasts and example genes
  idmap.json         ZH13 and Wm82.a2.v1 IDs -> Wm82.a4 IDs
  chunks/*.json      TPM values and DE results per gene, about 100 genes per file
  download files     TPM and count matrices, sample information, DEG table

Run from the website folder:
  python3 tools/build_expression.py <path to newly_processed_abiotic_RNAseq>
"""
import gzip
import json
import math
import re
import shutil
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "abiotic-rnaseq"
ID_TABLE = ROOT / "data" / "zh13-wm82a4" / "Wm82a4_to_ZH13.tsv"

GENE_PREFIX = "glyma.Wm82.gnm4.ann1."
STRESS_ORDER = ["Salt", "Drought", "Heat"]
BUCKET = 10000  # gene number range per chunk file


# ---------- Reading inputs ----------

def read_sra_sheet(path):
    """Read the SRA sample list. The file is an OpenDocument spreadsheet named .xlsx."""
    t = "urn:oasis:names:tc:opendocument:xmlns:table:1.0"
    x = "urn:oasis:names:tc:opendocument:xmlns:text:1.0"
    root = ET.fromstring(zipfile.ZipFile(path).read("content.xml"))
    table = next(root.iter(f"{{{t}}}table"))
    rows = []
    for row in table.iter(f"{{{t}}}table-row"):
        cells = []
        for c in row:
            if c.tag not in (f"{{{t}}}table-cell", f"{{{t}}}covered-table-cell"):
                continue
            rep = int(c.get(f"{{{t}}}number-columns-repeated", "1"))
            txt = " ".join("".join(p.itertext()) for p in c.findall(f"{{{x}}}p")).strip()
            cells.extend([txt] * min(rep, 20))
        rows.append(cells[:8])
    head = next(i for i, r in enumerate(rows) if r and r[0] == "Stress")
    data = [r + [""] * (8 - len(r)) for r in rows[head + 1:] if any(r)]
    df = pd.DataFrame(data, columns=rows[head]).replace("", pd.NA)
    df[["Stress", "Bioproject"]] = df[["Stress", "Bioproject"]].ffill()
    return df.dropna(subset=["SRA Id"]).reset_index(drop=True)


# ---------- Text clean-up ----------

def tidy_duration(s):
    m = re.match(r"^\s*(.*?)\s*(Hours?|Days?)\s*$", str(s))
    if not m:
        return str(s)
    unit = "h" if m.group(2).startswith("Hour") else "d"
    return f"{m.group(1)} {unit}"


def duration_hours(s):
    m = re.match(r"^\s*([\d.]+)\D*?(Hours?|Days?|h|d)\b", str(s))
    if not m:
        return 0.0
    v = float(m.group(1))
    return v * 24 if m.group(2).lower().startswith("d") else v


def tidy_intensity(s):
    s = str(s).strip()
    rules = [
        (r"^(\d+(?:\.\d+)?) mM/L$", r"\1 mM"),
        (r"^Soil with (\d+(?:\.\d+)?)% relative water$", r"soil water \1%"),
        (r"^(\d+) Degree Celcious$", r"\1 °C"),
        (r"^Withholding water supply$", "water withheld"),
        (r"^Normal water supply$", "normal watering"),
    ]
    for pat, rep in rules:
        s = re.sub(pat, rep, s)
    return s


def tidy_contrast_condition(raw):
    """'7.Days.Soil.with.24.3.relative.water' -> '7 d, soil water 24.3%'."""
    s = raw.replace(".", " ")
    m = re.match(r"^(.*?)\s+(Hours|Days)\s+(.*)$", s)
    when, unit, what = m.group(1), m.group(2), m.group(3)
    nums = when.replace(" and ", " & ").split()
    nums = [n for n in nums if n != "&"]
    when = nums[0] if len(nums) == 1 else ", ".join(nums[:-1]) + " and " + nums[-1]
    unit = "h" if unit == "Hours" else "d"
    rules = [
        (r"^Soil with (\d+) (\d) relative water$", r"soil water \1.\2%"),
        (r"^(\d+) mM L$", r"\1 mM"),
        (r"^(\d+) Degree Celcious$", r"\1 °C"),
        (r"^Withholding water supply$", "water withheld"),
        (r"^(\d+) PEG 6000$", r"\1% PEG-6000"),
    ]
    for pat, rep in rules:
        what = re.sub(pat, rep, what)
    return f"{when} {unit}, {what}"


def tidy_tissue(s):
    return "Leaf, root and stem" if s == "Leaf, Root and Stem" else s


# ---------- Number formatting ----------

def r_tpm(v):
    if v == 0:
        return 0
    if v < 100:
        return round(v, 2)
    if v < 10000:
        return round(v, 1)
    return round(v)


def r_lfc(v):
    return None if pd.isna(v) else round(float(v), 2)


def r_padj(v):
    if pd.isna(v):
        return None
    v = float(v)
    return 0 if v == 0 else float(f"{v:.2g}")


def chunk_key(gene):
    m = re.match(r"^Glyma\.(\d{2}|U)G?(\d+)$", gene)
    if not m:
        raise ValueError(f"Unexpected gene ID: {gene}")
    return f"{m.group(1)}_{int(m.group(2)) // BUCKET:02d}"


# ---------- Main ----------

def main(src):
    src = Path(src)
    mat = src / "matrices"
    tab = src / "downstream" / "tables"

    tpm = pd.read_csv(mat / "gene_TPM_matrix.tsv", sep="\t", index_col=0)
    counts = pd.read_csv(mat / "gene_counts_matrix.tsv", sep="\t", index_col=0)
    stats = pd.read_csv(mat / "alignment_stats.tsv", sep="\t").set_index("sample")
    for df in (tpm, counts):
        df.index = df.index.str.replace(GENE_PREFIX, "", regex=False)
        df.index.name = "gene"
    counts = counts.reindex(index=tpm.index, columns=tpm.columns)
    sra = read_sra_sheet(src / "SRR_bioproject_list.xlsx")
    assert len(sra) == tpm.shape[1], "SRA sheet and TPM matrix have different sample numbers"

    # Samples: TPM columns are in the same order as the SRA sheet rows (checked below).
    stress_name = {"Salt Stress": "Salt", "Drought Stress": "Drought", "Heat Stress": "Heat"}
    trt_name = {"Control": "Control", "Salt stress": "Salt", "Drought stress": "Drought", "Heat stress": "Heat"}
    samples = []
    for i, col in enumerate(tpm.columns):
        r = sra.iloc[i]
        code = col.split("_")
        assert r["Bioproject"].endswith(code[3]), (col, r["Bioproject"])
        assert (code[1] == "CK") == (r["Treatment"] == "Control"), (col, r["Treatment"])
        samples.append({
            "id": col,
            "srr": r["SRA Id"],
            "stress": stress_name[r["Stress"]],
            "project": r["Bioproject"],
            "treatment": trt_name[r["Treatment"]],
            "tissue": tidy_tissue(r["Tissue type"]),
            "stage": r["Developmental stages"],
            "duration": tidy_duration(r["Treatment duration"]),
            "intensity": tidy_intensity(r["Stress intensity"]),
            "reads": int(stats.loc[col, "total_reads"]),
            "aligned": float(stats.loc[col, "align_rate"]),
        })

    # Condition groups (biological replicates of one condition in one study).
    proj_order = list(dict.fromkeys(s["project"] for s in samples))
    tissue_order = list(dict.fromkeys(s["tissue"] for s in samples))
    keyed = {}
    for i, s in enumerate(samples):
        k = (s["project"], s["tissue"], s["treatment"], s["duration"],
             "" if s["treatment"] == "Control" else s["intensity"])
        keyed.setdefault(k, []).append(i)
    groups = []
    for k, idx in keyed.items():
        s = samples[idx[0]]
        label = f"Control, {s['duration']}" if s["treatment"] == "Control" else f"{s['intensity']}, {s['duration']}"
        groups.append({
            "stress": s["stress"], "project": s["project"], "tissue": s["tissue"], "stage": s["stage"],
            "treatment": s["treatment"], "duration": s["duration"], "label": label, "samples": idx,
            "_sort": (STRESS_ORDER.index(s["stress"]), proj_order.index(s["project"]),
                      tissue_order.index(s["tissue"]), s["treatment"] != "Control",
                      duration_hours(s["duration"]), label),
        })
    groups.sort(key=lambda g: g["_sort"])
    for gi, g in enumerate(groups):
        del g["_sort"]
        for i in g["samples"]:
            samples[i]["group"] = gi

    # Contrasts and DE results.
    con = pd.read_csv(tab / "contrasts.csv")
    con["tissue"] = con["tissue"].map(tidy_tissue)
    con["cond"] = con["condition"].map(tidy_contrast_condition)
    con["_sort"] = [(STRESS_ORDER.index(r.stress), proj_order.index(r.project),
                     tissue_order.index(r.tissue) if r.tissue in tissue_order else 99,
                     duration_hours(r.cond)) for r in con.itertuples()]
    con = con.sort_values("_sort").reset_index(drop=True)
    contrasts = [{"id": r.id, "stress": r.stress, "project": r.project, "tissue": r.tissue,
                  "condition": r.cond, "n_ctrl": int(r.n_ctrl), "n_trt": int(r.n_trt)}
                 for r in con.itertuples()]
    lfc = pd.DataFrame(index=tpm.index, dtype=float)
    padj = pd.DataFrame(index=tpm.index, dtype=float)
    for c in contrasts:
        de = pd.read_csv(tab / f"DE_{c['id']}.tsv", sep="\t").set_index("gene")
        lfc[c["id"]] = de["log2FoldChange"].reindex(tpm.index)
        padj[c["id"]] = de["padj"].reindex(tpm.index)
    sig = (padj < 0.05) & (lfc.abs() >= 1)

    # Gene descriptions and ZH13 IDs from the gene ID table.
    ids = pd.read_csv(ID_TABLE, sep="\t").set_index("Wm82.a4_gene_ID")
    desc = ids["Wm82.a4_description"].reindex(tpm.index).fillna("")
    zh = ids["ZH13_gene_ID"].reindex(tpm.index).fillna("")

    # Example genes: core multi-stress DEGs that respond in the most contrasts.
    core = [g for g in (tab / "core_multistress_genes.txt").read_text().split() if g in tpm.index]
    score = pd.DataFrame({"n": sig.loc[core].sum(axis=1),
                          "lfc": lfc.loc[core].abs().mean(axis=1),
                          "tpm": tpm.loc[core].mean(axis=1)})
    score = score[score["tpm"] >= 5].sort_values(["n", "lfc"], ascending=False)
    examples = {"single": score.index[0], "multiple": list(score.index[:30])}

    # Write output.
    if OUT.exists():
        shutil.rmtree(OUT / "chunks", ignore_errors=True)
    (OUT / "chunks").mkdir(parents=True, exist_ok=True)

    chunks = {}
    for gene in tpm.index:
        chunks.setdefault(chunk_key(gene), []).append(gene)
    tpm_rows = tpm.to_numpy()
    pos = {g: i for i, g in enumerate(tpm.index)}
    for key, genes in chunks.items():
        out = {}
        for g in genes:
            i = pos[g]
            out[g] = {
                "z": zh.iat[i],
                "d": desc.iat[i],
                "t": [r_tpm(float(v)) for v in tpm_rows[i]],
                "f": [r_lfc(v) for v in lfc.iloc[i]],
                "p": [r_padj(v) for v in padj.iloc[i]],
            }
        with open(OUT / "chunks" / f"{key}.json", "w", encoding="utf-8") as fh:
            json.dump(out, fh, ensure_ascii=False, separators=(",", ":"))

    meta = {
        "version": "1.0",
        "genes": int(len(tpm)),
        "samples": samples,
        "groups": groups,
        "contrasts": contrasts,
        "examples": examples,
        "bucket": BUCKET,
    }
    with open(OUT / "meta.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, separators=(",", ":"))

    idmap = {
        "zh": {z: g for g, z in ids["ZH13_gene_ID"].dropna().items()},
        "a2": {a: g for g, a in ids["Wm82.a2.v1_ID"].dropna().items() if a != g},
    }
    with open(OUT / "idmap.json", "w", encoding="utf-8") as fh:
        json.dump(idmap, fh, separators=(",", ":"))

    # Download files.
    with gzip.open(OUT / "abiotic_stress_TPM_matrix.tsv.gz", "wt", encoding="utf-8") as fh:
        tpm.to_csv(fh, sep="\t", float_format="%.4f")
    with gzip.open(OUT / "abiotic_stress_counts_matrix.tsv.gz", "wt", encoding="utf-8") as fh:
        counts.to_csv(fh, sep="\t")
    info = pd.DataFrame(samples).rename(columns={
        "id": "sample", "srr": "SRA_run", "project": "BioProject", "reads": "reads_after_trimming",
        "aligned": "alignment_rate_%"})
    info["condition_group"] = [groups[s["group"]]["label"] for s in samples]
    info = info[["sample", "SRA_run", "BioProject", "stress", "treatment", "tissue", "stage",
                 "duration", "intensity", "condition_group", "reads_after_trimming", "alignment_rate_%"]]
    info.to_csv(OUT / "abiotic_stress_sample_information.tsv", sep="\t", index=False)
    degs = pd.read_csv(tab / "all_DEGs_long.tsv", sep="\t")
    degs = degs[["gene", "id", "stress", "tissue", "project", "log2FoldChange", "padj", "dir"]]
    degs = degs.drop_duplicates().rename(columns={"id": "contrast", "dir": "direction"})
    with gzip.open(OUT / "abiotic_stress_DEGs.tsv.gz", "wt", encoding="utf-8") as fh:
        degs.to_csv(fh, sep="\t", index=False, float_format="%.6g")
    shutil.copy(tab / "deg_stress_membership.tsv", OUT / "abiotic_stress_DEG_membership.tsv")

    # Summary for the web pages.
    print("genes:", len(tpm), " samples:", len(samples), " groups:", len(groups),
          " contrasts:", len(contrasts), " chunk files:", len(chunks))
    print("studies:", len(proj_order), " by stress:",
          pd.Series([s["stress"] for s in samples]).value_counts().to_dict())
    print("DEGs per stress (any contrast):",
          {st: int(sig[[c["id"] for c in contrasts if c["stress"] == st]].any(axis=1).sum())
           for st in STRESS_ORDER})
    print("DEGs in any contrast:", int(sig.any(axis=1).sum()))
    print("example single:", examples["single"], desc[examples["single"]][:80])
    for f in sorted(OUT.iterdir()):
        if f.is_file():
            print(f"{f.stat().st_size / 1e6:8.2f} MB  {f.name}")
    sizes = [p.stat().st_size for p in (OUT / "chunks").iterdir()]
    print(f"chunks: total {sum(sizes) / 1e6:.1f} MB, max {max(sizes) / 1e3:.0f} KB")


if __name__ == "__main__":
    main(sys.argv[1])

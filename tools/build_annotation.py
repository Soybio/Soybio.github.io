#!/usr/bin/env python3
"""Build compact GO, KEGG, Pfam and stress-response files for the Gene Analysis page.

Inputs (in the website folder; the raw files are not published):
  Williams_82.GO.txt.gz              Gene, GO, Level, Description (parent terms included)
  Williams_82.KEGG.txt.gz            Gene, KO, Pathway, Description
  Williams82.a2.v1.pfam.domain.txt   GeneId (Wm82.a2.v1), PFAM
  data/zh13-wm82a4/Wm82a4_to_ZH13.tsv           to convert Wm82.a2.v1 IDs to Wm82.a4
  data/abiotic-rnaseq/abiotic_stress_DEGs.tsv.gz, meta.json   for stress responses

Output (data/annotation/):
  go_bp.json, go_mf.json, go_cc.json   {"terms": [[id, name], ...], "genes": {gene: [term index, ...]}}
  kegg.json    {"terms": [[map id, name, global?], ...], "genes": {gene: [index, ...]}, "ko": {gene: "K..,K.."}}
  pfam.json    {"terms": [[domain, ""], ...], "genes": {gene: [index, ...]}}
  genes.json   all Wm82.a4 genes in the expression data (the background for stress tests)
  Wm82a4_GO.tsv.gz, Wm82a4_KEGG.tsv.gz, Wm82a4_Pfam.tsv.gz   download files
  stress.json  {"genes": {gene: [salt up, salt down, drought up, drought down, heat up, heat down]},
                "comparisons": {"Salt": n, ...}, "universe": number of genes in the expression data}

Run from the website folder:  python3 tools/build_annotation.py
"""
import gzip
import json
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "annotation"
NAMESPACES = {"Biological Process": "go_bp", "Molecular Function": "go_mf", "Cellular Component": "go_cc"}
# KEGG global and overview maps (01100 Metabolic pathways, ...) are shown but not tested.
GLOBAL_MAPS = ("map011", "map012")


def write(name, obj):
    path = OUT / name
    path.write_text(json.dumps(obj, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    raw = path.stat().st_size
    print(f"{name:12s} {raw / 1e6:6.2f} MB (gzip {len(gzip.compress(path.read_bytes())) / 1e6:.2f} MB)")


def encode(pairs, labels):
    """pairs: DataFrame with columns gene, term. labels: {term: label}. Terms sorted by ID."""
    terms = sorted(pairs["term"].unique())
    index = {t: i for i, t in enumerate(terms)}
    genes = {g: sorted(index[t] for t in grp["term"]) for g, grp in pairs.groupby("gene")}
    return terms, index, genes


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    go = pd.read_csv(ROOT / "Williams_82.GO.txt.gz", sep="\t").drop_duplicates()
    for level, name in NAMESPACES.items():
        sub = go[go["Level"] == level].rename(columns={"Gene": "gene", "GO": "term"})
        names = sub.drop_duplicates("term").set_index("term")["Description"].to_dict()
        terms, _, genes = encode(sub, names)
        write(f"{name}.json", {"terms": [[t, names[t]] for t in terms], "genes": genes})

    kegg = pd.read_csv(ROOT / "Williams_82.KEGG.txt.gz", sep="\t").drop_duplicates()
    pw = kegg.rename(columns={"Gene": "gene", "Pathway": "term"})[["gene", "term"]].drop_duplicates()
    names = kegg.drop_duplicates("Pathway").set_index("Pathway")["Description"].to_dict()
    terms, _, genes = encode(pw, names)
    ko = kegg.groupby("Gene")["KO"].apply(lambda s: ",".join(sorted(set(s)))).to_dict()
    write("kegg.json", {"terms": [[t, names[t], int(t.startswith(GLOBAL_MAPS))] for t in terms],
                        "genes": genes, "ko": ko})

    ids = pd.read_csv(ROOT / "data" / "zh13-wm82a4" / "Wm82a4_to_ZH13.tsv", sep="\t")
    a4 = set(ids["Wm82.a4_gene_ID"])
    a2_to_a4 = {a2: g for g, a2 in zip(ids["Wm82.a4_gene_ID"], ids["Wm82.a2.v1_ID"]) if isinstance(a2, str)}
    pf = pd.read_csv(ROOT / "Williams82.a2.v1.pfam.domain.txt", sep="\t")
    pf["gene"] = [a2_to_a4.get(g, g if g in a4 else None) for g in pf["GeneId"]]
    print(f"Pfam: {pf['GeneId'].nunique()} Wm82.a2.v1 genes, {pf['gene'].nunique()} converted to Wm82.a4")
    pf = pf.dropna(subset=["gene"]).rename(columns={"PFAM": "term"})[["gene", "term"]].drop_duplicates()
    terms, _, genes = encode(pf, {})
    write("pfam.json", {"terms": [[t, ""] for t in terms], "genes": genes})

    # Download files: GO and KEGG as provided, Pfam with Wm82.a4 IDs.
    shutil.copy(ROOT / "Williams_82.GO.txt.gz", OUT / "Wm82a4_GO.tsv.gz")
    shutil.copy(ROOT / "Williams_82.KEGG.txt.gz", OUT / "Wm82a4_KEGG.tsv.gz")
    with gzip.open(OUT / "Wm82a4_Pfam.tsv.gz", "wt", encoding="utf-8") as fh:
        pf.rename(columns={"gene": "Gene", "term": "Pfam"}).sort_values(["Gene", "Pfam"]).to_csv(fh, sep="	", index=False)

    tpm_genes = pd.read_csv(ROOT / "data" / "abiotic-rnaseq" / "abiotic_stress_TPM_matrix.tsv.gz",
                            sep="	", usecols=[0]).iloc[:, 0]
    write("genes.json", sorted(tpm_genes))
    meta = json.loads((ROOT / "data" / "abiotic-rnaseq" / "meta.json").read_text(encoding="utf-8"))
    degs = pd.read_csv(ROOT / "data" / "abiotic-rnaseq" / "abiotic_stress_DEGs.tsv.gz", sep="\t")
    degs = degs.drop_duplicates(["gene", "contrast"])
    counts = degs.groupby(["gene", "stress", "direction"]).size().unstack(["stress", "direction"], fill_value=0)
    cols = [(s, d) for s in ("Salt", "Drought", "Heat") for d in ("Up", "Down")]
    counts = counts.reindex(columns=cols, fill_value=0)
    comparisons = pd.Series([c["stress"] for c in meta["contrasts"]]).value_counts().to_dict()
    write("stress.json", {"genes": {g: [int(v) for v in row] for g, row in zip(counts.index, counts.to_numpy())},
                          "comparisons": comparisons, "universe": meta["genes"]})


if __name__ == "__main__":
    main()

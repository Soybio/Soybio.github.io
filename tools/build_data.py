#!/usr/bin/env python3
"""Build the download files and the lookup data for the gene ID converter.

Reads the Excel table in data/zh13-wm82a4/ and writes, in the same folder:
  ZH13_to_Wm82a4.tsv   - sheet "ZH13_to_Wm82.a4" as tab-separated text
  Wm82a4_to_ZH13.tsv   - sheet "Wm82.a4_to_ZH13" as tab-separated text
  lookup.json          - compact copy of both sheets used by converter.html

Run from the website folder:  python3 tools/build_data.py
"""
import json
import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "zh13-wm82a4"
XLSX = DATA_DIR / "ZH13_to_Wm82a4_gene_ID_conversion.xlsx"

ZH_SHEET = "ZH13_to_Wm82.a4"
WM_SHEET = "Wm82.a4_to_ZH13"

# Columns kept in lookup.json. Wm82.a4 location, description and Wm82.a2.v1 ID
# are stored once, in the Wm82.a4 table, because they describe the Wm82.a4 gene.
ZH_FIELDS = {
    "ZH13_gene_ID": "zh",
    "Wm82.a4_gene_ID": "wm",
    "Confidence": "conf",
    "Evidence": "evidence",
    "Liftover_Wm82.a4_gene": "liftover",
    "Protein_identity_%": "pid",
    "CDS_lifted_onto_this_gene_%": "cds",
    "Synteny_block": "block",
    "ZH13_location": "zh_loc",
    "SoyOD_converter_answer": "soyod",
    "Agreement_with_SoyOD": "agree",
}
WM_FIELDS = {
    "Wm82.a4_gene_ID": "wm",
    "ZH13_gene_ID": "zh",
    "Confidence": "conf",
    "Wm82.a2.v1_ID": "a2",
    "Wm82.a4_location": "wm_loc",
    "Wm82.a4_description": "desc",
}


def clean(value):
    """Convert pandas values to plain JSON values ("" for empty text)."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    if isinstance(value, float):
        return round(value, 1)
    return str(value)


def encode_table(df, fields, same_as_wm=()):
    """Store a sheet as rows to keep the file small.

    Text columns with many repeated values become an index into a list
    (table["dicts"]). In the columns named in same_as_wm, a value equal to the
    Wm82.a4 gene ID of the same row is stored as 1.
    """
    table = {"fields": list(fields.values()), "dicts": {}, "rows": []}
    wm_ids = [clean(v) for v in df["Wm82.a4_gene_ID"].tolist()]
    columns = []
    for col, name in fields.items():
        values = [clean(v) for v in df[col].tolist()]
        if name in same_as_wm:
            values = [1 if v is not None and v == w else v for v, w in zip(values, wm_ids)]
        is_text = all(v is None or isinstance(v, str) for v in values)
        filled = [v for v in values if v is not None]
        unique = sorted(set(filled)) if is_text else []
        if is_text and len(unique) < 0.7 * len(filled):
            index = {v: i for i, v in enumerate(unique)}
            table["dicts"][name] = unique
            values = [None if v is None else index[v] for v in values]
        elif is_text:
            values = ["" if v is None else v for v in values]
        columns.append(values)
    table["rows"] = [list(row) for row in zip(*columns)]
    return table


def main():
    sheets = pd.read_excel(XLSX, sheet_name=None)
    zh = sheets[ZH_SHEET]
    wm = sheets[WM_SHEET]

    zh.to_csv(DATA_DIR / "ZH13_to_Wm82a4.tsv", sep="\t", index=False, float_format="%g")
    wm.to_csv(DATA_DIR / "Wm82a4_to_ZH13.tsv", sep="\t", index=False, float_format="%g")

    lookup = {
        "source": XLSX.name,
        "zh": encode_table(zh, ZH_FIELDS, same_as_wm=("liftover", "soyod")),
        "wm": encode_table(wm, WM_FIELDS),
    }
    with open(DATA_DIR / "lookup.json", "w", encoding="utf-8") as fh:
        json.dump(lookup, fh, ensure_ascii=False, separators=(",", ":"))

    # Numbers used on the website pages.
    print("ZH13 genes:", len(zh))
    print("Wm82.a4 genes:", len(wm))
    print("Gene pairs:", zh["Wm82.a4_gene_ID"].notna().sum())
    print(zh["Confidence"].value_counts().to_string(), "\n")
    print(wm["Confidence"].value_counts().to_string(), "\n")
    print(zh["Agreement_with_SoyOD"].value_counts().to_string(), "\n")
    print(pd.crosstab(zh["Evidence"].str.replace(r"Glyma\.\w+|\d+(\.\d+)? ?%", "#", regex=True),
                      zh["Confidence"]).to_string())
    for f in sorted(DATA_DIR.iterdir()):
        print(f"{f.stat().st_size / 1e6:8.2f} MB  {f.name}")


if __name__ == "__main__":
    main()

# Soybio: Soybean Bioinformatics and Genomics

Source files of the website <https://soybio.github.io>.

## Folder layout

| Path | Contents |
|---|---|
| `index.html` | Home page |
| `data.html` | List of datasets |
| `data/zh13-wm82a4/` | ZH13 to Wm82.a4 gene ID conversion: dataset page and download files |
| `data/abiotic-rnaseq/` | Abiotic stress RNA-seq: dataset page, download files, and the data used by the Expression page (`meta.json`, `idmap.json`, `chunks/`) |
| `expression.html` | Expression page (bar plots and heatmaps) |
| `analysis.html` | Gene Analysis page (enrichment, co-expression) |
| `reference.html` | qPCR Reference Gene Finder page |
| `tools.html` | Tools page (sequence and qPCR tools) |
| `converter.html` | Gene ID Converter page |
| `resources.html` | Links to databases and software |
| `about.html` | About page |
| `assets/css/style.css` | Styles for all pages |
| `assets/js/converter.js` | Code of the Gene ID Converter |
| `assets/js/genedata.js` | Shared loading of the RNA-seq data and gene ID lookup |
| `assets/js/expression.js` | Code of the Expression page |
| `assets/js/reference.js` | Code of the Reference Genes page |
| `assets/js/analysis.js` | Code of the Gene Analysis page |
| `assets/js/seqdata.js` | Constants for the Tools page, generated from Biopython by `tools/build_seqdata.py` |
| `assets/js/seqtools.js` | Calculations of the Tools page |
| `assets/js/tools.js` | Forms and plots of the Tools page |
| `assets/js/export.js` | Shared CSV, SVG and PNG download helpers |
| `assets/js/site.js` | Drop-down menus and the small-screen menu button |
| `data/annotation/` | GO, KEGG and Pfam files (dataset page, downloads, data for Gene Analysis) |
| `tools/build_data.py` | Builds the TSV files and `lookup.json` from the Excel file |
| `tools/build_expression.py` | Builds the files in `data/abiotic-rnaseq/` from the RNA-seq results |
| `tools/build_annotation.py` | Builds `data/annotation/` from the GO, KEGG and Pfam files in the top folder |
| `tools/build_coexpression.py` | Builds `data/abiotic-rnaseq/coexp/` (co-expressed genes) |
| `tools/build_reference.py` | Builds `data/abiotic-rnaseq/reference_genes.json` (run after `build_expression.py`) |
| `tools/build_nav.py` | Writes the header, main menu and footer into every page. To change the menu, edit it there and run `python3 tools/build_nav.py` |

## Updating the ZH13 to Wm82.a4 table

1. Replace `data/zh13-wm82a4/ZH13_to_Wm82a4_gene_ID_conversion.xlsx` with the new file
   (same sheet and column names).
2. Run `python3 tools/build_data.py` (needs `pandas` and `openpyxl`).
3. Update the numbers, file sizes and version history in `data/zh13-wm82a4/index.html`
   and `index.html`.
4. Commit and push. GitHub Pages updates the website within a few minutes.

## Updating the RNA-seq expression data

1. Unzip the RNA-seq results (folder `newly_processed_abiotic_RNAseq` with `matrices/`,
   `downstream/tables/` and `SRR_bioproject_list.xlsx`) outside this folder.
2. Run `python3 tools/build_expression.py <path to newly_processed_abiotic_RNAseq>`.
   The TPM matrix columns must be in the same order as the rows of the SRA list.
3. Update the numbers in `data/abiotic-rnaseq/index.html`, `expression.html` and `index.html`.
4. Commit and push. The zip file itself is not uploaded (it is listed in `.gitignore`).

## Checking the Tools page

The calculations in `assets/js/seqtools.js` are compared with Biopython and SciPy on random test cases:

    python3 tools/validate_tools.py make cases.json
    node tools/validate_tools.js cases.json results.json
    python3 tools/validate_tools.py check cases.json results.json

Run this after any change to `seqtools.js` or `seqdata.js`; every line should show 0 failed.

## Adding a new dataset

1. Create a folder `data/<dataset-name>/` with the files and an `index.html`
   (copy `data/zh13-wm82a4/index.html` as a starting point).
2. Add a row for the dataset in `data.html`.

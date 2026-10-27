# Soybio – Soybean Bioinformatics and Genomics

Source files of the website <https://soybio.github.io>.

## Folder layout

| Path | Contents |
|---|---|
| `index.html` | Home page |
| `data.html` | List of datasets |
| `data/zh13-wm82a4/` | ZH13 to Wm82.a4 gene ID conversion: dataset page and download files |
| `converter.html` | Gene ID Converter page |
| `resources.html` | Links to databases and software |
| `about.html` | About page |
| `assets/css/style.css` | Styles for all pages |
| `assets/js/converter.js` | Code of the Gene ID Converter |
| `tools/build_data.py` | Builds the TSV files and `lookup.json` from the Excel file |

## Updating the ZH13 to Wm82.a4 table

1. Replace `data/zh13-wm82a4/ZH13_to_Wm82a4_gene_ID_conversion.xlsx` with the new file
   (same sheet and column names).
2. Run `python3 tools/build_data.py` (needs `pandas` and `openpyxl`).
3. Update the numbers, file sizes and version history in `data/zh13-wm82a4/index.html`
   and `index.html`.
4. Commit and push. GitHub Pages updates the website within a few minutes.

## Adding a new dataset

1. Create a folder `data/<dataset-name>/` with the files and an `index.html`
   (copy `data/zh13-wm82a4/index.html` as a starting point).
2. Add a row for the dataset in `data.html`.

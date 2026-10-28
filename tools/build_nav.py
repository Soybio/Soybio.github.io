"""Write the shared header, main menu and footer into every page.

The menu is defined once here. After changing it, run from the website folder:

    python3 tools/build_nav.py

Each page keeps its own content; only the blocks between the
<!-- site header --> and <!-- site footer --> markers are replaced.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# (label, link) for single items; (label, [entries]) for drop-down menus.
# In a drop-down, an entry starting with "#" is a group heading.
MENU = [
    ("Home", "/"),
    ("Data", [
        ("Data catalogue", "/data.html"),
        "-",
        ("ZH13 to Wm82.a4 gene IDs", "/data/zh13-wm82a4/"),
        ("Abiotic stress RNA-seq", "/data/abiotic-rnaseq/"),
        ("Functional annotation (GO, KEGG, Pfam)", "/data/annotation/"),
    ]),
    ("Analysis", [
        ("Expression browser", "/expression.html"),
        ("Gene analysis and enrichment", "/analysis.html"),
        ("Reference gene finder", "/reference.html"),
        ("Gene ID converter", "/converter.html"),
    ]),
    ("Tools", [
        ("Overview of all tools", "/tools.html"),
        "#Sequence analysis",
        ("Protein properties", "/tools.html#protein"),
        ("DNA statistics", "/tools.html#dna"),
        ("ORF finder", "/tools.html#orf"),
        ("Codon usage", "/tools.html#codon"),
        ("Pairwise alignment", "/tools.html#align"),
        "#Primers and cloning",
        ("Primer check and in-silico PCR", "/tools.html#primer"),
        ("Restriction sites", "/tools.html#restriction"),
        "#qPCR",
        ("Relative expression", "/tools.html#qpcr"),
        ("Standard curve", "/tools.html#curve"),
    ]),
    ("Resources", "/resources.html"),
    ("About", "/about.html"),
]

# Page file -> its link in the menu
PAGES = {
    "index.html": "/",
    "data.html": "/data.html",
    "data/zh13-wm82a4/index.html": "/data/zh13-wm82a4/",
    "data/abiotic-rnaseq/index.html": "/data/abiotic-rnaseq/",
    "data/annotation/index.html": "/data/annotation/",
    "expression.html": "/expression.html",
    "analysis.html": "/analysis.html",
    "reference.html": "/reference.html",
    "converter.html": "/converter.html",
    "tools.html": "/tools.html",
    "resources.html": "/resources.html",
    "about.html": "/about.html",
    "404.html": None,
}


def header(current):
    out = [
        '  <!-- site header -->',
        '  <header class="site-header">',
        '    <div class="container masthead">',
        '      <a class="brand" href="/">',
        '        <img src="/assets/img/logo.svg" alt="" width="44" height="44">',
        '        <span>',
        '          <span class="brand-name">Soybio</span>',
        '          <span class="brand-tag">Soybean Bioinformatics and Genomics</span>',
        '        </span>',
        '      </a>',
        '      <p class="masthead-meta">Genome versions: Wm82.a4, Wm82.a2.v1 and ZH13<br>Free data access, no account needed</p>',
        '    </div>',
        '  </header>',
        '  <nav class="site-nav" aria-label="Main">',
        '    <div class="container">',
        '      <button type="button" class="nav-menu-btn" aria-expanded="false" aria-controls="nav-list">Menu</button>',
        '      <ul class="nav-list" id="nav-list">',
    ]
    for label, target in MENU:
        if isinstance(target, str):
            cur = ' aria-current="page"' if target == current else ''
            out.append(f'        <li><a href="{target}"{cur}>{label}</a></li>')
            continue
        links = [e[1] for e in target if isinstance(e, tuple)]
        cls = ' current' if current in links else ''
        out.append('        <li class="has-sub">')
        out.append(f'          <button type="button" class="sub-toggle{cls}" aria-expanded="false">{label}</button>')
        out.append('          <ul class="submenu">')
        for e in target:
            if e == "-":
                out.append('            <li class="sub-sep" role="presentation"></li>')
            elif isinstance(e, str):
                out.append(f'            <li class="sub-head">{e[1:]}</li>')
            else:
                cur = ' aria-current="page"' if e[1] == current else ''
                out.append(f'            <li><a href="{e[1]}"{cur}>{e[0]}</a></li>')
        out.append('          </ul>')
        out.append('        </li>')
    out += [
        '      </ul>',
        '    </div>',
        '  </nav>',
        '  <!-- /site header -->',
    ]
    return "\n".join(out)


def footer_list(title, entries):
    rows = [f'        <p class="footer-title">{title}</p>', '        <ul>']
    for label, link in entries:
        rows.append(f'          <li><a href="{link}">{label}</a></li>')
    rows.append('        </ul>')
    return rows


def footer():
    data = [e for e in dict(MENU)["Data"] if isinstance(e, tuple)]
    analysis = dict(MENU)["Analysis"]
    site = [("Tools", "/tools.html"), ("Resources", "/resources.html"),
            ("About", "/about.html"), ("GitHub", "https://github.com/Soybio")]
    out = [
        '  <!-- site footer -->',
        '  <footer class="site-footer">',
        '    <div class="container footer-grid">',
        '      <div class="footer-about">',
        '        <p class="footer-title">Soybio</p>',
        '        <p>Data, analysis and tools for soybean genome research. All files are free to download, and the tools run in your web browser.</p>',
        '      </div>',
        '      <div>',
        *footer_list("Data", data),
        '      </div>',
        '      <div>',
        *footer_list("Analysis", analysis),
        '      </div>',
        '      <div>',
        *footer_list("Website", site),
        '      </div>',
        '    </div>',
        '    <div class="footer-bottom">',
        '      <div class="container">Soybio · Soybean Bioinformatics and Genomics · Last updated: October 2026</div>',
        '    </div>',
        '  </footer>',
        '  <!-- /site footer -->',
    ]
    return "\n".join(out)


def replace_block(html, marker, tag_pattern, new):
    marked = re.compile(rf'  <!-- {marker} -->.*?<!-- /{marker} -->', re.S)
    if marked.search(html):
        return marked.sub(lambda m: new, html, count=1)
    plain = re.compile(tag_pattern, re.S)
    if not plain.search(html):
        raise SystemExit(f"no {marker} block found")
    return plain.sub(lambda m: new, html, count=1)


def main():
    foot = footer()
    for page, current in PAGES.items():
        path = ROOT / page
        html = path.read_text(encoding="utf-8")
        html = replace_block(html, "site header", r'  <header class="site-header">.*?</header>', header(current))
        html = replace_block(html, "site footer", r'  <footer class="site-footer">.*?</footer>', foot)
        if '/assets/js/site.js' not in html:
            html = html.replace('  <!-- /site footer -->', '  <!-- /site footer -->\n\n  <script src="/assets/js/site.js"></script>', 1)
        path.write_text(html, encoding="utf-8")
        print("updated", page)


if __name__ == "__main__":
    main()

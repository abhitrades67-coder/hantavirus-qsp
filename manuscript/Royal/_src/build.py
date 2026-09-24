"""Build the Journal of the Royal Society Interface submission package.

Converts the Markdown sources to .docx with pandoc, so LaTeX maths becomes
native Word equations, copies the figures out of outputs/ under their
submission names, and reports anything missing.

The manuscript and supplement are produced double spaced with continuous line
numbers, which is what the journal asks for at review. That is done with a
patched pandoc reference document rather than by hand, so a rebuild cannot
lose the setting.

Run build_tables.py first: it writes manuscript.md and supplementary.md from
the .md.in sources with every data table filled in from outputs/.

Usage, from the repository root:
    python "manuscript/Royal/_src/build_tables.py"
    python "manuscript/Royal/_src/build.py"
"""
import io
import os
import re
import shutil
import datetime
import subprocess
import zipfile
from xml.sax.saxutils import escape
import zipfile

import pypandoc

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)                       # manuscript/Royal
REPO = os.path.dirname(os.path.dirname(PKG))      # repository root
OUTPUTS = os.path.join(REPO, "outputs")

# Document properties written into every built file, in place of the ones
# pandoc copies from the reference document.
AUTHORS = ("Abhishek Anil, Sayantan Shankar Roy, Archana A. S., "
           "Debasish Hota, Anand Srinivasan")
STAMP = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
SUPP_FIGS = os.path.join(REPO, "manuscript", "supplementary")
FIGDIR = os.path.join(PKG, "Figures")

# (source markdown, output .docx, double spaced and line numbered?)
DOCS = [
    ("title_page.md",    "Title_Page.docx",                 False),
    ("manuscript.md",    "Manuscript.docx",                 True),
    ("supplementary.md", "Supplementary_Material.docx",     True),
    ("cover_letter.md",  "Cover_Letter.docx",               False),
    ("readme.md",        "00_Submission_Guide.docx",        False),
]

FIGURES = {
    "Figure_1.png":  os.path.join(OUTPUTS, "viral_kinetics.png"),
    "Figure_2.png":  os.path.join(OUTPUTS, "treatment_window.png"),
    "Figure_3.png":  os.path.join(OUTPUTS, "organ_injury_heatmap.png"),
    "Figure_4.png":  os.path.join(OUTPUTS, "duration_start_interaction.png"),
    "Figure_5.png":  os.path.join(OUTPUTS, "identifiability_profile.png"),
    "Figure_S1.png": os.path.join(OUTPUTS, "pk_profiles.png"),
    "Figure_S2.png": os.path.join(OUTPUTS, "adaptive_immunity_heatmap.png"),
    "Figure_S3.png": os.path.join(SUPP_FIGS, "FigureS4_ribavirin_ablation.png"),
    "Figure_S4.png": os.path.join(OUTPUTS, "preexposure_viral_peak.png"),
    "Figure_S5.png": os.path.join(OUTPUTS, "gsa_prcc_tornado.png"),
    "Figure_S6.png": os.path.join(OUTPUTS, "figure_S7_external_validation.png"),
}

LN_NUM = '<w:lnNumType w:countBy="1" w:restart="continuous" w:distance="360"/>'
LINE = 'w:line="480" w:lineRule="auto"'   # 480 twentieths of a point is double
SPACED_STYLES = ("Normal", "BodyText", "FirstParagraph", "Compact")


def _force_line_spacing(pPr_block):
    if "<w:spacing" in pPr_block:
        def fix(m):
            tag = m.group(0)
            tag = re.sub(r'\s*w:line="[^"]*"', "", tag)
            tag = re.sub(r'\s*w:lineRule="[^"]*"', "", tag)
            return tag.rstrip("/>").rstrip() + " " + LINE + "/>"
        return re.sub(r"<w:spacing[^>]*/>", fix, pPr_block, count=1)
    return pPr_block.replace("<w:pPr>", "<w:pPr><w:spacing " + LINE + "/>", 1)


def set_double_spacing(styles_xml):
    """Pandoc's body text is not the Normal style: it uses Compact, BodyText
    and FirstParagraph, each setting only before/after. Setting the line height
    on Normal alone therefore does nothing visible, so it is written into the
    document defaults and into each of those styles."""
    m = re.search(r"<w:pPrDefault>.*?</w:pPrDefault>", styles_xml, re.S)
    if m:
        styles_xml = styles_xml.replace(m.group(0), _force_line_spacing(m.group(0)), 1)
    for sid in SPACED_STYLES:
        m = re.search(r'<w:style [^>]*w:styleId="%s".*?</w:style>' % sid, styles_xml, re.S)
        if not m:
            continue
        block = m.group(0)
        if "<w:pPr>" in block:
            patched = _force_line_spacing(block)
        else:
            anchor = re.search(r"<w:qFormat\s*/>|<w:name [^>]*/>", block)
            if not anchor:
                continue
            patched = block.replace(
                anchor.group(0),
                anchor.group(0) + "<w:pPr><w:spacing " + LINE + "/></w:pPr>", 1)
        styles_xml = styles_xml.replace(block, patched, 1)
    return styles_xml


def make_reference_docx(path):
    raw = pypandoc.get_pandoc_path()
    base = os.path.join(os.path.dirname(path), "_pandoc_default.docx")
    # subprocess, not os.system: the package path contains a space, which
    # os.system mangles on Windows.
    res = subprocess.run([raw, "--print-default-data-file", "reference.docx"],
                         capture_output=True)
    if res.returncode != 0 or not res.stdout:
        raise RuntimeError("could not extract pandoc's default reference.docx: %s"
                           % res.stderr.decode("utf-8", "replace")[:200])
    with open(base, "wb") as fh:
        fh.write(res.stdout)
    with zipfile.ZipFile(base) as zin:
        items = {n: zin.read(n) for n in zin.namelist()}
    doc = items["word/document.xml"].decode("utf-8")
    if "<w:lnNumType" not in doc:
        # lnNumType must precede w:cols inside w:sectPr to satisfy the schema.
        if "<w:cols" in doc:
            doc = doc.replace("<w:cols", LN_NUM + "<w:cols", 1)
        else:
            doc = doc.replace("</w:sectPr>", LN_NUM + "</w:sectPr>", 1)
    items["word/document.xml"] = doc.encode("utf-8")
    items["word/styles.xml"] = set_double_spacing(
        items["word/styles.xml"].decode("utf-8")).encode("utf-8")
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zout:
        for name, data in items.items():
            zout.writestr(name, data)
    os.remove(base)
    return path


def strip_frontmatter(text):
    if text.startswith("---"):
        end = text.find("\n---", 3)
        if end != -1:
            nl = text.find("\n", end + 1)
            return text[nl + 1:] if nl != -1 else ""
    return text


def clean_metadata(path, title, authors):
    """Replace the document properties pandoc inherits from the reference file.

    The reference document is a 2008-era template, so every built file would
    otherwise carry its word count, page count, application name and template
    name, none of which describe the file they are attached to. A submission
    system that reads those fields would be told the manuscript is 83 words
    long. Author and title are set from the paper itself; the counts are
    dropped rather than recomputed, and the word processor regenerates them.
    """
    core = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<cp:coreProperties'
        ' xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/'
        'core-properties"'
        ' xmlns:dc="http://purl.org/dc/elements/1.1/"'
        ' xmlns:dcterms="http://purl.org/dc/terms/"'
        ' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
        '<dc:title>%s</dc:title>'
        '<dc:creator>%s</dc:creator>'
        '<cp:lastModifiedBy>%s</cp:lastModifiedBy>'
        '<dcterms:created xsi:type="dcterms:W3CDTF">%s</dcterms:created>'
        '<dcterms:modified xsi:type="dcterms:W3CDTF">%s</dcterms:modified>'
        '</cp:coreProperties>'
    ) % (escape(title), escape(authors), escape(authors), STAMP, STAMP)

    app = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Properties'
        ' xmlns="http://schemas.openxmlformats.org/officeDocument/2006/'
        'extended-properties"'
        ' xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/'
        'docPropsVTypes">'
        '<DocSecurity>0</DocSecurity><ScaleCrop>false</ScaleCrop>'
        '<LinksUpToDate>false</LinksUpToDate><SharedDoc>false</SharedDoc>'
        '<HyperlinksChanged>false</HyperlinksChanged>'
        '</Properties>'
    )

    replace = {"docProps/core.xml": core, "docProps/app.xml": app}
    with zipfile.ZipFile(path) as z:
        items = [(i, z.read(i.filename)) for i in z.infolist()]
    tmp = path + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for info, data in items:
            if info.filename in replace:
                data = replace[info.filename].encode("utf-8")
            z.writestr(info.filename, data)
    os.replace(tmp, path)


def paper_title():
    """The title, taken from the manuscript itself rather than typed again."""
    path = os.path.join(HERE, "manuscript.md")
    for line in io.open(path, encoding="utf-8"):
        if line.startswith("# "):
            return line[2:].strip()
    raise RuntimeError("no title heading in manuscript.md")


def main():
    TITLE = paper_title()
    ref = make_reference_docx(os.path.join(HERE, "_reference_rsif.docx"))
    print("reference document: double spacing, continuous line numbers")

    for src, dst, spaced in DOCS:
        spath = os.path.join(HERE, src)
        if not os.path.exists(spath):
            print("  SKIP %-30s (no source)" % dst)
            continue
        text = strip_frontmatter(io.open(spath, encoding="utf-8").read())
        args = ["--standalone"]
        if spaced:
            args += ["--reference-doc=" + ref]
        out = os.path.join(PKG, dst)
        pypandoc.convert_text(text, "docx", format="markdown",
                              outputfile=out, extra_args=args)
        clean_metadata(out, TITLE, AUTHORS)
        print("  OK   %-30s %6.0f KB%s"
              % (dst, os.path.getsize(out) / 1024,
                 "  (double spaced, line numbered)" if spaced else ""))

    if not os.path.isdir(FIGDIR):
        os.makedirs(FIGDIR)
    print("figures:")
    missing = []
    for name, src in sorted(FIGURES.items(), key=lambda kv: (len(kv[0]), kv[0])):
        if not os.path.exists(src):
            missing.append((name, src))
            print("  MISSING %-14s <- %s" % (name, os.path.relpath(src, REPO)))
            continue
        shutil.copy2(src, os.path.join(FIGDIR, name))
        print("  OK      %-14s <- %-46s %6.0f KB"
              % (name, os.path.relpath(src, REPO), os.path.getsize(src) / 1024))

    if missing:
        print("\n%d figure(s) missing; run the generating script before "
              "submitting:" % len(missing))
        for name, src in missing:
            print("  %s <- %s" % (name, os.path.relpath(src, REPO)))
    else:
        print("\nAll %d figures present." % len(FIGURES))


if __name__ == "__main__":
    main()

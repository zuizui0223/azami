#!/usr/bin/env python3
"""Build editable, anonymous GEB review documents from the current frozen sources.

Run with the bundled runtime Python.  No source Markdown or statistical value is
edited.  Navigation substitutions are enumerated in REVIEW_NAVIGATION.json.
Figures are embedded only when an actual matching source file exists; absent
figures are reported on stderr and in the build audit, never fabricated.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import math
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from zipfile import ZipFile

from docx import Document
from docx.enum.section import WD_SECTION_START, WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.opc.constants import RELATIONSHIP_TYPE
from docx.shared import Inches, Pt, RGBColor
from lxml import etree
from PIL import Image


ROOT = Path(__file__).resolve().parent
FONT = "Times New Roman"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main", "m": "http://schemas.openxmlformats.org/officeDocument/2006/math"}
APPENDICES = [
    ("SI_S2_FOCAL_SPECIES_SELECTION_AUDIT.md", "Source eligibility of the two focal thistles"),
    ("SI_PHOTOGRAPHIC_TRAIT_MEASUREMENT_DEFINITIONS_20261029.md", "Photographic trait definitions and measurement provenance"),
    ("SI_V4_FULL_COVARIATION_AND_CONFIGURATION_CONTEXT.md", "Full trait covariation and configuration context"),
    ("SI_S1_BIO18_ANGLE_RAW_BH_HOLM_CONTEXT.md", "Native taxon angle and precipitation estimates"),
    ("SI_S1_BIO18_SPATIAL_SENSITIVITY_20261029.md", "Native taxon spatial sensitivity"),
    ("SI_S2_COLOUR_NATIVE_INTRO_FIXED_HOLM_CONTEXT.md", "Native and introduced colour coefficients"),
    ("SI_S2_COLOUR_SPATIAL_BLOCK_SENSITIVITY_20261029.md", "Native and introduced colour spatial sensitivity"),
    ("SI_H3_CLIMATE_MORPHOSPACE_BREADTH_20261029.md", "Climatic breadth and morphological breadth"),
    ("SI_H3_ENVIRONMENT_MODERATED_TRAIT_COVARIANCE_20261029.md", "Trait breadth and climate moderated covariance"),
    ("SI_H3_CLIMATE_ANALOGUE_CONDITIONAL_OVERLAP_20261029.md", "Local climatic analogues and floral resemblance"),
    ("SI_H3C_POSTPOSITIVE_ROBUSTNESS_20261029.md", "Climate analogue sensitivity and random sequence stability"),
]
APPENDIX_IDS = {name: i for i, (name, _) in enumerate(APPENDICES, 1)}
RESOURCE_IDS = {
    "Table_S1_full_208_BIO18_tests.csv": "R001",
    "Table_S2_all_28_RV_pairs.csv": "R002",
    "Table_S2_S2_original_native_introduced_colour_coefficients.csv": "R003",
    "SI_S1_all_80_original_raw_BH_Holm.csv": "R004",
    "SI_S1_BIO18_image_angle_both_scales.csv": "R005",
    "20261029_H3C_CLIMATE_ANALOGUE_TRAIT_OVERLAP_SOURCE_VERIFIED.json": "R006",
    "20261029_H3_CLIMATE_TRAIT_HYPERVOLUME_SOURCE_VERIFIED.json": "R007",
    "20261029_H3D_ENVIRONMENT_MODERATED_TRAIT_COVARIANCE_SOURCE_VERIFIED.json": "R008",
    "20261029_S1_BIO18_GEOGRAPHY_SOURCE_VERIFIED.json": "R009",
    "20261029_S1_BIO18_SPATIAL_BLOCK_SOURCE_VERIFIED.json": "R010",
    "20261029_S2_COLOUR_SPATIAL_BLOCK_SOURCE_VERIFIED.json": "R011",
    "S2_ALL_350_SPECIES_SOURCE_SCREEN_VERIFIED_20261021.json": "R012",
    "H3C_ROBUSTNESS_SENSITIVITY_AGGREGATE.json": "R013",
    "H3C_MONTE_CARLO_SEED_STABILITY_AGGREGATE.json": "R014",
}
CONTEXT_RESOURCE = {1: "R012", 2: "R001", 3: "R002", 4: "R004", 5: "R009", 6: "R003", 7: "R011", 8: "R007", 9: "R008", 10: "R006", 11: "R013"}
RUN_RESOURCE = {
    "40062126006": "R001", "40837107465": "R009", "40855350461": "R010",
    "40855508122": "R011", "40841007115": "R007", "40859868003": "R013",
}
FIGURE_STEMS = {
    1: "figure1", 2: "figure2", 3: "figure3_BIO18_among_taxa_image_angle",
    4: "figure4_S2_native_introduced_colour_coefficients",
    5: "figure5_H3_climatic_breadth_floral_morphospace",
}
SI_FIGURES = {
    3: [
        ("Figure_S1_S1_all_28_RV", "All 28 photographic trait pair associations in the common native Cirsium cohort. Original unadjusted and environment adjusted associations are retained; the common cohort comprises 156,683 photographs and 174 taxa. See Source R002 for the complete numerical table."),
        ("SI_module_RV_medians", "Original within domain and across domain photographic association medians in native Cirsium. These are descriptive associations from the same measured cohort, not tests of independent genetic modules."),
        ("SI_S2_trait_domain_RV", "Photographic trait association medians for native and introduced Cirsium arvense and C. vulgare strata. Strata have different photographic sampling and do not establish evolutionary changes in trait integration."),
        ("figure2_cross_functional_trait_pairs", "Three exploratory cross domain photographic trait pairs from the common native Cirsium cohort. These historically selected pairs remain supplementary descriptive results and are not the main sampling map."),
    ],
    4: [("SI_Figure_BIO18_image_angle_both_scales", "Original within taxon and among taxon BIO18 coefficients for image referenced capitulum angle in native Cirsium. Intervals are pointwise; the complete original family and the two rows are retained in Sources R004 and R005.")],
    5: [("figure3_BIO18_geographic_sensitivity", "Exploratory geographical conditioning of the native Cirsium BIO18 angle association. These pointwise model sensitivities are separate from the original 80 test Holm family. Source R009 retains the aggregate estimates.")],
    10: [("SI_H3_context_distance_overlap_breadth", "Supplementary climatic breadth, distance matching and Gaussian KDE density similarity context for native Cirsium. Original exploratory correction families and the distinction between density similarity and thresholded hypervolume intersection are retained. See Sources R006 and R007.")],
}
MATH_TEXT = {
    r"\sqrt{a^{*2}+b^{*2}}": "√(a*²+b*²)", r"4\pi A/P^2": "4πA/P²",
    r"A": "A", r"P": "P", r"10^{-6}": "10⁻⁶", r"\operatorname{atan2}(b^*,a^*)": "atan2(b*,a*)",
}
DISPLAY_EQUATION = r"\Delta D = \overline{d}_{\mathrm{all\ 64\ cell\ pairs}} - \overline{d}_{\mathrm{climate\ analogues}},"
DISPLAY_EQUATION_TEXT = "ΔD = dall 64 cell pairs − dclimate analogues,"
TABLE_LABEL_MAP = {"Table S3": "Table S1.1", "Table S4": "Table S2.1"}
TABLE_CAPTIONS = {
    (1, 1): "Source eligibility of Cirsium names recorded in both native and introduced regions.",
    (2, 1): "Eight photographic constructs and their archived definitions for Cirsium capitula.",
    (2, 2): "Archived implementation and evidence for Cirsium image measurements and ecological analyses.",
    (3, 1): "Original environment adjusted association medians within and across image domains in native Cirsium.",
    (3, 2): "Exploratory signed correlations among three photographic trait pairs in native Cirsium.",
    (4, 1): "Original native Cirsium image angle and BIO18 coefficients with unadjusted, BH and Holm evidence.",
    (5, 1): "Exploratory geographical conditioning of the native Cirsium BIO18 coefficient.",
    (6, 1): "Original colour and environment coefficients in native and introduced Cirsium arvense and C. vulgare strata.",
    (7, 1): "Original colour coefficients and spatial block sensitivity in native and introduced Cirsium arvense and C. vulgare strata.",
    (8, 1): "Original climatic breadth and trait breadth comparisons in sampled native Cirsium.",
    (9, 1): "Individual photographic trait breadth against climatic breadth in native Cirsium.",
    (9, 2): "Photographic trait covariance conditional on local climate in native Cirsium.",
    (10, 1): "Original local climate analogue comparisons of native Cirsium capitulum resemblance.",
    (11, 1): "Post-positive scenario sensitivity of climate matched Cirsium photographic distance.",
    (11, 2): "Random sequence stability of Cirsium climate analogue comparisons.",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@dataclass
class Block:
    kind: str
    text: str = ""
    rows: list[list[str]] | None = None
    level: int = 0


def blocks(text: str) -> list[Block]:
    """Parse the deliberately small Markdown vocabulary used by these sources."""
    lines = text.splitlines()
    out: list[Block] = []
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if not s:
            i += 1
            continue
        if s.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                row = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-+:?", c.replace(" ", "")) for c in row):
                    rows.append(row)
                i += 1
            if len({len(r) for r in rows}) != 1:
                raise ValueError("Inconsistent Markdown table width")
            out.append(Block("table", rows=rows))
            continue
        if s == r"\[":
            i += 1; formula = []
            while i < len(lines) and lines[i].strip() != r"\]":
                formula.append(lines[i].strip()); i += 1
            out.append(Block("equation", " ".join(formula))); i += 1
            continue
        h = re.match(r"^(#{1,6})\s+(.*)$", s)
        if h:
            out.append(Block("heading", h[2], level=len(h[1])))
            i += 1
            continue
        if s.startswith("```"):
            i += 1
            chunk = []
            while i < len(lines) and not lines[i].strip().startswith("```"):
                chunk.append(lines[i]); i += 1
            out.append(Block("code", "\n".join(chunk))); i += 1
            continue
        is_list = bool(re.match(r"^(?:[-*] |\d+\. )", s))
        chunk = [s]
        i += 1
        while i < len(lines) and lines[i].strip():
            next_s = lines[i].strip()
            if next_s.startswith(("#", "|", "```")) or re.match(r"^(?:[-*] |\d+\. )", next_s):
                break
            chunk.append(next_s); i += 1
        out.append(Block("list" if is_list else "paragraph", " ".join(chunk)))
    return out


def plain(text: str) -> str:
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1", text)
    text = re.sub(r"\\\((.*?)\\\)", lambda m: MATH_TEXT.get(m[1], m[1]), text)
    text = text.replace("**", "").replace("`", "")
    text = re.sub(r"(?<!\w)\*([^*]+)\*(?!\w)", r"\1", text)
    return text


class Navigation:
    def __init__(self):
        self.records: list[dict] = []

    def transform(self, text: str, appendix: int | None = None) -> str:
        original = text
        # These two source sentences describe superseded display placement only.
        text = text.replace(
            "the main new Figure 3 shows the regional-sensitivity forest.",
            "the geographical-sensitivity forest is Figure S5.1; main Figure 3 retains the original among-taxon coefficient.")
        text = text.replace(
            "a not-yet-built source-grounded geographic coverage map.",
            "the source-grounded geographic coverage map.")
        text = text.replace("The EAzami repository's verified nuclear 20-tip panel", "The archived verified nuclear 20-tip panel")
        def replace_link(match):
            label, target = match[1], match[2]
            name = Path(target).name
            if name in APPENDIX_IDS:
                sid = APPENDIX_IDS[name]
                if label in TABLE_LABEL_MAP:
                    return TABLE_LABEL_MAP[label]
                if re.fullmatch(r"(?:Supporting Information|.*supplement|.*SI|preserved separately)", label, re.I):
                    return f"Appendix S{sid}"
                return f"{label} (Appendix S{sid})"
            if name in RESOURCE_IDS:
                return f"{label} (Source {RESOURCE_IDS[name]})"
            if "github.com" in target.lower():
                run = re.search(r"/runs/(\d+)", target)
                rid = RUN_RESOURCE.get(run[1] if run else "", CONTEXT_RESOURCE.get(appendix or 0, "R001"))
                if run and run[1] == "40859868003":
                    rid = "R013 and Source R014"
                return f"{label} (provenance for Source {rid})"
            if not re.match(r"https?://", target):
                rid = CONTEXT_RESOURCE.get(appendix or 0, "R001")
                return f"{label} (associated provenance for Source {rid})"
            return match[0]
        text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", replace_link, text)
        for old, new in TABLE_LABEL_MAP.items():
            text = re.sub(re.escape(old) + r"(?!\d|\.\d)", new, text)
        # Filenames remain human readable, with working anonymous manifest IDs.
        for name, rid in RESOURCE_IDS.items():
            text = text.replace(f"`si/{name}`", f"`{name}` (Source {rid})")
        if text != original:
            # Do not copy author-identifying original targets into the review map.
            self.records.append({"appendix": appendix, "source_text_sha256": hashlib.sha256(original.encode()).hexdigest(), "export_text": plain(text)})
        return text


def font(run, size=None, bold=None, italic=None):
    run.font.name = FONT
    run.font.color.rgb = RGBColor(0, 0, 0)
    if size is not None: run.font.size = Pt(size)
    if bold is not None: run.bold = bold
    if italic is not None: run.italic = italic
    run._element.get_or_add_rPr().get_or_add_rFonts().set(qn("w:eastAsia"), FONT)


def inline(paragraph, text: str, size=None):
    """Write editable runs and safe external literature/photograph links."""
    tokens = re.split(r"(\[[^\]]+\]\([^)]+\)|\\\(.*?\\\)|\*\*|(?<!\w)\*|\*(?!\w)|`)", text)
    bold = italic = code = False
    for token in tokens:
        link = re.fullmatch(r"\[([^\]]+)\]\(([^)]+)\)", token)
        if link:
            relation = paragraph.part.relate_to(link[2], RELATIONSHIP_TYPE.HYPERLINK, is_external=True)
            h = OxmlElement("w:hyperlink"); h.set(qn("r:id"), relation)
            r = OxmlElement("w:r"); rp = OxmlElement("w:rPr")
            for tag, value in [("color", "000000"), ("sz", str(round((size or 11) * 2)))]:
                child = OxmlElement("w:" + tag); child.set(qn("w:val"), value); rp.append(child)
            if bold: rp.append(OxmlElement("w:b"))
            if italic: rp.append(OxmlElement("w:i"))
            r.append(rp); t = OxmlElement("w:t"); t.text = plain(link[1]); r.append(t); h.append(r); paragraph._p.append(h)
        elif token.startswith(r"\(") and token.endswith(r"\)"):
            run = paragraph.add_run(MATH_TEXT.get(token[2:-2], token[2:-2]))
            font(run, size)
        elif token == "**": bold = not bold
        elif token == "*": italic = not italic
        elif token == "`": code = not code
        elif token:
            run = paragraph.add_run(token)
            font(run, size, True if bold else None, True if italic else None)
            if code:
                # Keep long checksums wrappable without changing their contents.
                run.font.size = Pt(min(size or 11, 9.5))


def keep_with_next(p, value=True):
    p.paragraph_format.keep_with_next = value


def bookmark(p, name, number):
    start = OxmlElement("w:bookmarkStart"); start.set(qn("w:id"), str(number)); start.set(qn("w:name"), name)
    end = OxmlElement("w:bookmarkEnd"); end.set(qn("w:id"), str(number))
    p._p.insert(0, start); p._p.append(end)


def internal_link(p, text, anchor):
    h = OxmlElement("w:hyperlink"); h.set(qn("w:anchor"), anchor)
    r = OxmlElement("w:r"); rp = OxmlElement("w:rPr")
    color = OxmlElement("w:color"); color.set(qn("w:val"), "000000"); rp.append(color)
    r.append(rp); t = OxmlElement("w:t"); t.text = text; r.append(t); h.append(r); p._p.append(h)


def setup_document(main: bool):
    doc = Document()
    sec = doc.sections[0]
    sec.page_width = Inches(8.5); sec.page_height = Inches(11)
    sec.top_margin = sec.bottom_margin = Inches(0.8)
    sec.left_margin = sec.right_margin = Inches(0.9)
    sec.header_distance = sec.footer_distance = Inches(0.35)
    base = 12 if main else 11
    for name in ["Normal", "Title", "Subtitle", "Heading 1", "Heading 2", "Heading 3", "Caption", "List Paragraph", "Header", "Footer"]:
        st = doc.styles[name]
        st.font.name = FONT; st.font.color.rgb = RGBColor(0, 0, 0); st.font.size = Pt(base)
        st.paragraph_format.line_spacing = 2 if main else 1
        st.paragraph_format.space_after = Pt(5 if main else 6)
        st.paragraph_format.widow_control = True
    doc.styles["Title"].font.size = Pt(15 if main else 16)
    doc.styles["Title"].font.bold = True
    doc.styles["Title"].paragraph_format.space_after = Pt(12)
    for h in ["Heading 1", "Heading 2", "Heading 3"]:
        doc.styles[h].font.bold = True
        doc.styles[h].paragraph_format.space_before = Pt(9)
        doc.styles[h].paragraph_format.keep_with_next = True
    for h in ["Caption", "Footer"]:
        doc.styles[h].paragraph_format.line_spacing = 1
    doc.styles["Caption"].font.size = Pt(11)
    doc.styles["Caption"].font.italic = False
    doc.styles["Caption"].font.bold = False
    doc.styles["Caption"].paragraph_format.space_after = Pt(7)
    doc.core_properties.author = ""; doc.core_properties.last_modified_by = ""
    doc.core_properties.comments = ""; doc.core_properties.keywords = ""
    doc.core_properties.subject = "Anonymous review manuscript" if main else "Anonymous supporting information"
    doc.core_properties.title = "Cirsium climatic phenotypes" if main else "Cirsium supporting information"
    for element in doc.styles.element.xpath(".//w:pBdr"):
        element.getparent().remove(element)
    if main:
        ln = OxmlElement("w:lnNumType")
        for key, value in {"countBy": "1", "start": "0", "distance": "240", "restart": "continuous"}.items():
            ln.set(qn("w:" + key), value)
        sec._sectPr.append(ln)
    foot = sec.footer.paragraphs[0]; foot.alignment = WD_ALIGN_PARAGRAPH.CENTER
    field = OxmlElement("w:fldSimple"); field.set(qn("w:instr"), "PAGE")
    r = OxmlElement("w:r"); t = OxmlElement("w:t"); t.text = "1"; r.append(t); field.append(r); foot._p.append(field)
    return doc


def add_text(doc, b: Block, text: str, main=False, initial_title=False):
    if b.kind == "equation":
        if text != DISPLAY_EQUATION:
            raise ValueError("Unrecognised equation: preserve the source and add an explicit OMML conversion")
        p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        equation = OxmlElement("m:oMath")
        def mr(value, normal=False):
            r = OxmlElement("m:r")
            if normal:
                rp = OxmlElement("m:rPr"); rp.append(OxmlElement("m:nor")); r.append(rp)
            t = OxmlElement("m:t"); t.text = value; r.append(t); return r
        equation.append(mr("ΔD = "))
        def mean_d(label):
            sub = OxmlElement("m:sSub"); e = OxmlElement("m:e")
            bar = OxmlElement("m:bar"); bp = OxmlElement("m:barPr"); pos = OxmlElement("m:pos"); pos.set(qn("m:val"), "top"); bp.append(pos); bar.append(bp)
            be = OxmlElement("m:e"); be.append(mr("d")); bar.append(be); e.append(bar); sub.append(e)
            se = OxmlElement("m:sub"); se.append(mr(label, normal=True)); sub.append(se); return sub
        equation.append(mean_d("all 64 cell pairs")); equation.append(mr(" − "))
        equation.append(mean_d("climate analogues")); equation.append(mr(",")); p._p.append(equation)
        return p
    if b.kind == "heading":
        style = "Title" if initial_title else "Heading " + str(min(max(b.level - (1 if main else 0), 1), 3))
    else:
        style = "Normal"
    p = doc.add_paragraph(style=style)
    if b.kind == "list":
        p.paragraph_format.left_indent = Inches(0.15)
        p.paragraph_format.first_line_indent = Inches(-0.15)
    inline(p, text)
    if b.kind == "code": p.paragraph_format.line_spacing = 1
    return p


def table_widths(rows, width):
    n = len(rows[0])
    weights = []
    for c in range(n):
        vals = [plain(row[c]) for row in rows]
        body = vals[1:]
        numeric = all(re.fullmatch(r"[\d\s.,+−\-\[\]/%°×()–:]+", v) for v in body)
        longest_word = max((len(w) for v in vals for w in v.split()), default=3)
        avg = sum(min(len(v), 100) for v in vals) / len(vals)
        weights.append(max(7, min(26, 0.4 * avg + 0.75 * longest_word)) * (0.88 if numeric else 1))
    total = sum(weights)
    return [width * v / total for v in weights]


def add_table(doc, rows, main=False, appendix=None):
    width = 9.2 if not main and len(rows[0]) >= 7 else 6.7
    n = len(rows[0]); size = 10 if main else (9 if n >= 7 else 10)
    widths = table_widths(rows, width)
    if main and n == 5: widths = [1.45, 1.28, 1.32, 1.30, 1.35]
    if main and n == 7: widths = [1.20, 1.00, 0.78, 0.68, 1.00, 1.22, 0.82]
    if not main and n >= 7:
        size = 10
        if appendix == 1: widths = [1.5, 1.05, 1.05, 1.05, 1.15, 1.05, 1.05, 1.3]
        if appendix == 4: widths = [1.5, 0.8, 1.2, 1.7, 1.05, 1.05, 1.05, 0.85]
        if appendix == 6: widths = [0.95, 1.65, 0.8, 0.75, 1.0, 1.35, 0.9, 0.9, 0.9]
        if appendix == 7: widths = [1.5, 1.7, 0.8, 1.05, 0.75, 1.95, 1.45]
    tab = doc.add_table(rows=1, cols=n); tab.alignment = WD_TABLE_ALIGNMENT.CENTER; tab.autofit = False
    pr = tab._tbl.tblPr
    borders = OxmlElement("w:tblBorders")
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        e = OxmlElement("w:" + edge); e.set(qn("w:val"), "single"); e.set(qn("w:sz"), "4"); e.set(qn("w:color"), "D9D9D9"); borders.append(e)
    pr.append(borders)
    margins = OxmlElement("w:tblCellMar")
    for edge, value in [("top", "70"), ("bottom", "70"), ("left", "65"), ("right", "65")]:
        e = OxmlElement("w:" + edge); e.set(qn("w:w"), value); e.set(qn("w:type"), "dxa"); margins.append(e)
    pr.append(margins)
    for c, w in zip(tab.columns, widths): c.width = Inches(w)
    for ridx, row in enumerate(rows):
        cells = tab.rows[0].cells if ridx == 0 else tab.add_row().cells
        trpr = cells[0]._tc.getparent().get_or_add_trPr()
        cant = OxmlElement("w:cantSplit"); trpr.append(cant)
        if ridx == 0:
            repeat = OxmlElement("w:tblHeader"); trpr.append(repeat)
        for cidx, val in enumerate(row):
            cell = cells[cidx]; cell.width = Inches(widths[cidx]); cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            p = cell.paragraphs[0]
            p.paragraph_format.line_spacing = 1
            p.paragraph_format.space_before = p.paragraph_format.space_after = Pt(1)
            p.paragraph_format.keep_with_next = ridx == 0
            p.paragraph_format.widow_control = True
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER if cidx and len(plain(val)) < 33 else WD_ALIGN_PARAGRAPH.LEFT
            inline(p, val, size)
            if ridx == 0:
                shade = OxmlElement("w:shd"); shade.set(qn("w:fill"), "E7EDF2"); cell._tc.get_or_add_tcPr().append(shade)
                for run in p.runs: run.bold = True
    after = doc.add_paragraph(); after.paragraph_format.space_after = Pt(3); after.paragraph_format.line_spacing = 1
    return tab


def image_file(stem: str, dirs: list[Path], number=None):
    for d in dirs:
        exact = d / (stem + ".png")
        if exact.exists(): return exact
        if number in (1, 2):
            matches = sorted(d.glob(f"figure{number}_*.png"))
            matches = [p for p in matches if "cross_functional" not in p.name]
            if len(matches) == 1: return matches[0]
            if len(matches) > 1: raise ValueError(f"Ambiguous Figure {number}: {matches}")
    return None


def add_figure(doc, caption: str | list[str], path: Path | None, main=False, title=None, tall_paper="letter"):
    captions = [caption] if isinstance(caption, str) else caption
    full_caption = " ".join(captions)
    estimated_lines = max(2, math.ceil(len(plain(full_caption)) / 90)) + len(captions) - 1
    max_figure_width = 168.0 / 25.4 if main else 6.7
    page_height = 11
    if main and path is not None and tall_paper == "legal":
        with Image.open(path) as im:
            if max_figure_width * im.height / im.width + estimated_lines * 0.17 > 8.9:
                page_height = 14
    section_start = False
    if main and round(doc.sections[-1].page_height / Inches(1)) != page_height:
        sec = doc.add_section(WD_SECTION_START.NEW_PAGE)
        sec.page_width = Inches(8.5); sec.page_height = Inches(page_height)
        section_start = True
    for i, part in enumerate(captions):
        p = doc.add_paragraph(style="Caption")
        p.paragraph_format.page_break_before = i == 0 and not section_start
        p.paragraph_format.keep_with_next = path is not None or i < len(captions) - 1
        inline(p, part, 11)
    if path is None:
        return
    with Image.open(path) as im: ratio = im.height / im.width
    # Reserve caption space using a conservative single spaced line estimate.
    max_height = min(9.0 if page_height == 14 else 7.3, page_height - 1.6 - estimated_lines * 0.17 - 0.55)
    width = min(max_figure_width, max_height / ratio)
    ip = doc.add_paragraph(); ip.alignment = WD_ALIGN_PARAGRAPH.CENTER
    ip.paragraph_format.line_spacing = 1; ip.paragraph_format.space_after = Pt(0)
    pic = ip.add_run().add_picture(str(path), width=Inches(width))
    pic._inline.docPr.set("descr", plain(full_caption))
    pic._inline.docPr.set("name", title or path.stem)


def extract_docx(path: Path):
    with ZipFile(path) as z:
        root = etree.fromstring(z.read("word/document.xml"))
        paragraphs = ["".join(p.xpath(".//w:t/text() | .//m:t/text()", namespaces=NS)) for p in root.xpath(".//w:p", namespaces=NS)]
        tables = [[["".join(c.xpath(".//w:t/text()", namespaces=NS)) for c in row.xpath("./w:tc", namespaces=NS)] for row in tab.xpath("./w:tr", namespaces=NS)] for tab in root.xpath(".//w:tbl", namespaces=NS)]
        xml_text = "\n".join(z.read(n).decode("utf-8", errors="ignore") for n in z.namelist() if n.endswith((".xml", ".rels")))
        if re.search(r"github\.com/|zuizui|EAzami|/workspace/|/tmp/", xml_text, re.I):
            raise AssertionError(f"Author/repository/local-path identifier in {path.name}")
        core = etree.fromstring(z.read("docProps/core.xml"))
        props = [x.text for x in core if etree.QName(x).localname in {"creator", "lastModifiedBy"}]
        if any(props): raise AssertionError("Nonempty author properties")
    return paragraphs, tables


def displayed(text):
    if text == DISPLAY_EQUATION: return DISPLAY_EQUATION_TEXT
    return plain(text)


def verify_content(path, expected_paras, expected_tables):
    paragraphs, tables = extract_docx(path)
    available = Counter(paragraphs)
    missing = Counter(displayed(x) for x in expected_paras) - available
    if missing: raise AssertionError(f"Missing source paragraphs in {path.name}: {list(missing.items())[:3]}")
    expected = [[[displayed(c) for c in row] for row in tab] for tab in expected_tables]
    if tables != expected: raise AssertionError(f"Editable table text differs in {path.name}")
    return {"source_paragraphs_verified": len(expected_paras), "editable_tables_verified": len(tables), "table_cells_verified": sum(len(r) for t in tables for r in t), "anonymity_metadata_scan": "passed"}


def make_main(path, nav, image_dirs, tall_paper):
    source = path.read_text()
    source_blocks = blocks(source)
    doc = setup_document(True)
    paragraphs, table_rows = [], []
    deferred_tables, figure_captions = [], []
    i = 0; display_plan = False
    while i < len(source_blocks):
        b = source_blocks[i]
        if b.kind == "heading" and re.match(r"(?:Main display pieces|Main figures|Figure legends|Figures$)", b.text, re.I):
            display_plan = True; i += 1; continue
        if display_plan:
            if re.match(r"^\*\*Figure\s+\d", b.text):
                figure_captions.append([b.text])
            elif figure_captions and b.kind == "paragraph" and not re.match(r"^\*\*(?:Tables?\s|Supporting Information)", b.text):
                figure_captions[-1].append(b.text)
            i += 1; continue
        if b.kind == "paragraph" and re.match(r"^\*\*Table\s+\d", b.text) and i + 1 < len(source_blocks) and source_blocks[i+1].kind == "table":
            deferred_tables.append((b.text, source_blocks[i+1].rows)); i += 2; continue
        if b.kind == "table": raise ValueError("Main table has no adjacent caption")
        text = nav.transform(b.text)
        add_text(doc, b, text, main=True, initial_title=i == 0)
        paragraphs.append(text); i += 1
    for caption, rows in deferred_tables:
        p = doc.add_paragraph(style="Caption"); p.paragraph_format.page_break_before = True; keep_with_next(p)
        text = nav.transform(caption); inline(p, text, 11); paragraphs.append(text)
        rows = [[nav.transform(c) for c in row] for row in rows]
        add_table(doc, rows, True); table_rows.append(rows)
    missing = []
    included = []
    for caption in figure_captions:
        number = int(re.match(r"^\*\*Figure\s+(\d+)", caption[0])[1])
        fp = image_file(FIGURE_STEMS.get(number, f"figure{number}"), image_dirs, number)
        cap = [nav.transform(part) for part in caption]
        add_figure(doc, cap, fp, True, f"Figure {number}", tall_paper)
        paragraphs.extend(cap)
        if fp: included.append({"figure": number, "filename": fp.name, "sha256": digest(fp)})
        else: missing.append(number)
    return doc, paragraphs, table_rows, {"source_sha256": digest(path), "figures_included": included, "missing_figures": missing}


def make_si(nav, dirs):
    doc = setup_document(False)
    p = doc.add_paragraph("Supporting information", style="Title")
    doc.add_paragraph((ROOT / "GEB_RESEARCH_ARTICLE.md").read_text().splitlines()[0].lstrip("# "))
    doc.add_paragraph("The appendices retain the original scientific estimates, denominators and correction families. Source R001 through Source R014 refer to review_manifest.json in the accompanying aggregate review bundle. Internal execution links are represented by their associated source records; they are not hosted anonymous access links.")
    doc.add_paragraph("Contents", style="Heading 1")
    for i, (_, title) in enumerate(APPENDICES, 1):
        p = doc.add_paragraph(); p.paragraph_format.space_after = Pt(3)
        internal_link(p, f"Appendix S{i}  {title}", f"AppendixS{i}")
    paras = []; tables = []; figures = []; source_receipts = []
    for i, (name, title) in enumerate(APPENDICES, 1):
        source = ROOT / name
        p = doc.add_paragraph(f"Appendix S{i}  {title}", style="Heading 1")
        p.paragraph_format.page_break_before = True; bookmark(p, f"AppendixS{i}", i)
        tnum = 0
        source_blocks = blocks(source.read_text())
        for j, b in enumerate(source_blocks):
            if j == 0 and b.kind == "heading":
                # The original heading remains below the appendix heading.
                b = Block("heading", b.text, level=2)
            if b.kind == "table":
                tnum += 1
                wide = len(b.rows[0]) >= 7
                if wide:
                    preceding_heading = None
                    if doc.paragraphs and doc.paragraphs[-1].style.name.startswith("Heading"):
                        preceding_heading = deepcopy(doc.paragraphs[-1]._p)
                        old = doc.paragraphs[-1]._p; old.getparent().remove(old)
                    section = doc.add_section(WD_SECTION_START.NEW_PAGE)
                    section.orientation = WD_ORIENT.LANDSCAPE
                    section.page_width = Inches(11); section.page_height = Inches(8.5)
                    if preceding_heading is not None:
                        doc._body._body.insert(len(doc._body._body) - 1, preceding_heading)
                cap = doc.add_paragraph(f"Table S{i}.{tnum}. {TABLE_CAPTIONS[(i, tnum)]}", style="Caption"); keep_with_next(cap)
                rows = [[nav.transform(c, i) for c in row] for row in b.rows]
                add_table(doc, rows, appendix=i); tables.append(rows)
                if wide:
                    section = doc.add_section(WD_SECTION_START.NEW_PAGE)
                    section.orientation = WD_ORIENT.PORTRAIT
                    section.page_width = Inches(8.5); section.page_height = Inches(11)
            else:
                text = nav.transform(b.text, i)
                add_text(doc, b, text)
                paras.append(text)
        for n, (stem, cap) in enumerate(SI_FIGURES.get(i, []), 1):
            fp = image_file(stem, dirs)
            if fp:
                add_figure(doc, f"**Figure S{i}.{n}.** {cap}", fp, title=f"Figure S{i}.{n}")
                figures.append({"figure": f"S{i}.{n}", "filename": fp.name, "sha256": digest(fp)})
        source_receipts.append({"appendix": f"S{i}", "filename": name, "sha256": digest(source), "tables": tnum})
    return doc, paras, tables, {"sources": source_receipts, "figures_included": figures}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "submission")
    parser.add_argument("--asset-dir", type=Path, action="append", default=[])
    parser.add_argument("--review-manifest", type=Path)
    parser.add_argument("--tall-figure-paper", choices=["letter", "legal"], default="letter", help="Use a taller review page only when needed to keep a tall figure and all caption paragraphs together")
    parser.add_argument("--qa-report", type=Path)
    args = parser.parse_args()
    dirs = [ROOT / "figures", ROOT / "si", *args.asset_dir]
    nav = Navigation()
    main_doc, mp, mt, ma = make_main(ROOT / "GEB_RESEARCH_ARTICLE.md", nav, dirs, args.tall_figure_paper)
    si_doc, sp, st, sa = make_si(nav, dirs)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    main_path = args.output_dir / "GEB_RESEARCH_ARTICLE_REVIEW.docx"
    si_path = args.output_dir / "GEB_SUPPORTING_INFORMATION_REVIEW.docx"
    main_doc.save(main_path); si_doc.save(si_path)
    ma.update(verify_content(main_path, mp, mt)); sa.update(verify_content(si_path, sp, st))
    audit = {"main": ma, "supporting_information": sa, "visual_qa": "Must be rendered and visually inspected after the final source change"}
    manifest_path = args.review_manifest
    if manifest_path is None:
        manifests = sorted(args.output_dir.glob("aggregate_*/review_aggregate_candidate/review_manifest.json"))
        manifest_path = manifests[-1] if manifests else None
    manifest_receipt = {"filename": "review_manifest.json", "status": "not located; pair these documents with the aggregate review bundle"}
    if manifest_path is not None:
        manifest = json.loads(manifest_path.read_text())
        resources = {r["id"]: r for r in manifest["resources"]}
        for rid in RESOURCE_IDS.values():
            if rid not in resources: raise AssertionError(f"Required review resource {rid} missing from manifest")
            rec = resources[rid]; resource = manifest_path.parent / rec["path"]
            if not resource.is_file() or digest(resource) != rec["sha256"]:
                raise AssertionError(f"Review resource {rid} differs from its manifest")
        manifest_receipt = {"filename": manifest_path.name, "sha256": digest(manifest_path), "status": "located", "author_anonymity_certified": manifest.get("author_anonymity_certified", False)}
    mapping = {"main_source_sha256": ma["source_sha256"], "appendices": sa["sources"], "legacy_table_labels": TABLE_LABEL_MAP, "source_manifest": manifest_receipt, "source_resources": RESOURCE_IDS, "navigation_substitutions": nav.records, "review_figures": ma["figures_included"], "supplementary_figures": sa["figures_included"]}
    (args.output_dir / "REVIEW_NAVIGATION.json").write_text(json.dumps(mapping, ensure_ascii=False, indent=2) + "\n")
    if args.qa_report:
        args.qa_report.parent.mkdir(parents=True, exist_ok=True)
        args.qa_report.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"main": str(main_path), "si": str(si_path), "missing_figures": ma["missing_figures"], "main_tables": len(mt), "si_tables": len(st), "si_figures": len(sa["figures_included"])}, indent=2))


if __name__ == "__main__":
    main()

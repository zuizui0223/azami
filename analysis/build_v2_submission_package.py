#!/usr/bin/env python3
"""Build blinded DOCX manuscripts and a hashed v2 submission-candidate bundle."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

from docx import Document
from docx.enum.section import WD_SECTION_START
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "submission"
FIGURES = ROOT / "manuscript" / "figures" / "v2_submission"
RESULTS = ROOT / "analysis_outputs" / "v2_full27_environment_atlas_2026-08-27"

MAIN_SOURCES = [
    ROOT / "manuscript" / "00_title_abstract.md",
    ROOT / "manuscript" / "01_introduction.md",
    ROOT / "manuscript" / "02_methods.md",
    ROOT / "manuscript" / "03_results.md",
    ROOT / "manuscript" / "04_discussion.md",
    ROOT / "manuscript" / "05_conclusion_and_declarations.md",
    ROOT / "manuscript" / "06_references.md",
    ROOT / "manuscript" / "07_data_code_availability.md",
]

SUPPLEMENT_SOURCE = ROOT / "manuscript" / "supplement" / "SUPPLEMENTARY_INFORMATION.md"
TITLE = "Global continuous phenomics reveals trait- and scale-specific environmental geography of thistle capitula"

INLINE = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`|\[[^\]]+\]\([^\)]+\))")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--finalize", action="store_true", help="require PDFs and build final ZIP/manifest")
    return parser.parse_args()


def set_run_font(run, name: str = "Times New Roman", size: float = 11.5) -> None:
    run.font.name = name
    run.font.size = Pt(size)
    run._element.get_or_add_rPr().rFonts.set(qn("w:ascii"), name)
    run._element.get_or_add_rPr().rFonts.set(qn("w:hAnsi"), name)
    run._element.get_or_add_rPr().rFonts.set(qn("w:eastAsia"), name)


def add_inline(paragraph, text: str) -> None:
    position = 0
    for match in INLINE.finditer(text):
        if match.start() > position:
            set_run_font(paragraph.add_run(text[position:match.start()]))
        token = match.group(0)
        if token.startswith("**"):
            run = paragraph.add_run(token[2:-2])
            run.bold = True
        elif token.startswith("*"):
            run = paragraph.add_run(token[1:-1])
            run.italic = True
        elif token.startswith("`"):
            run = paragraph.add_run(token[1:-1])
            run.font.name = "Consolas"
            run.font.size = Pt(9.5)
        else:
            label = token[1:token.index("]")]
            run = paragraph.add_run(label)
            run.font.color.rgb = None
        set_run_font(run, run.font.name or "Times New Roman", run.font.size.pt if run.font.size else 11.5)
        position = match.end()
    if position < len(text):
        set_run_font(paragraph.add_run(text[position:]))


def style_document(document: Document) -> None:
    section = document.sections[0]
    section.top_margin = Inches(0.9)
    section.bottom_margin = Inches(0.9)
    section.left_margin = Inches(1.0)
    section.right_margin = Inches(1.0)
    for style_name, size, bold in (
        ("Normal", 11.5, False),
        ("Heading 1", 15, True),
        ("Heading 2", 13, True),
        ("Heading 3", 11.5, True),
    ):
        style = document.styles[style_name]
        style.font.name = "Times New Roman"
        style.font.size = Pt(size)
        style.font.bold = bold
        style._element.rPr.rFonts.set(qn("w:ascii"), "Times New Roman")
        style._element.rPr.rFonts.set(qn("w:hAnsi"), "Times New Roman")
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    normal = document.styles["Normal"].paragraph_format
    normal.line_spacing = 1.35
    normal.space_after = Pt(5)
    document.core_properties.author = ""
    document.core_properties.last_modified_by = ""
    document.core_properties.comments = "Blinded submission candidate; scientific HOLD retained."


def add_page_number(paragraph) -> None:
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run()
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = "PAGE"
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    run._r.extend([begin, instr, end])


def add_markdown(document: Document, text: str, *, skip_title_scaffold: bool = False) -> None:
    lines = text.splitlines()
    paragraph_buffer: list[str] = []

    def flush() -> None:
        if paragraph_buffer:
            paragraph = document.add_paragraph()
            add_inline(paragraph, " ".join(value.strip() for value in paragraph_buffer))
            paragraph_buffer.clear()

    skipped_preferred_title = False
    for raw in lines:
        line = raw.rstrip()
        if skip_title_scaffold and line in {"# Title and abstract", "## Preferred title"}:
            continue
        if skip_title_scaffold and not skipped_preferred_title and line.startswith("**Global continuous phenomics"):
            skipped_preferred_title = True
            continue
        if not line.strip():
            flush()
            continue
        heading = re.match(r"^(#{1,3})\s+(.+)$", line)
        if heading:
            flush()
            level = len(heading.group(1))
            title = heading.group(2)
            if skip_title_scaffold and title == "Abstract":
                level = 1
            elif skip_title_scaffold and level == 3:
                level = 2
            document.add_heading(title, level=level)
            continue
        if line.startswith("- "):
            flush()
            paragraph = document.add_paragraph(style="List Bullet")
            add_inline(paragraph, line[2:].strip())
            continue
        if re.match(r"^\d+\.\s+", line):
            flush()
            paragraph = document.add_paragraph(style="List Number")
            add_inline(paragraph, re.sub(r"^\d+\.\s+", "", line))
            continue
        if line.startswith("> "):
            flush()
            paragraph = document.add_paragraph()
            paragraph.paragraph_format.left_indent = Inches(0.35)
            run = paragraph.add_run(line[2:])
            run.italic = True
            set_run_font(run)
            continue
        if line.startswith("|"):
            # No active submission source relies on Markdown tables; fail closed if one leaks in.
            raise ValueError("Markdown table detected; convert it explicitly before DOCX build")
        paragraph_buffer.append(line)
    flush()


def set_image_alt_text(shape, description: str) -> None:
    doc_pr = shape._inline.docPr
    doc_pr.set("descr", description)
    doc_pr.set("title", description.split(".", 1)[0])


def add_figure(document: Document, path: Path, caption: str, alt_text: str) -> None:
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run()
    shape = run.add_picture(str(path), width=Inches(6.3))
    set_image_alt_text(shape, alt_text)
    caption_paragraph = document.add_paragraph()
    caption_paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
    caption_run = caption_paragraph.add_run(caption)
    caption_run.bold = True
    set_run_font(caption_run, size=10)


def build_main(output: Path) -> Path:
    document = Document()
    style_document(document)
    title = document.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title_run = title.add_run(TITLE)
    title_run.bold = True
    set_run_font(title_run, size=16)
    subtitle = document.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    subrun = subtitle.add_run("Blinded submission candidate — scientific HOLD retained")
    subrun.italic = True
    set_run_font(subrun, size=10)
    document.add_page_break()
    for index, source in enumerate(MAIN_SOURCES):
        add_markdown(document, source.read_text(encoding="utf-8"), skip_title_scaffold=index == 0)
    document.add_page_break()
    document.add_heading("Figures", level=1)
    figures = [
        ("Figure_1_v2_method_flow.png", "Figure 1. From public photographs to a gated continuous spatial phenotype atlas.", "Flow from public photographs through 27 registered endpoints, 22 measured endpoints, present geography and robustness gates to two candidates."),
        ("Figure_2_v2_taxon_mean_information_loss.png", "Figure 2. Visible image-phenotype variation compressed by taxon means.", "Horizontal bars compare raw and equal-replication variation below taxon means across 22 measured endpoints."),
        ("Figure_3_v2_within_among_scale_atlas.png", "Figure 3. Endpoint-specific environmental alignment at within- and among-taxon scales.", "Heat map of both-scale, within-only, among-only, neither and not-comparable endpoint-gradient rows."),
        ("Figure_4_v2_candidate_robustness.png", "Figure 4. Sequential robustness of the two final adaptive-pattern candidates under current controls.", "Two panels compare global and broad-space coefficients and minimum sampling effect-magnitude ratios for chroma-radiation and orientation-precipitation."),
    ]
    for filename, caption, alt in figures:
        add_figure(document, FIGURES / filename, caption, alt)
    add_page_number(document.sections[0].footer.paragraphs[0])
    destination = output / "Azami_Chapter1_v2_Main.docx"
    document.save(destination)
    return destination


def build_supplement(output: Path) -> Path:
    document = Document()
    style_document(document)
    title = document.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = title.add_run("Supporting Information — Azami Chapter 1 v2")
    run.bold = True
    set_run_font(run, size=16)
    add_markdown(document, SUPPLEMENT_SOURCE.read_text(encoding="utf-8"))
    document.add_page_break()
    document.add_heading("Supporting figures", level=1)
    add_figure(
        document,
        FIGURES / "Figure_2_v2_taxon_mean_information_loss.png",
        "Figure S1. Taxon-mean information loss across all 22 measured endpoints.",
        "Raw and equal-replication fractions of visible variation below taxon means for all measured endpoints.",
    )
    add_figure(
        document,
        FIGURES / "Figure_3_v2_within_among_scale_atlas.png",
        "Figure S2. Complete within- and among-taxon endpoint-gradient classification.",
        "Complete heat map retaining supported, unsupported and not-comparable endpoint-gradient rows.",
    )
    add_figure(
        document,
        FIGURES / "Figure_4_v2_candidate_robustness.png",
        "Figure S3. Sampling, broad-space and historical-placement robustness of the two candidates.",
        "Robustness summary for chroma-radiation and orientation-annual-precipitation, both passing 52 placement trees.",
    )
    add_page_number(document.sections[0].footer.paragraphs[0])
    destination = output / "Azami_Chapter1_v2_Supplement.docx"
    document.save(destination)
    return destination


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


TEXT_SUFFIXES = {".csv", ".json", ".md", ".py", ".txt"}


def package_bytes(path: Path, name: str) -> bytes:
    """Return cross-platform package bytes for tracked text and exact binary bytes."""
    if Path(name).suffix.lower() in TEXT_SUFFIXES:
        text = path.read_text(encoding="utf-8").replace("\r\n", "\n").replace("\r", "\n")
        return text.encode("utf-8")
    return path.read_bytes()


def bundle_sources(output: Path) -> list[tuple[Path, str]]:
    files: list[tuple[Path, str]] = []
    files.append((output / "README.md", "README.md"))
    for path in sorted(output.glob("Azami_Chapter1_v2_*.*")):
        if path.suffix.lower() in {".docx", ".pdf"}:
            files.append((path, path.name))
    for path in sorted(FIGURES.glob("Figure_*.*")):
        if path.suffix.lower() in {".png", ".pdf"}:
            files.append((path, f"figures/{path.name}"))
    for path in sorted(RESULTS.rglob("*")):
        if path.is_file():
            files.append((path, f"data/{path.relative_to(RESULTS).as_posix()}"))
    for path in MAIN_SOURCES + [
        SUPPLEMENT_SOURCE,
        ROOT / "manuscript/final_claims.json",
        ROOT / "manuscript/EXTERNAL_COMPLETION_GATES.md",
    ]:
        files.append((path, f"source/{path.relative_to(ROOT).as_posix()}"))
    for path in (
        ROOT / "analysis/ch1/hypothesis_recovery_status.json",
        ROOT / "analysis/ch1/pipeline.json",
        ROOT / "analysis/ch1/submission_contract.json",
        ROOT / "analysis/ch1/v2_full27_environment_atlas_contract.json",
        ROOT / "analysis/ch1/v2_full27_sampling_composition_sensitivity_contract.json",
    ):
        files.append((path, f"source/{path.relative_to(ROOT).as_posix()}"))
    for path in (
        ROOT / "analysis/run_geb_v2_full27_environment_atlas.py",
        ROOT / "analysis/run_geb_v2_full27_sampling_composition_sensitivity.py",
        ROOT / "analysis/run_geb_v2_full27_spatial_sensitivity.py",
        ROOT / "analysis/run_geb_v2_full27_historical_sensitivity.py",
        ROOT / "analysis/validate_geb_v2_full27_environment_atlas.py",
        ROOT / "analysis/validate_geb_v2_full27_sampling_composition.py",
        ROOT / "analysis/validate_geb_v2_full27_sensitivities.py",
        ROOT / "analysis/validate_v2_submission_artifacts.py",
        ROOT / "analysis/validate_manuscript_citations.py",
        ROOT / "analysis/build_v2_submission_figures.py",
    ):
        files.append((path, f"code/{path.relative_to(ROOT / 'analysis').as_posix()}"))
    return files


def write_deterministic_zip(destination: Path, files: list[tuple[Path, str]]) -> None:
    with ZipFile(destination, "w", compression=ZIP_DEFLATED, compresslevel=9) as archive:
        for source, name in sorted(files, key=lambda item: item[1]):
            info = ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            archive.writestr(info, package_bytes(source, name))


def finalize(output: Path) -> None:
    required = [
        output / "Azami_Chapter1_v2_Main.docx",
        output / "Azami_Chapter1_v2_Main.pdf",
        output / "Azami_Chapter1_v2_Supplement.docx",
        output / "Azami_Chapter1_v2_Supplement.pdf",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("Missing rendered submission files:\n- " + "\n- ".join(missing))
    files = bundle_sources(output)
    manifest = {
        "schema_version": 1,
        "release_label": "Azami Chapter 1 v2 submission candidate",
        "status": "scientific_hold_pending_independent_validation",
        "source_freeze_tag": "azami-ch1-v2-2026-08-27",
        "source_lane": "full27_full_environment_only",
        "claim_ceiling": "adaptive-pattern candidates under current sequential controls; not demonstrated adaptation or mechanism",
        "files": [
            {
                "path": name,
                "bytes": len(package_bytes(source, name)),
                "sha256": sha256_bytes(package_bytes(source, name)),
            }
            for source, name in sorted(files, key=lambda item: item[1])
        ],
    }
    manifest_path = output / "SUBMISSION_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    sums_path = output / "SHA256SUMS.txt"
    sums_path.write_text(
        "\n".join(f'{row["sha256"]}  {row["path"]}' for row in manifest["files"]) + "\n",
        encoding="utf-8",
    )
    zip_path = output / "Azami_Chapter1_v2_submission_candidate.zip"
    zip_files = files + [
        (manifest_path, "SUBMISSION_MANIFEST.json"),
        (sums_path, "SHA256SUMS.txt"),
    ]
    write_deterministic_zip(zip_path, zip_files)
    with sums_path.open("a", encoding="utf-8") as handle:
        handle.write(f"{sha256(zip_path)}  {zip_path.name}\n")
    print(json.dumps({
        "status": "ok",
        "main_docx": str(required[0]),
        "supplement_docx": str(required[2]),
        "bundle": str(zip_path),
        "bundle_sha256": sha256(zip_path),
        "n_manifest_files": len(manifest["files"]),
    }, indent=2))


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for path in MAIN_SOURCES + [SUPPLEMENT_SOURCE]:
        text = path.read_text(encoding="utf-8")
        if re.search(r"\b(TODO|TBD|PLACEHOLDER)\b", text, flags=re.IGNORECASE):
            raise SystemExit(f"Placeholder token in {path}")
    if args.finalize:
        # Preserve the already privacy-scrubbed, rendered and visually
        # inspected DOCX files. Rebuilding here would invalidate that review.
        finalize(args.output_dir)
    else:
        build_main(args.output_dir)
        build_supplement(args.output_dir)
        print(json.dumps({"status": "docx_built", "output": str(args.output_dir)}, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Convert PDF book to Markdown format.

This script extracts text from a PDF file and converts it to Markdown,
preserving structure like headings, lists, and tables where possible.

Usage:
    uv run -m sp_pharkin.tools.pdf_to_markdown <pdf_path> [--output <md_path>]

Example:
    uv run -m sp_pharkin.tools.pdf_to_markdown "book/Rowe P. - Pharmacokinetics - libgen.li.pdf" --output book.md
"""

import argparse
import sys
from pathlib import Path


def extract_pdf_text(pdf_path: str) -> str:
    """
    Extract text from PDF file.

    Uses pdfplumber for better structure preservation.
    Falls back to pypdf if pdfplumber not available.
    """
    try:
        import pdfplumber

        print(f"📖 Reading PDF: {pdf_path}")
        markdown_content = []

        with pdfplumber.open(pdf_path) as pdf:
            total_pages = len(pdf.pages)
            print(f"📄 Total pages: {total_pages}")

            for page_num, page in enumerate(pdf.pages, 1):
                if page_num % 10 == 0:
                    print(f"  Processing page {page_num}/{total_pages}...")

                # Extract text with basic structure
                text = page.extract_text()
                if text:
                    markdown_content.append(text)
                    markdown_content.append("\n---\n")  # Page separator

                # Try to extract tables
                tables = page.extract_tables()
                if tables:
                    for table in tables:
                        markdown_content.append(
                            "\n| "
                            + " | ".join(str(cell or "") for cell in table[0])
                            + " |\n"
                        )
                        markdown_content.append(
                            "|" + "|".join(["---"] * len(table[0])) + "|\n"
                        )
                        for row in table[1:]:
                            markdown_content.append(
                                "| "
                                + " | ".join(str(cell or "") for cell in row)
                                + " |\n"
                            )
                        markdown_content.append("\n")

        return "".join(markdown_content)

    except ImportError:
        print("⚠️  pdfplumber not installed, trying pypdf...")
        try:
            from pypdf import PdfReader  # type: ignore[import-not-found]

            print(f"📖 Reading PDF: {pdf_path}")
            markdown_content = []

            reader = PdfReader(pdf_path)
            total_pages = len(reader.pages)
            print(f"📄 Total pages: {total_pages}")

            for page_num, page in enumerate(reader.pages, 1):
                if page_num % 10 == 0:
                    print(f"  Processing page {page_num}/{total_pages}...")

                text = page.extract_text()
                if text:
                    markdown_content.append(text)
                    markdown_content.append("\n---\n")

            return "".join(markdown_content)

        except ImportError:
            print("❌ Error: Neither pdfplumber nor pypdf installed.")
            print("Install with: uv pip install pdfplumber")
            sys.exit(1)


def simple_markdown_conversion(text: str) -> str:
    """
    Apply basic heuristics to convert extracted text to Markdown.

    This is a simple approach - more sophisticated parsing would require
    ML or template matching for the specific book format.
    """
    lines = text.split("\n")
    markdown_lines = []

    for line in lines:
        stripped = line.strip()

        # Skip empty lines (but preserve some for readability)
        if not stripped:
            if markdown_lines and markdown_lines[-1] != "\n":
                markdown_lines.append("\n")
            continue

        # Try to detect headings (all caps, short lines)
        if len(stripped) < 80 and stripped.isupper() and len(stripped.split()) < 10:
            markdown_lines.append(f"\n## {stripped}\n")

        # Try to detect subheadings (Title Case, indented or short)
        elif (
            len(stripped) < 60
            and stripped[0].isupper()
            and not stripped.endswith(".")
            and len(stripped.split()) < 8
        ):
            markdown_lines.append(f"\n### {stripped}\n")

        # Regular paragraphs
        else:
            markdown_lines.append(stripped)

    return "\n".join(markdown_lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert PDF book to Markdown format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  uv run -m sp_pharkin.tools.pdf_to_markdown book.pdf
  uv run -m sp_pharkin.tools.pdf_to_markdown book.pdf --output output.md
  uv run -m sp_pharkin.tools.pdf_to_markdown "book/Rowe P. - Pharmacokinetics - libgen.li.pdf"
        """,
    )

    parser.add_argument("pdf_path", help="Path to PDF file")
    parser.add_argument(
        "--output", "-o", help="Output Markdown file path (default: pdf_name.md)"
    )
    parser.add_argument(
        "--raw",
        action="store_true",
        help="Output raw extracted text without Markdown conversion",
    )

    args = parser.parse_args()

    # Validate PDF exists
    pdf_path = Path(args.pdf_path)
    if not pdf_path.exists():
        print(f"❌ Error: PDF file not found: {pdf_path}")
        sys.exit(1)

    if not pdf_path.suffix.lower() == ".pdf":
        print(f"⚠️  Warning: File doesn't have .pdf extension: {pdf_path}")

    # Determine output path
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = pdf_path.with_suffix(".md")

    # Extract text
    print(f"\n🔄 Extracting text from PDF...")
    text = extract_pdf_text(str(pdf_path))

    # Convert to Markdown (unless --raw)
    if not args.raw:
        print(f"📝 Converting to Markdown...")
        text = simple_markdown_conversion(text)

    # Write output
    print(f"💾 Writing to: {output_path}")
    output_path.write_text(text, encoding="utf-8")

    # Stats
    num_lines = len(text.split("\n"))
    num_chars = len(text)
    print(f"\n✅ Done!")
    print(f"  Lines: {num_lines:,}")
    print(f"  Characters: {num_chars:,}")
    print(f"  File size: {output_path.stat().st_size / 1024:.1f} KB")


if __name__ == "__main__":
    main()

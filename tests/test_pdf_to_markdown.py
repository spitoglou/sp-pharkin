import sys
import types
from pathlib import Path

import sp_pharkin.tools.pdf_to_markdown as pdfmod


class FakePage:
    def __init__(self, text: str, tables=None):
        self._text = text
        self._tables = tables or []

    def extract_text(self):
        return self._text

    def extract_tables(self):
        return self._tables


class FakePDF:
    def __init__(self, pages):
        self.pages = pages

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def test_simple_markdown_conversion_headings_and_paragraphs():
    text = """TITLE\nSubheading Line\nParagraph one.\n\nSecond paragraph."""
    out = pdfmod.simple_markdown_conversion(text)
    assert "## TITLE" in out
    assert "### Subheading Line" in out
    assert "Paragraph one." in out
    assert "Second paragraph." in out


def test_extract_pdf_text_prefers_pdfplumber(monkeypatch):
    pages = [
        FakePage("Alpha", tables=[[["c1", "c2"], ["r1c1", "r1c2"]]]),
        FakePage("Beta"),
    ]

    fake_module = types.SimpleNamespace(open=lambda path: FakePDF(pages))
    monkeypatch.setitem(sys.modules, "pdfplumber", fake_module)

    text = pdfmod.extract_pdf_text("dummy.pdf")
    assert "Alpha" in text
    assert "Beta" in text
    assert "| c1 | c2 |" in text
    assert "---" in text

    monkeypatch.delitem(sys.modules, "pdfplumber", raising=False)


def test_extract_pdf_text_fallback_to_pypdf(monkeypatch):
    monkeypatch.delitem(sys.modules, "pdfplumber", raising=False)

    import builtins

    real_import = builtins.__import__

    def _no_pdfplumber(name, *args, **kwargs):
        if name == "pdfplumber":
            raise ImportError
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _no_pdfplumber)

    class FakeReader:
        def __init__(self, _path):
            self.pages = [FakePage("Gamma"), FakePage("Delta")]

    fake_pypdf = types.SimpleNamespace(PdfReader=FakeReader)
    monkeypatch.setitem(sys.modules, "pypdf", fake_pypdf)

    text = pdfmod.extract_pdf_text("dummy.pdf")
    assert "Gamma" in text
    assert "Delta" in text

    monkeypatch.delitem(sys.modules, "pypdf", raising=False)
    monkeypatch.setattr(builtins, "__import__", real_import)


def test_main_writes_converted_output(monkeypatch, tmp_path):
    pdf_path = tmp_path / "input.pdf"
    pdf_path.write_text("placeholder", encoding="utf-8")

    monkeypatch.setitem(sys.modules, "pdfplumber", None)
    monkeypatch.setattr(pdfmod, "extract_pdf_text", lambda path: "TITLE\nBody line")
    monkeypatch.setattr(sys, "argv", ["pdf_to_markdown", str(pdf_path)])

    pdfmod.main()

    out_path = pdf_path.with_suffix(".md")
    content = out_path.read_text(encoding="utf-8")
    assert "## TITLE" in content
    assert "Body line" in content

    monkeypatch.delitem(sys.modules, "pdfplumber", raising=False)

"""Page template, styles and helpers for the technical report (ReportLab)."""
from __future__ import annotations

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether, PageTemplate,
                                Paragraph, Spacer, Table, TableStyle)

FONT_DIR = "/usr/share/fonts/truetype/dejavu/"
for name, file in [("Serif", "DejaVuSerif.ttf"), ("Serif-B", "DejaVuSerif-Bold.ttf"),
                   ("Serif-I", "DejaVuSerif-Italic.ttf"), ("Serif-BI", "DejaVuSerif-BoldItalic.ttf"),
                   ("Sans", "DejaVuSans.ttf"), ("Sans-B", "DejaVuSans-Bold.ttf"),
                   ("Mono", "DejaVuSansMono.ttf")]:
    pdfmetrics.registerFont(TTFont(name, FONT_DIR + file))
pdfmetrics.registerFontFamily("Serif", normal="Serif", bold="Serif-B", italic="Serif-I", boldItalic="Serif-BI")
pdfmetrics.registerFontFamily("Sans", normal="Sans", bold="Sans-B", italic="Sans", boldItalic="Sans-B")

INK = colors.HexColor("#17202B")
MUTED = colors.HexColor("#5B6472")
ACCENT = colors.HexColor("#A8461A")
BLUE = colors.HexColor("#2A5788")
RULE = colors.HexColor("#CFCBBE")
BAND = colors.HexColor("#F3F1EA")

PAGE_W, PAGE_H = A4
MARGIN_X = 20 * mm
TEXT_W = PAGE_W - 2 * MARGIN_X

S = {
    "body": ParagraphStyle("body", fontName="Serif", fontSize=9.6, leading=13.6, alignment=TA_JUSTIFY,
                           textColor=INK, spaceAfter=5.5),
    "bodyl": ParagraphStyle("bodyl", fontName="Serif", fontSize=9.6, leading=13.6, alignment=TA_LEFT,
                            textColor=INK, spaceAfter=4),
    "h1": ParagraphStyle("h1", fontName="Sans-B", fontSize=13.5, leading=17, textColor=INK,
                         spaceBefore=12, spaceAfter=6, keepWithNext=1),
    "h2": ParagraphStyle("h2", fontName="Sans-B", fontSize=10.6, leading=14, textColor=INK,
                         spaceBefore=9, spaceAfter=4, keepWithNext=1),
    "eq": ParagraphStyle("eq", fontName="Serif", fontSize=9.8, leading=14, alignment=TA_CENTER,
                         textColor=INK, spaceBefore=3, spaceAfter=6),
    "cap": ParagraphStyle("cap", fontName="Sans", fontSize=8, leading=10.5, textColor=MUTED,
                          alignment=TA_LEFT, spaceBefore=3, spaceAfter=9),
    "cell": ParagraphStyle("cell", fontName="Sans", fontSize=8, leading=10.2, textColor=INK),
    "cellb": ParagraphStyle("cellb", fontName="Sans-B", fontSize=8, leading=10.2, textColor=INK),
    "cellr": ParagraphStyle("cellr", fontName="Sans", fontSize=8, leading=10.2, textColor=INK, alignment=2),
    "code": ParagraphStyle("code", fontName="Mono", fontSize=7.6, leading=10, textColor=INK,
                           backColor=BAND, borderPadding=(4, 5, 4, 5), spaceBefore=3, spaceAfter=8),
    "bullet": ParagraphStyle("bullet", fontName="Serif", fontSize=9.6, leading=13.4, textColor=INK,
                             leftIndent=13, bulletIndent=3, spaceAfter=2.5, alignment=TA_LEFT),
    "ref": ParagraphStyle("ref", fontName="Serif", fontSize=8.6, leading=11.6, textColor=INK,
                          leftIndent=16, firstLineIndent=-16, spaceAfter=3),
    "title": ParagraphStyle("title", fontName="Sans-B", fontSize=25, leading=31, textColor=INK, alignment=TA_LEFT),
    "subtitle": ParagraphStyle("subtitle", fontName="Sans", fontSize=12.5, leading=17, textColor=MUTED, alignment=TA_LEFT),
    "kpi_v": ParagraphStyle("kpi_v", fontName="Sans-B", fontSize=15, leading=18, textColor=INK, alignment=TA_CENTER),
    "kpi_l": ParagraphStyle("kpi_l", fontName="Sans", fontSize=7.4, leading=9.4, textColor=MUTED, alignment=TA_CENTER),
}


class Report(BaseDocTemplate):
    def __init__(self, path, title, short_title, author, date):
        super().__init__(path, pagesize=A4, leftMargin=MARGIN_X, rightMargin=MARGIN_X,
                         topMargin=22 * mm, bottomMargin=20 * mm, title=title, author=author,
                         subject="Moreno cross-guide directional coupler, WR-22, 40 GHz",
                         creator="moreno-coupler/docs/report/build_report.py")
        self.short_title, self.author, self.date = short_title, author, date
        frame = Frame(MARGIN_X, 20 * mm, TEXT_W, PAGE_H - 42 * mm, id="f", leftPadding=0, rightPadding=0,
                      topPadding=0, bottomPadding=0)
        self.addPageTemplates([PageTemplate("cover", [frame], onPage=self._cover),
                               PageTemplate("body", [frame], onPage=self._page)])
        self._heading_count = 0

    def _cover(self, c, doc):
        c.saveState()
        c.setFillColor(INK)
        c.rect(0, PAGE_H - 9 * mm, PAGE_W, 9 * mm, stroke=0, fill=1)
        c.setFillColor(ACCENT)
        c.rect(0, PAGE_H - 10.2 * mm, PAGE_W, 1.2 * mm, stroke=0, fill=1)
        self._footer(c, doc)
        c.restoreState()

    def _page(self, c, doc):
        c.saveState()
        c.setFont("Sans", 7.6)
        c.setFillColor(MUTED)
        c.drawString(MARGIN_X, PAGE_H - 13 * mm, self.short_title)
        c.drawRightString(PAGE_W - MARGIN_X, PAGE_H - 13 * mm, self.date)
        c.setStrokeColor(RULE)
        c.setLineWidth(0.6)
        c.line(MARGIN_X, PAGE_H - 15 * mm, PAGE_W - MARGIN_X, PAGE_H - 15 * mm)
        self._footer(c, doc)
        c.restoreState()

    def _footer(self, c, doc):
        c.setFont("Sans", 7.6)
        c.setFillColor(MUTED)
        c.setStrokeColor(RULE)
        c.line(MARGIN_X, 14 * mm, PAGE_W - MARGIN_X, 14 * mm)
        c.drawString(MARGIN_X, 10 * mm, self.author)
        c.drawRightString(PAGE_W - MARGIN_X, 10 * mm, "Page %d" % doc.page)

    def afterFlowable(self, f):
        if isinstance(f, Paragraph) and f.style.name in ("h1", "h2"):
            level = 0 if f.style.name == "h1" else 1
            text = f.getPlainText()
            key = "h%d" % self._heading_count
            self._heading_count += 1
            self.canv.bookmarkPage(key)
            self.canv.addOutlineEntry(text, key, level=level, closed=False)


def P(text, style="body"):
    if style == "h1":
        from reportlab.platypus import CondPageBreak
        # never leave a section heading at the foot of a page
        return [CondPageBreak(110), Paragraph(text, S[style])]
    if style == "h2":
        from reportlab.platypus import CondPageBreak
        return [CondPageBreak(80), Paragraph(text, S[style])]
    return Paragraph(text, S[style])


def bullets(items, style="bullet"):
    return [Paragraph(t, S[style], bulletText="•") for t in items]


def table(rows, widths, header=True, align_right_cols=(), zebra=True, bold_first_col=False):
    data = []
    for i, row in enumerate(rows):
        out = []
        for j, cell in enumerate(row):
            if isinstance(cell, str):
                st = "cellb" if (header and i == 0) or (bold_first_col and j == 0) else (
                    "cellr" if j in align_right_cols else "cell")
                cell = Paragraph(cell, S[st])
            out.append(cell)
        data.append(out)
    t = Table(data, colWidths=[w * TEXT_W for w in widths], repeatRows=1 if header else 0)
    style = [("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
             ("TOPPADDING", (0, 0), (-1, -1), 3.2), ("BOTTOMPADDING", (0, 0), (-1, -1), 3.2),
             ("LEFTPADDING", (0, 0), (-1, -1), 5), ("RIGHTPADDING", (0, 0), (-1, -1), 5),
             ("LINEBELOW", (0, -1), (-1, -1), 0.8, INK)]
    if header:
        style += [("LINEABOVE", (0, 0), (-1, 0), 0.8, INK), ("LINEBELOW", (0, 0), (-1, 0), 0.6, INK)]
    else:
        style += [("LINEABOVE", (0, 0), (-1, 0), 0.8, INK)]
    if zebra:
        for i in range(1 if header else 0, len(rows)):
            if (i - (1 if header else 0)) % 2 == 1:
                style.append(("BACKGROUND", (0, i), (-1, i), BAND))
    t.setStyle(TableStyle(style))
    return t


def figure(path, caption, width_frac=1.0):
    from PIL import Image as PILImage

    w, h = PILImage.open(path).size
    width = TEXT_W * width_frac
    img = Image(str(path), width=width, height=width * h / w)
    return KeepTogether([img, P(caption, "cap")])


def kpis(items):
    cells = [[P(v, "kpi_v") for v, _ in items], [P(lbl, "kpi_l") for _, lbl in items]]
    t = Table(cells, colWidths=[TEXT_W / len(items)] * len(items))
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), BAND),
                           ("LINEABOVE", (0, 0), (-1, 0), 1.2, ACCENT),
                           ("TOPPADDING", (0, 0), (-1, 0), 8), ("BOTTOMPADDING", (0, -1), (-1, -1), 8),
                           ("LINEAFTER", (0, 0), (-2, -1), 0.5, RULE)]))
    return t


def spacer(h=6):
    return Spacer(1, h)

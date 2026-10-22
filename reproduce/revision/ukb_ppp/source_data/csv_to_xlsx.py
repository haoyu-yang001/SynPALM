"""Assemble the Source Data CSV files (one per sheet) into one Excel workbook.

    python3 csv_to_xlsx.py <csv_dir> <sheet_order.txt> <out.xlsx>

Called by build_source_data.R. Numbers are written as numbers, TRUE/FALSE as booleans,
empty fields as empty cells; NaN / Inf (undefined ratios) are kept as text.
"""
import csv
import sys

from openpyxl import Workbook
from openpyxl.cell import WriteOnlyCell
from openpyxl.styles import Alignment, Font

csv_dir, order_file, out = sys.argv[1:4]
sheets = [s.strip() for s in open(order_file) if s.strip()]
csv.field_size_limit(sys.maxsize)
INT_RE_OK = set("+-0123456789")


def conv(v):
    if v == "" or v == "NA":
        return None
    if v in ("TRUE", "FALSE"):
        return v == "TRUE"
    if v[0] in INT_RE_OK or v[0] == ".":
        try:
            return int(v) if v.lstrip("+-").isdigit() and len(v) < 16 else _float(v)
        except ValueError:
            return v
    return v


def _float(v):
    f = float(v)
    if f != f or f in (float("inf"), float("-inf")):
        return v
    return f


wb = Workbook(write_only=True)
bold = Font(bold=True)
wrap = Alignment(wrap_text=True, vertical="top")
for name in sheets:
    ws = wb.create_sheet(title=name)
    with open(f"{csv_dir}/{name}.csv", newline="", encoding="utf-8") as fh:
        rd = csv.reader(fh)
        header = next(rd)
        if name == "README":
            for col, w in zip("ABCDE", (22, 34, 90, 9, 110)):
                ws.column_dimensions[col].width = w
        else:
            ws.freeze_panes = "A2"
            for i, h in enumerate(header):
                ws.column_dimensions[chr(65 + i) if i < 26 else "A" + chr(65 + i - 26)].width = max(10, min(40, len(h) + 2))
        hdr = []
        for h in header:
            c = WriteOnlyCell(ws, value=h)
            c.font = bold
            hdr.append(c)
        ws.append(hdr)
        for row in rd:
            if name == "README":
                cells = []
                for i, v in enumerate(row):
                    c = WriteOnlyCell(ws, value=conv(v) if i == 3 else (v if v != "" else None))
                    c.alignment = wrap
                    if i == 0:
                        c.font = bold
                    cells.append(c)
                ws.append(cells)
            else:
                ws.append([conv(v) for v in row])
wb.save(out)
print("wrote", out, "with", len(sheets), "sheets")

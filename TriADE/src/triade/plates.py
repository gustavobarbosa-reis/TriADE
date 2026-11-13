"""Plate formats (96, 384 and 1536 wells), well names and ejection orders."""
from __future__ import annotations

import re
import string
from dataclasses import dataclass, field


def _row_labels(n: int) -> list[str]:
    """A, B, ..., Z, AA, AB, ... (1536-well plates use rows A–AF)."""
    letters = string.ascii_uppercase
    out = []
    for i in range(n):
        out.append(letters[i] if i < 26 else letters[i // 26 - 1] + letters[i % 26])
    return out


@dataclass(frozen=True)
class PlateFormat:
    name: str
    n_rows: int
    n_cols: int
    rows: tuple = field(init=False)

    def __post_init__(self):
        object.__setattr__(self, "rows", tuple(_row_labels(self.n_rows)))

    @property
    def n_wells(self) -> int:
        return self.n_rows * self.n_cols

    def wells(self) -> list[str]:
        """Natural plate order: A1, A2, ..., then the next row."""
        return [f"{r}{c}" for r in self.rows for c in range(1, self.n_cols + 1)]

    def row_index(self, row: str) -> int:
        return self.rows.index(row)

    def parse(self, text) -> str | None:
        """'A1', 'a01', 'A:1', ' AF 48 ' -> canonical well name, or None if it is not a well of this plate."""
        m = re.fullmatch(r"\s*([A-Za-z]{1,2})\s*:?\s*0*(\d{1,3})\s*", str(text))
        if not m:
            return None
        row, col = m.group(1).upper(), int(m.group(2))
        if row not in self.rows or not 1 <= col <= self.n_cols:
            return None
        return f"{row}{col}"

    def split(self, well: str) -> tuple[str, int]:
        """'AB12' -> ('AB', 12)."""
        m = re.fullmatch(r"([A-Z]{1,2})(\d{1,3})", well)
        if not m:
            raise ValueError(f"not a well name: {well!r}")
        return m.group(1), int(m.group(2))

    def position(self, well: str) -> tuple[int, int]:
        """Zero-based (row index, column index)."""
        r, c = self.split(well)
        return self.rows.index(r), c - 1

    def ejection_orders(self) -> dict[str, list[str]]:
        """Common Echo ejection patterns over the whole plate."""
        R, C = list(self.rows), list(range(1, self.n_cols + 1))
        first_row, last_row, last_col = R[0], R[-1], C[-1]
        return {
            f"row-serpentine (A1→A{last_col}, B{last_col}→B1, ...)":
                [f"{r}{c}" for i, r in enumerate(R) for c in (C if i % 2 == 0 else C[::-1])],
            f"row-major (A1→A{last_col}, B1→B{last_col}, ...)": [f"{r}{c}" for r in R for c in C],
            f"column-serpentine ({first_row}1→{last_row}1, {last_row}2→{first_row}2, ...)":
                [f"{r}{c}" for j, c in enumerate(C) for r in (R if j % 2 == 0 else R[::-1])],
            f"column-major ({first_row}1→{last_row}1, {first_row}2→{last_row}2, ...)": [f"{r}{c}" for c in C for r in R],
        }


PLATE_FORMATS = {
    "96": PlateFormat("96", 8, 12),
    "384": PlateFormat("384", 16, 24),
    "1536": PlateFormat("1536", 32, 48),
}


def get_plate_format(name) -> PlateFormat:
    key = str(name).strip().replace("-well", "").replace(" wells", "")
    if key not in PLATE_FORMATS:
        raise ValueError(f"unknown plate format {name!r}; choose one of {', '.join(PLATE_FORMATS)}")
    return PLATE_FORMATS[key]

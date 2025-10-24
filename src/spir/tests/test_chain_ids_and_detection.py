from __future__ import annotations

from ..detect import detect_from_pathnames
from ..utils.ids import spreadsheet_ids


def test_spreadsheet_ids():
    assert spreadsheet_ids(1) == ["A"]
    assert spreadsheet_ids(3) == ["A", "B", "C"]
    assert spreadsheet_ids(28)[-2:] == ["AA", "AB"]
    assert spreadsheet_ids(29)[-2:] == ["AB", "AC"]

from __future__ import annotations

from typing import List


def spreadsheet_ids(count: int) -> List[str]:
    return [_index_to_letters(i) for i in range(count)]


def _index_to_letters(i: int) -> str:
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    s = ""
    i += 1
    while i > 0:
        i, r = divmod(i - 1, 26)
        s = letters[r] + s
    return s

from __future__ import annotations

from typing import Iterable, List, Tuple


def read_fasta(path: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header = None
    seq_lines: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        records.append((header, "".join(seq_lines)))
    return records


def write_fasta(path: str, records: Iterable[Tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            f.write(f"{seq}\n")

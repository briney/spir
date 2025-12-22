from __future__ import annotations

from dataclasses import dataclass

from foldir.dialects import get_dialect
from foldir.ir.models import DocumentIR
from foldir.ir.normalize import normalize_document


@dataclass(frozen=True)
class ConvertOptions:
    default_seed: int = 1
    default_glycan_parent_atom: str = "O4"
    default_glycan_child_atom: str = "C1"
    default_asn_atom: str = "ND2"
    default_ser_atom: str = "OG"
    default_thr_atom: str = "OG1"


def convert(
    in_path: str,
    in_dialect: str,
    out_path: str,
    out_dialect: str,
    opts: ConvertOptions,
) -> None:
    src = get_dialect(in_dialect)
    dst = get_dialect(out_dialect)

    doc = src.parse(in_path)
    doc = normalize_document(doc, opts=opts)
    dst.render(doc, out_path)

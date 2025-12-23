from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional


from spir.dialects import get_dialect
from spir.ir.models import DocumentIR
from spir.ir.normalize import normalize_document


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
    out_prefix: str,
    out_dialect: str,
    opts: ConvertOptions,
    restraints_path: Optional[str] = None,
) -> None:
    src = get_dialect(in_dialect)
    dst = get_dialect(out_dialect)

    if restraints_path:
        if in_dialect.lower() != "chai1":
            raise ValueError("--restraints is only supported for chai1 input.")
        doc = src.parse(in_path, restraints_path)
    else:
        doc = src.parse(in_path)
    doc = normalize_document(doc, opts=opts)
    render_target, output_paths = _resolve_output_paths(out_prefix, out_dialect)
    for path in output_paths:
        _ensure_parent_dir(path)
    dst.render(doc, render_target)


_OUTPUT_EXTENSIONS = {
    "alphafold3": ".json",
    "alphafold3server": ".json",
    "alphafoldserver": ".json",
    "boltz2": ".yaml",
    "protenix": ".json",
}

_MULTI_OUTPUT_EXTENSIONS = {
    "chai1": (".fasta", ".constraints.csv"),
}


def _resolve_output_paths(out_prefix: str, out_dialect: str) -> tuple[str, list[str]]:
    key = out_dialect.lower()
    if key in _MULTI_OUTPUT_EXTENSIONS:
        exts = _MULTI_OUTPUT_EXTENSIONS[key]
        _ensure_no_extension(out_prefix, exts, out_dialect)
        return out_prefix, [out_prefix + ext for ext in exts]
    ext = _OUTPUT_EXTENSIONS.get(key)
    if not ext:
        raise ValueError(f"Unknown output dialect: {out_dialect}")
    _ensure_no_extension(out_prefix, (ext,), out_dialect)
    out_path = out_prefix + ext
    return out_path, [out_path]


def _ensure_no_extension(out_prefix: str, exts: tuple[str, ...], out_dialect: str) -> None:
    for ext in exts:
        if out_prefix.endswith(ext):
            raise ValueError(
                f"Output prefix for {out_dialect} should not include '{ext}'. "
                f"Pass a prefix without an extension."
            )


def _ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)

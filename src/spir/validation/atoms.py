from __future__ import annotations

from typing import Iterable

from ..ir.model import Bond


def validate_bonds_have_atoms(bonds: Iterable[Bond]) -> None:
    for b in bonds:
        if b.atom1.atom_name is None and b.atom1.atom_index is None:
            # Allow unspecified in intermediate IR, but most writers need at least a name or index
            continue
        if b.atom2.atom_name is None and b.atom2.atom_index is None:
            continue

from __future__ import annotations

from typing import Dict, List, Tuple

from ..ir.model import AtomRef, Bond, ComplexInput, Ligand


def merge_duplicate_bonds(bonds: List[Bond]) -> List[Bond]:
    seen = set()
    out: List[Bond] = []
    for b in bonds:
        key = (
            b.atom1.chain_id,
            b.atom1.residue_index,
            b.atom1.component_index,
            b.atom1.atom_name,
            b.atom1.atom_index,
            b.atom2.chain_id,
            b.atom2.residue_index,
            b.atom2.component_index,
            b.atom2.atom_name,
            b.atom2.atom_index,
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(b)
    return out


def normalize_multi_ccd_for_boltz(ci: ComplexInput) -> ComplexInput:
    """Transform multi-CCD glycans into per-component ligands with rewired bonds.

    Boltz constraints address atoms by [CHAIN_ID, RES_IDX, ATOM_NAME], where
    ligands are typically single-CCD with RES_IDX=1. This normalizer:
      - For each ligand with multiple CCDs, creates N new ligands (one per CCD)
        with ids like G1..GN.
      - Rewrites intra-glycan bonds to reference these new ligand ids with
        residue_index=1 on both sides and the same atom names.
      - Rewrites any protein↔glycan anchor to target the new root glycan id "G1"
        at @C1 (residue_index=1).
    """
    # Build mapping for multi-CCD ligands
    new_ligands: List[Ligand] = []
    split_map: Dict[str, List[str]] = {}

    def is_multi_ccd(l: Ligand) -> bool:
        return l.ccd_codes is not None and len(l.ccd_codes) > 1

    for l in ci.ligands:
        if not is_multi_ccd(l):
            new_ligands.append(l)
            continue
        # Derive a stable base from the first id
        base = l.ids[0]
        ids_for_components: List[str] = []
        for i, ccd in enumerate(l.ccd_codes or []):
            new_id = f"{base}{i+1}"
            ids_for_components.append(new_id)
            new_ligands.append(Ligand(ids=[new_id], ccd_codes=[ccd]))
        split_map[base] = ids_for_components

    # If no splits, return original ci
    if not split_map:
        return ci

    # Helper: find original ligand base for a chain id
    base_by_chain: Dict[str, str] = {}
    for l in ci.ligands:
        if l.ccd_codes is None:
            continue
        base = l.ids[0]
        base_by_chain[base] = base

    # Rewire bonds
    new_bonds: List[Bond] = []
    for b in ci.bonds:
        a1 = b.atom1
        a2 = b.atom2

        def rewire(a: AtomRef) -> AtomRef:
            # If this chain was split, redirect to the appropriate component id
            chain = a.chain_id
            # chain is exactly the base id for multi-CCD ligands we created
            if chain in split_map:
                # Choose component index (default to 1 if absent)
                comp_idx = a.component_index or a.residue_index or 1
                comp_ids = split_map[chain]
                comp_idx = max(1, min(comp_idx, len(comp_ids)))
                new_chain = comp_ids[comp_idx - 1]
                return AtomRef(
                    chain_id=new_chain,
                    residue_index=1,
                    atom_name=a.atom_name,
                )
            # Unchanged
            return a

        ra1 = rewire(a1)
        ra2 = rewire(a2)
        new_bonds.append(Bond(atom1=ra1, atom2=ra2))

    # Build transformed ComplexInput
    return ComplexInput(
        name=ci.name,
        seeds=ci.seeds,
        proteins=ci.proteins,
        rnas=ci.rnas,
        dnas=ci.dnas,
        ligands=new_ligands,
        ions=ci.ions,
        bonds=merge_duplicate_bonds(new_bonds),
        user_ccd=ci.user_ccd,
    )

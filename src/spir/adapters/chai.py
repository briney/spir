from __future__ import annotations

import csv
from typing import Dict, List, Optional, Tuple

from ..ir.glycan import (
    GlycanBond,
    GlycanComponent,
    GlycanGraph,
    parse_chai_glycan,
    render_chai_glycan,
)
from ..ir.model import AtomRef, Bond, ComplexInput, Ligand, PolymerChain
from ..utils.ids import spreadsheet_ids


class ChaiReaderWriter:
    @staticmethod
    def load_fasta_and_restraints(
        fasta_text: str, restraints_csv_text: Optional[str] = None
    ) -> ComplexInput:
        entries = _parse_fasta(fasta_text)
        proteins: List[PolymerChain] = []
        ligands: List[Ligand] = []
        bonds: List[Bond] = []
        # Assign chain ids A, B, ... in order
        chain_ids = spreadsheet_ids(len(entries))
        for idx, (etype, name, payload) in enumerate(entries):
            cid = chain_ids[idx]
            if etype == "protein":
                proteins.append(
                    PolymerChain(type="protein", ids=[cid], sequence=payload)
                )
            elif etype == "glycan":
                graph = parse_chai_glycan(payload)
                ligands.append(
                    Ligand(ids=[cid], ccd_codes=graph.components_as_ccd_list())
                )
                # Intra-glycan bonds captured via bonds connecting components; here we cannot map component indices to separate chains
                # Representation: treat the multi-CCD ligand as one chain and use component_index in bonds when anchor is defined in CSV
            else:
                # Chai supports smiles ligands in examples; support if present
                ligands.append(Ligand(ids=[cid], smiles=payload))

        if restraints_csv_text:
            bonds.extend(_parse_restraints_csv(restraints_csv_text, chain_ids))

        return ComplexInput(proteins=proteins, ligands=ligands, bonds=bonds)

    @staticmethod
    def dump(
        ci: ComplexInput, include_restraints: bool = True
    ) -> Tuple[str, Optional[str]]:
        # Emit proteins as protein entries; ligands with ccd_codes as glycan entries, others as ligand entries
        fasta_lines: List[str] = []
        chain_ids: List[str] = []
        # Map protein chain id to its sequence for residue-letter lookup
        protein_sequence_by_chain: Dict[str, str] = {}

        def add_entry(header: str, seq: str):
            fasta_lines.append(">" + header)
            fasta_lines.append(seq)

        for p in ci.proteins:
            for cid in p.ids:
                add_entry(f"protein|{cid}", p.sequence)
                chain_ids.append(cid)
                protein_sequence_by_chain[cid] = p.sequence
        for l in ci.ligands:
            for cid in l.ids:
                if l.ccd_codes is not None:
                    # Render full glycan string when multi-CCD using intra-glycan bonds
                    if len(l.ccd_codes) == 1:
                        spec = l.ccd_codes[0]
                    else:
                        # Build a GlycanGraph from CCD codes and bonds scoped to this ligand id
                        components = [GlycanComponent(ccd=c) for c in l.ccd_codes]
                        gbonds: List[GlycanBond] = []
                        for b in ci.bonds:
                            # Intra-glycan bond: both atoms reference this ligand and component indices exist
                            if (
                                b.atom1.chain_id == cid
                                and b.atom2.chain_id == cid
                                and b.atom1.component_index is not None
                                and b.atom2.component_index is not None
                            ):
                                # Determine parent (O*) and child (C1)
                                a1 = (b.atom1.component_index, b.atom1.atom_name or "")
                                a2 = (b.atom2.component_index, b.atom2.atom_name or "")
                                if a1[1].startswith("O") and a2[1].startswith("C"):
                                    gbonds.append(
                                        GlycanBond(
                                            parent_index=a1[0],
                                            parent_atom=a1[1],
                                            child_index=a2[0],
                                            child_atom=a2[1],
                                        )
                                    )
                                elif a2[1].startswith("O") and a1[1].startswith("C"):
                                    gbonds.append(
                                        GlycanBond(
                                            parent_index=a2[0],
                                            parent_atom=a2[1],
                                            child_index=a1[0],
                                            child_atom=a1[1],
                                        )
                                    )
                        if gbonds:
                            graph = GlycanGraph(components=components, bonds=gbonds)
                            spec = render_chai_glycan(graph)
                        else:
                            # Fallbacks: prefer explicit AF3 residues if preserved; otherwise single CCD
                            spec = l.af3_residues or l.ccd_codes[0]
                    add_entry(f"glycan|{cid}", spec)
                else:
                    add_entry(f"ligand|{cid}", l.smiles or "")
                chain_ids.append(cid)

        fasta_text = "\n".join(fasta_lines) + ("\n" if fasta_lines else "")

        restraints_text: Optional[str] = None
        if include_restraints and ci.bonds:
            # Only include the protein↔glycan root (C1) anchor bonds
            protein_ids = {cid for p in ci.proteins for cid in p.ids}
            glycan_ids = {
                cid for l in ci.ligands if l.ccd_codes is not None for cid in l.ids
            }
            rows: List[List[str]] = []
            header = [
                "chainA",
                "res_idxA",
                "chainB",
                "res_idxB",
                "connection_type",
                "confidence",
                "min_distance_angstrom",
                "max_distance_angstrom",
                "comment",
                "restraint_id",
            ]
            rows.append(header)
            seen: set = set()
            count = 0
            for b in ci.bonds:
                # Identify protein and glycan sides, enforce orientation: protein first
                if b.atom1.chain_id in protein_ids and b.atom2.chain_id in glycan_ids:
                    p_atom, g_atom = b.atom1, b.atom2
                elif b.atom2.chain_id in protein_ids and b.atom1.chain_id in glycan_ids:
                    p_atom, g_atom = b.atom2, b.atom1
                else:
                    continue
                # Root glycan anchor must target component 1; prefer explicit @C1
                if not (
                    (g_atom.component_index == 1)
                    or (g_atom.component_index is None and g_atom.residue_index == 1)
                ):
                    continue
                key = (p_atom.chain_id, p_atom.residue_index, g_atom.chain_id)
                if key in seen:
                    continue
                seen.add(key)
                count += 1
                # Format as ResidueLetterPosition@SingleAtomLetter (e.g., N192@N)
                resA = _format_protein_anchor(p_atom, protein_sequence_by_chain)
                # Force glycan side to @C1 for clarity
                resB = "@C1"
                rows.append(
                    [
                        p_atom.chain_id,
                        resA,
                        g_atom.chain_id,
                        resB,
                        "covalent",
                        "1.0",
                        "0.0",
                        "0.0",
                        "protein-glycan",
                        f"bond{count}",
                    ]
                )
            if len(rows) > 1:
                restraints_text = "\n".join([",".join(r) for r in rows]) + "\n"

        return fasta_text, restraints_text


def _format_res_idx(a: AtomRef) -> str:
    if a.residue_index is not None and a.atom_name is not None:
        return f"{a.residue_index}@{a.atom_name}"
    if a.residue_index is not None:
        return str(a.residue_index)
    if a.atom_name is not None:
        return f"@{a.atom_name}"
    return ""


def _format_protein_anchor(a: AtomRef, seq_by_chain: Dict[str, str]) -> str:
    """Format protein-side glycan anchor for Chai restraints.

    Expected: ResidueLetterPosition@SingleAtomLetter, e.g., N192@N or T10@O.
    Falls back to legacy formatting if sequence/atom not available.
    """
    pos = a.residue_index
    atom = (a.atom_name or "").strip()
    if not pos:
        return _format_res_idx(a)
    seq = seq_by_chain.get(a.chain_id, "")
    res_letter = seq[pos - 1].upper() if 1 <= pos <= len(seq) else ""
    atom_letter: str = ""
    if atom:
        atom_letter = atom[0].upper()
    else:
        if res_letter == "N":
            atom_letter = "N"
        elif res_letter in {"S", "T"}:
            atom_letter = "O"
    # If we could not determine residue letter, fall back to index-only
    left = f"{res_letter}{pos}" if res_letter else str(pos)
    # If no atom letter, omit @
    right = f"@{atom_letter}" if atom_letter else ""
    return f"{left}{right}"


def _parse_fasta(text: str) -> List[Tuple[str, str, str]]:
    entries: List[Tuple[str, str, str]] = []
    etype: Optional[str] = None
    name: Optional[str] = None
    buf: List[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if etype is not None:
                entries.append((etype, name or "", "".join(buf)))
            buf = []
            header = line[1:]
            parts = header.split("|")
            etype = parts[0]
            name = parts[1] if len(parts) > 1 else ""
        else:
            buf.append(line)
    if etype is not None:
        entries.append((etype, name or "", "".join(buf)))
    return entries


def _expand_chain_ids(n: int) -> List[str]:
    # Deprecated: use spreadsheet_ids instead
    return spreadsheet_ids(n)


def _parse_restraints_csv(text: str, valid_chain_ids: List[str]) -> List[Bond]:
    bonds: List[Bond] = []
    reader = csv.DictReader(text.splitlines())
    for row in reader:
        if row.get("connection_type", "").lower() != "covalent":
            continue
        chainA = row.get("chainA", "").strip()
        chainB = row.get("chainB", "").strip()
        if not chainA or not chainB:
            continue
        a = _parse_res_token(row.get("res_idxA", ""))
        b = _parse_res_token(row.get("res_idxB", ""))
        bonds.append(
            Bond(
                atom1=AtomRef(chain_id=chainA, residue_index=a[0], atom_name=a[1]),
                atom2=AtomRef(chain_id=chainB, residue_index=b[0], atom_name=b[1]),
            )
        )
    return bonds


def _parse_res_token(token: str) -> Tuple[Optional[int], Optional[str]]:
    token = token.strip()
    if not token:
        return None, None
    if "@" in token:
        # e.g., N436@N or @C1
        left, right = token.split("@", 1)
        if left and left[0].isalpha():
            # strip residue letter prefix (e.g., N436 -> 436)
            digits = "".join([c for c in left if c.isdigit()])
            idx = int(digits) if digits else None
        else:
            idx = int(left) if left.isdigit() else None
        return idx, right or None
    # no atom name
    digits = "".join([c for c in token if c.isdigit()])
    return (int(digits) if digits else None), None

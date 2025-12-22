from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

from foldir.io.csv import read_csv, write_csv
from foldir.io.fasta import read_fasta, write_fasta
from foldir.ir.glycans.parse_chai import parse_chai_glycan_string
from foldir.ir.glycans.render_chai import render_chai_glycan_string
from foldir.ir.models import (
    AtomRef,
    ContactConstraint,
    CovalentBond,
    DocumentIR,
    Glycan,
    GlycanAttachment,
    JobIR,
    Ligand,
    LigandReprType,
    PocketConstraint,
    PolymerChain,
    PolymerType,
)


class Chai1Dialect:
    name = "chai1"

    def parse(self, path: str) -> DocumentIR:
        fasta_path, restraints_path = _resolve_inputs(path)
        records = read_fasta(fasta_path)
        job, glycans, glycan_by_chain = _parse_fasta(records)
        if restraints_path and os.path.exists(restraints_path):
            bonds, constraints, glycan_attachments = _parse_restraints(
                restraints_path, glycan_by_chain, glycans
            )
            job = job.model_copy(
                update={
                    "covalent_bonds": bonds,
                    "constraints": constraints,
                    "glycans": glycan_attachments,
                }
            )
        else:
            job = job.model_copy(update={"glycans": glycans})
        return DocumentIR(jobs=[job])

    def render(self, doc: DocumentIR, out_path: str) -> None:
        if len(doc.jobs) != 1:
            raise ValueError("Chai-1 expects exactly one job per directory.")
        job = doc.jobs[0]
        out_dir, fasta_path, restraints_path = _resolve_outputs(out_path)
        os.makedirs(out_dir, exist_ok=True)

        records, chain_letters = _render_fasta(job)
        write_fasta(fasta_path, records)

        rows = _render_restraints(job, chain_letters)
        if rows:
            fieldnames = [
                "restraint_id",
                "chainA",
                "res_idxA",
                "chainB",
                "res_idxB",
                "connection_type",
                "confidence",
                "min_distance_angstrom",
                "max_distance_angstrom",
                "comment",
            ]
            write_csv(restraints_path, rows, fieldnames)


def _resolve_inputs(path: str) -> Tuple[str, Optional[str]]:
    if os.path.isdir(path):
        fasta_path = _find_first(path, (".fasta", ".fa"))
        restraints_path = _find_first(path, (".csv", ".restraints"))
        if not fasta_path:
            raise ValueError("Chai input directory must contain a FASTA file.")
        return fasta_path, restraints_path
    if path.endswith((".fasta", ".fa")):
        return path, None
    raise ValueError("Chai input must be a directory or FASTA path.")


def _resolve_outputs(out_path: str) -> Tuple[str, str, str]:
    if out_path.endswith((".fasta", ".fa")):
        out_dir = os.path.dirname(out_path) or "."
        fasta_path = out_path
        restraints_path = os.path.join(out_dir, "restraints.csv")
        return out_dir, fasta_path, restraints_path
    out_dir = out_path
    fasta_path = os.path.join(out_dir, "chains.fasta")
    restraints_path = os.path.join(out_dir, "restraints.csv")
    return out_dir, fasta_path, restraints_path


def _find_first(root: str, exts: Tuple[str, ...]) -> Optional[str]:
    for fname in os.listdir(root):
        if fname.endswith(exts):
            return os.path.join(root, fname)
    return None


def _parse_fasta(
    records: List[Tuple[str, str]]
) -> Tuple[JobIR, List[Glycan], Dict[str, Glycan]]:
    polymers: List[PolymerChain] = []
    ligands: List[Ligand] = []
    glycans: List[Glycan] = []
    glycan_by_chain: Dict[str, Glycan] = {}

    for idx, (header, seq) in enumerate(records, start=1):
        chain_id = _chain_id(idx)
        parts = header.split("|", 1)
        kind = parts[0].strip().lower()
        if kind == "protein":
            polymers.append(PolymerChain(id=chain_id, type=PolymerType.protein, sequence=seq))
        elif kind == "dna":
            polymers.append(PolymerChain(id=chain_id, type=PolymerType.dna, sequence=seq))
        elif kind == "rna":
            polymers.append(PolymerChain(id=chain_id, type=PolymerType.rna, sequence=seq))
        elif kind == "glycan":
            glycan_id = f"{chain_id}_glycan"
            glycan = parse_chai_glycan_string(glycan_id, seq)
            glycans.append(glycan)
            glycan_by_chain[chain_id] = glycan
        elif kind == "ligand":
            ligands.append(_parse_ligand(chain_id, seq))

    job = JobIR(name="job", polymers=polymers, ligands=ligands, glycans=glycans)
    return job, glycans, glycan_by_chain


def _parse_ligand(chain_id: str, seq: str) -> Ligand:
    if seq.isalnum() and seq.upper() == seq and len(seq) <= 6:
        return Ligand(id=chain_id, repr_type=LigandReprType.ccd, ccd_codes=[seq])
    return Ligand(id=chain_id, repr_type=LigandReprType.smiles, smiles=seq)


def _parse_restraints(
    path: str,
    glycan_by_chain: Dict[str, Glycan],
    glycans: List[Glycan],
) -> Tuple[List[CovalentBond], List[object], List[Glycan]]:
    rows = read_csv(path)
    bonds: List[CovalentBond] = []
    constraints: List[object] = []
    attachments_by_glycan: Dict[str, List[GlycanAttachment]] = {}

    for row in rows:
        connection_type = (row.get("connection_type") or "").strip().lower()
        chain_a = row.get("chainA", "").strip()
        chain_b = row.get("chainB", "").strip()
        res_a = row.get("res_idxA", "").strip()
        res_b = row.get("res_idxB", "").strip()

        if connection_type == "covalent":
            if chain_a in glycan_by_chain or chain_b in glycan_by_chain:
                g_chain = chain_a if chain_a in glycan_by_chain else chain_b
                p_chain = chain_b if g_chain == chain_a else chain_a
                glycan = glycan_by_chain[g_chain]
                p_token = res_b if g_chain == chain_a else res_a
                g_token = res_a if g_chain == chain_a else res_b
                p_pos, p_atom = _parse_residue_token(p_token)
                _, g_atom = _parse_residue_token(g_token)
                if p_pos is None:
                    continue
                attachment = GlycanAttachment(
                    polymer_id=p_chain,
                    polymer_residue_index=p_pos,
                    polymer_atom=p_atom,
                    root_node=f"{glycan.glycan_id}.n0",
                    root_atom=g_atom or "C1",
                )
                attachments_by_glycan.setdefault(glycan.glycan_id, []).append(attachment)
            else:
                a_ref = _atomref_from_token(chain_a, res_a)
                b_ref = _atomref_from_token(chain_b, res_b)
                bonds.append(CovalentBond(a=a_ref, b=b_ref))
        elif connection_type == "contact":
            constraints.append(
                ContactConstraint(
                    token1=_atomref_from_token(chain_a, res_a),
                    token2=_atomref_from_token(chain_b, res_b),
                    max_distance_angstrom=float(row.get("max_distance_angstrom", 0) or 0),
                )
            )
        elif connection_type == "pocket":
            if res_a:
                contacts = [_atomref_from_token(chain_a, res_a)]
                binder = chain_b
            else:
                contacts = [_atomref_from_token(chain_b, res_b)]
                binder = chain_a
            constraints.append(
                PocketConstraint(
                    binder_entity_id=binder,
                    contacts=contacts,
                    max_distance_angstrom=float(row.get("max_distance_angstrom", 0) or 0),
                )
            )

    if attachments_by_glycan:
        new_glycans: List[Glycan] = []
        for g in glycans:
            if g.glycan_id in attachments_by_glycan:
                new_glycans.append(
                    g.model_copy(update={"attachments": attachments_by_glycan[g.glycan_id]})
                )
            else:
                new_glycans.append(g)
        glycans = new_glycans

    return bonds, constraints, glycans


def _render_fasta(job: JobIR) -> Tuple[List[Tuple[str, str]], Dict[str, str]]:
    records: List[Tuple[str, str]] = []
    chain_letters: Dict[str, str] = {}
    chain_idx = 0

    for p in job.polymers:
        chain_idx += 1
        letter = _chain_id(chain_idx)
        chain_letters[p.id] = letter
        records.append((f"{p.type.value}|{p.id}", p.sequence))

    for lig in job.ligands:
        chain_idx += 1
        letter = _chain_id(chain_idx)
        chain_letters[lig.id] = letter
        records.append(("ligand|%s" % lig.id, _ligand_sequence(lig)))

    for g in job.glycans:
        chain_idx += 1
        letter = _chain_id(chain_idx)
        chain_letters[g.glycan_id] = letter
        root_node = g.attachments[0].root_node if g.attachments else g.nodes[0].node_id
        glycan_string = render_chai_glycan_string(g, root_node)
        records.append(("glycan|%s" % g.glycan_id, glycan_string))

    return records, chain_letters


def _render_restraints(job: JobIR, chain_letters: Dict[str, str]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    idx = 1

    for g in job.glycans:
        glycan_chain = chain_letters.get(g.glycan_id)
        if not glycan_chain:
            continue
        for att in g.attachments:
            polymer_chain = chain_letters.get(att.polymer_id)
            if not polymer_chain:
                continue
            res_letter = _residue_letter(job, att.polymer_id, att.polymer_residue_index)
            atom = att.polymer_atom or "N"
            rows.append(
                {
                    "restraint_id": f"bond{idx}",
                    "chainA": polymer_chain,
                    "res_idxA": f"{res_letter}{att.polymer_residue_index}@{atom}",
                    "chainB": glycan_chain,
                    "res_idxB": f"@{att.root_atom or 'C1'}",
                    "connection_type": "covalent",
                    "confidence": "1.0",
                    "min_distance_angstrom": "0.0",
                    "max_distance_angstrom": "0.0",
                    "comment": "protein-glycan",
                }
            )
            idx += 1

    for bond in job.covalent_bonds:
        chain_a = chain_letters.get(bond.a.entity_id)
        chain_b = chain_letters.get(bond.b.entity_id)
        if not chain_a or not chain_b:
            continue
        rows.append(
            {
                "restraint_id": f"bond{idx}",
                "chainA": chain_a,
                "res_idxA": _format_res_idx(job, bond.a),
                "chainB": chain_b,
                "res_idxB": _format_res_idx(job, bond.b),
                "connection_type": "covalent",
                "confidence": "1.0",
                "min_distance_angstrom": "0.0",
                "max_distance_angstrom": "0.0",
                "comment": "covalent",
            }
        )
        idx += 1

    for constraint in job.constraints:
        if isinstance(constraint, ContactConstraint):
            rows.append(
                {
                    "restraint_id": f"restraint{idx}",
                    "chainA": chain_letters.get(constraint.token1.entity_id, ""),
                    "res_idxA": _format_res_idx(job, constraint.token1),
                    "chainB": chain_letters.get(constraint.token2.entity_id, ""),
                    "res_idxB": _format_res_idx(job, constraint.token2),
                    "connection_type": "contact",
                    "confidence": "1.0",
                    "min_distance_angstrom": "0.0",
                    "max_distance_angstrom": str(constraint.max_distance_angstrom),
                    "comment": "contact",
                }
            )
            idx += 1
        elif isinstance(constraint, PocketConstraint):
            for contact in constraint.contacts:
                rows.append(
                    {
                        "restraint_id": f"restraint{idx}",
                        "chainA": chain_letters.get(constraint.binder_entity_id, ""),
                        "res_idxA": "",
                        "chainB": chain_letters.get(contact.entity_id, ""),
                        "res_idxB": _format_res_idx(job, contact),
                        "connection_type": "pocket",
                        "confidence": "1.0",
                        "min_distance_angstrom": "0.0",
                        "max_distance_angstrom": str(constraint.max_distance_angstrom),
                        "comment": "pocket",
                    }
                )
                idx += 1

    return rows


def _format_res_idx(job: JobIR, atom: AtomRef) -> str:
    if not any(p.id == atom.entity_id for p in job.polymers):
        return f"@{atom.atom}"
    res_letter = _residue_letter(job, atom.entity_id, atom.position)
    if isinstance(atom.atom, str):
        return f"{res_letter}{atom.position}@{atom.atom}"
    return f"{res_letter}{atom.position}@{atom.atom}"


def _residue_letter(job: JobIR, entity_id: str, position: int) -> str:
    for p in job.polymers:
        if p.id == entity_id and 0 < position <= len(p.sequence):
            return p.sequence[position - 1]
    return "X"


def _parse_residue_token(token: str) -> Tuple[Optional[int], Optional[str]]:
    if not token:
        return None, None
    if "@" in token:
        res_part, atom = token.split("@", 1)
    else:
        res_part, atom = token, None
    if not res_part:
        return None, atom
    digits = "".join(ch for ch in res_part if ch.isdigit())
    if not digits:
        return None, atom
    return int(digits), atom


def _atomref_from_token(chain_id: str, token: str) -> AtomRef:
    pos, atom = _parse_residue_token(token)
    if pos is None:
        pos = 1
    return AtomRef(entity_id=chain_id, position=pos, atom=atom or "CA")


def _ligand_sequence(lig: Ligand) -> str:
    if lig.repr_type.value == "ccd":
        return lig.ccd_codes[0]
    if lig.repr_type.value == "smiles":
        return lig.smiles or ""
    return lig.file_path or ""


def _chain_id(idx: int) -> str:
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    out = ""
    n = idx
    while n > 0:
        n -= 1
        out = alphabet[n % 26] + out
        n //= 26
    return out

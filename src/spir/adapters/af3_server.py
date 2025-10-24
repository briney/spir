from __future__ import annotations

import json
from typing import Any, Dict, List, Optional, Tuple

from ..ir.glycan import parse_af3_server_glycan
from ..ir.model import AtomRef, Bond, ComplexInput, Ion, Ligand, PolymerChain
from ..utils.ids import spreadsheet_ids


class AF3ServerReaderWriter:
    @staticmethod
    def load_str(data: str) -> ComplexInput:
        arr = json.loads(data)
        if not isinstance(arr, list):
            raise ValueError("AF3 Server JSON expects a top-level list of jobs")
        # Support first job only for now
        obj = arr[0]
        name = obj.get("name")
        sequences = obj.get("sequences", [])
        proteins: List[PolymerChain] = []
        rnas: List[PolymerChain] = []
        dnas: List[PolymerChain] = []
        ligands: List[Ligand] = []
        ions: List[Ion] = []
        bonds: List[Bond] = []

        # Assign chain ids in order of appearance: A, B, ...
        chain_alloc_index: int = 0

        def new_ids(count: int) -> List[str]:
            nonlocal chain_alloc_index
            ids = spreadsheet_ids(chain_alloc_index + count)[chain_alloc_index:]
            chain_alloc_index += count
            return ids

        for entry in sequences:
            if "proteinChain" in entry:
                p = entry["proteinChain"]
                ids = new_ids(int(p.get("count", 1)))
                proteins.append(
                    PolymerChain(type="protein", ids=ids, sequence=p["sequence"])
                )
                # Glycans field (UI) is a simplified spec; we record as separate glycan ligands if present
                for g in p.get("glycans", []) or []:
                    # residues: like NAG(NAG)(BMA); position: anchor residue index
                    comp_spec = g.get("residues", "")
                    graph = parse_af3_server_glycan(comp_spec)
                    gid = new_ids(1)
                    # Multi-CCD ligand with original residues string preserved
                    ligands.append(
                        Ligand(
                            ids=gid,
                            ccd_codes=graph.components_as_ccd_list(),
                            af3_residues=comp_spec,
                        )
                    )
                    # Intra-glycan bonds: parent O[n] -> child C1 addressed by component indices
                    for b in graph.bonds:
                        bonds.append(
                            Bond(
                                atom1=AtomRef(
                                    chain_id=gid[0],
                                    residue_index=b.parent_index,
                                    component_index=b.parent_index,
                                    atom_name=b.parent_atom,
                                ),
                                atom2=AtomRef(
                                    chain_id=gid[0],
                                    residue_index=b.child_index,
                                    component_index=b.child_index,
                                    atom_name=b.child_atom,
                                ),
                            )
                        )
                    # Anchor bond from protein residue to glycan root C1
                    pos = int(g.get("position", 1))
                    anchor_atom = _infer_anchor_atom_name(p.get("sequence", ""), pos)
                    bonds.append(
                        Bond(
                            atom1=AtomRef(
                                chain_id=ids[0],
                                residue_index=pos,
                                atom_name=anchor_atom,
                            ),
                            atom2=AtomRef(
                                chain_id=gid[0],
                                residue_index=1,
                                component_index=1,
                                atom_name="C1",
                            ),
                        )
                    )
            elif "dnaSequence" in entry:
                p = entry["dnaSequence"]
                ids = new_ids(int(p.get("count", 1)))
                dnas.append(PolymerChain(type="dna", ids=ids, sequence=p["sequence"]))
            elif "rnaSequence" in entry:
                p = entry["rnaSequence"]
                ids = new_ids(int(p.get("count", 1)))
                rnas.append(PolymerChain(type="rna", ids=ids, sequence=p["sequence"]))
            elif "ligand" in entry:
                p = entry["ligand"]
                ids = new_ids(int(p.get("count", 1)))
                # AF3 Server restricts ligands to a curated CCD set; map directly
                ligands.append(
                    Ligand(ids=ids, ccd_codes=[_strip_ccd_prefix(p["ligand"])])
                )
            elif "ion" in entry:
                p = entry["ion"]
                ids = new_ids(int(p.get("count", 1)))
                ions.append(Ion(ids=ids, code=p["ion"]))
            else:
                raise ValueError("Unknown sequence entry in AF3 Server JSON")

        seeds = [int(s) for s in obj.get("modelSeeds", [])]
        return ComplexInput(
            name=name,
            seeds=seeds,
            proteins=proteins,
            rnas=rnas,
            dnas=dnas,
            ligands=ligands,
            ions=ions,
            bonds=bonds,
        )

    @staticmethod
    def dump(ci: ComplexInput) -> str:
        sequences: List[Dict[str, Any]] = []

        # Build quick lookups
        id_to_protein_index: Dict[str, int] = {}
        for idx, p in enumerate(ci.proteins):
            for cid in p.ids:
                id_to_protein_index[cid] = idx
        id_to_ligand: Dict[str, Ligand] = {}
        for ligand in ci.ligands:
            for cid in ligand.ids:
                id_to_ligand[cid] = ligand

        # Collect protein-attached glycans: map protein index -> list of (position, ligand_id)
        protein_glycans: Dict[int, List[Tuple[int, str]]] = {
            i: [] for i in range(len(ci.proteins))
        }
        seen_pairs: Dict[int, set] = {i: set() for i in range(len(ci.proteins))}
        anchored_ligand_ids: set = set()
        for b in ci.bonds:
            p_idx: Optional[int] = None
            pos: Optional[int] = None
            lig_id: Optional[str] = None
            # Case: atom1 is protein, atom2 is ligand
            if (
                b.atom1.chain_id in id_to_protein_index
                and b.atom2.chain_id in id_to_ligand
            ):
                p_idx = id_to_protein_index[b.atom1.chain_id]
                pos = int(b.atom1.residue_index or 1)
                lig_id = b.atom2.chain_id
            # Case: atom2 is protein, atom1 is ligand
            elif (
                b.atom2.chain_id in id_to_protein_index
                and b.atom1.chain_id in id_to_ligand
            ):
                p_idx = id_to_protein_index[b.atom2.chain_id]
                pos = int(b.atom2.residue_index or 1)
                lig_id = b.atom1.chain_id
            if p_idx is None or pos is None or lig_id is None:
                continue
            # Only consider CCD-based ligands
            lig = id_to_ligand.get(lig_id)
            if lig is None or lig.ccd_codes is None or len(lig.ccd_codes) == 0:
                continue
            key = (pos, lig_id)
            if key not in seen_pairs[p_idx]:
                protein_glycans[p_idx].append((pos, lig_id))
                seen_pairs[p_idx].add(key)
                anchored_ligand_ids.add(lig_id)

        # Proteins with optional glycans
        for idx, p in enumerate(ci.proteins):
            entry: Dict[str, Any] = {
                "proteinChain": {"sequence": p.sequence, "count": len(p.ids)}
            }
            glycans_out: List[Dict[str, Any]] = []
            for pos, lig_id in protein_glycans.get(idx, []) or []:
                lig = id_to_ligand[lig_id]
                residues = _render_af3_glycan(lig)
                glycans_out.append({"residues": residues, "position": pos})
            if glycans_out:
                entry["proteinChain"]["glycans"] = glycans_out
            sequences.append(entry)

        # DNA chains
        for d in ci.dnas:
            sequences.append(
                {"dnaSequence": {"sequence": d.sequence, "count": len(d.ids)}}
            )

        # RNA chains
        for r in ci.rnas:
            sequences.append(
                {"rnaSequence": {"sequence": r.sequence, "count": len(r.ids)}}
            )

        # Top-level ligands (single CCD only), excluding anchored ones
        unsupported: List[str] = []
        for ligand in ci.ligands:
            # Skip if any of its ids are anchored
            if any((cid in anchored_ligand_ids) for cid in ligand.ids):
                continue
            if ligand.ccd_codes is not None:
                if len(ligand.ccd_codes) == 1:
                    sequences.append(
                        {
                            "ligand": {
                                "ligand": f"CCD_{ligand.ccd_codes[0]}",
                                "count": len(ligand.ids),
                            }
                        }
                    )
                else:
                    unsupported.append("multi-CCD unanchored glycan")
            else:
                # SMILES or other ligand cannot be expressed in AF3 Server inputs
                unsupported.append("SMILES ligand not supported in AF3 Server")

        # Ions
        for i in ci.ions:
            sequences.append({"ion": {"ion": i.code, "count": len(i.ids)}})

        if unsupported:
            raise ValueError(
                "Cannot write AF3 Server for unsupported entities: "
                + ", ".join(sorted(set(unsupported)))
            )

        job: Dict[str, Any] = {
            "name": ci.name or "spir-job",
            "modelSeeds": [str(s) for s in (ci.seeds or [])],
            "sequences": sequences,
            "dialect": "alphafoldserver",
            "version": 1,
        }
        return json.dumps([job], indent=2)


def _strip_ccd_prefix(s: str) -> str:
    return s[4:] if s.startswith("CCD_") else s


def _infer_anchor_atom_name(sequence: str, position_1based: int) -> Optional[str]:
    # Conservative inference for common glycosylation link atoms
    try:
        if position_1based >= 1 and position_1based <= len(sequence):
            aa = sequence[position_1based - 1].upper()
            if aa == "N":
                return "ND2"
            if aa == "S":
                return "OG"
            if aa == "T":
                return "OG1"
    except Exception:
        pass
    return None


def _best_effort_components(glycan_residues_spec: str) -> List[str]:
    # Very light extraction: take tokens that look like CCD codes (alnum/_/-) ignoring parentheses and links
    comp: List[str] = []
    buf = []
    i = 0
    s = glycan_residues_spec.strip()
    while i < len(s):
        ch = s[i]
        if ch.isalnum() or ch in {"_", "-"}:
            buf.append(ch)
            i += 1
            continue
        if buf:
            comp.append("".join(buf))
            buf = []
        i += 1
    if buf:
        comp.append("".join(buf))
    return [c for c in comp if c]


def _render_af3_glycan(lig: Ligand) -> str:
    if lig.af3_residues:
        return lig.af3_residues
    if lig.ccd_codes is None or len(lig.ccd_codes) == 0:
        raise ValueError(
            "Ligand must have ccd_codes or af3_residues for AF3 Server glycan"
        )
    # Left-nested rendering: [A,B,C] -> A(B(C))
    codes = list(lig.ccd_codes)
    s = codes[-1]
    for c in reversed(codes[:-1]):
        s = f"{c}({s})"
    return s


# Backward-compat alias
AF3ServerReader = AF3ServerReaderWriter

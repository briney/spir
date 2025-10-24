from __future__ import annotations

import json
from typing import Any, Dict, List, Tuple, Union

from ..ir.model import AtomRef, Bond, ComplexInput, Ion, Ligand, PolymerChain
from ..utils.ids import spreadsheet_ids


def _split_ccd_concat(code: str) -> List[str]:
    # Expect strings like CCD_NAG_BMA_BGC
    if not code.startswith("CCD_"):
        return [code]
    return [p for p in code.split("_")[1:] if p]


class ProtenixReaderWriter:
    @staticmethod
    def load_str(data: str) -> ComplexInput:
        arr = json.loads(data)
        if not isinstance(arr, list):
            raise ValueError("Protenix JSON expects a top-level list")
        if len(arr) != 1:
            # For simplicity, only support a single job per file for now
            obj = arr[0]
        else:
            obj = arr[0]
        name = obj.get("name")
        seqs = obj.get("sequences", [])
        proteins: List[PolymerChain] = []
        rnas: List[PolymerChain] = []
        dnas: List[PolymerChain] = []
        ligands: List[Ligand] = []
        ions: List[Ion] = []
        bonds: List[Bond] = []

        # Track implicit entity indices by order
        chain_ids_by_entity_index: List[List[str]] = []

        for entry in seqs:
            if "proteinChain" in entry:
                p = entry["proteinChain"]
                count = int(p.get("count", 1))
                ids = spreadsheet_ids(count)
                proteins.append(
                    PolymerChain(type="protein", ids=ids, sequence=p["sequence"])
                )
                chain_ids_by_entity_index.append(ids)
            elif "dnaSequence" in entry:
                p = entry["dnaSequence"]
                count = int(p.get("count", 1))
                ids = spreadsheet_ids(count)
                dnas.append(PolymerChain(type="dna", ids=ids, sequence=p["sequence"]))
                chain_ids_by_entity_index.append(ids)
            elif "rnaSequence" in entry:
                p = entry["rnaSequence"]
                count = int(p.get("count", 1))
                ids = spreadsheet_ids(count)
                rnas.append(PolymerChain(type="rna", ids=ids, sequence=p["sequence"]))
                chain_ids_by_entity_index.append(ids)
            elif "ligand" in entry:
                p = entry["ligand"]
                count = int(p.get("count", 1))
                ids = spreadsheet_ids(count)
                lig = p.get("ligand")
                if isinstance(lig, str) and lig.startswith("CCD_"):
                    ccds = _split_ccd_concat(lig)
                    ligands.append(Ligand(ids=ids, ccd_codes=ccds))
                else:
                    # Treat as SMILES or FILE_
                    ligands.append(Ligand(ids=ids, smiles=str(lig)))
                chain_ids_by_entity_index.append(ids)
            elif "ion" in entry:
                p = entry["ion"]
                count = int(p.get("count", 1))
                ids = spreadsheet_ids(count)
                ions.append(Ion(ids=ids, code=p["ion"]))
                chain_ids_by_entity_index.append(ids)
            else:
                raise ValueError("Unknown sequence entry in Protenix JSON")

        for cb in obj.get("covalent_bonds", []) or []:
            e1 = int(cb["entity1"]) - 1
            e2 = int(cb["entity2"]) - 1
            c1 = int(cb.get("copy1", 1))
            c2 = int(cb.get("copy2", 1))
            pos1 = int(cb["position1"]) if "position1" in cb else 1
            pos2 = int(cb["position2"]) if "position2" in cb else 1
            atom1 = cb.get("atom1")
            atom2 = cb.get("atom2")
            chain1 = chain_ids_by_entity_index[e1][c1 - 1]
            chain2 = chain_ids_by_entity_index[e2][c2 - 1]
            bonds.append(
                Bond(
                    atom1=AtomRef(
                        chain_id=chain1,
                        residue_index=pos1,
                        component_index=pos1,
                        atom_name=atom1 if atom1 and not _is_int(atom1) else None,
                        atom_index=int(atom1) if _is_int(atom1) else None,
                    ),
                    atom2=AtomRef(
                        chain_id=chain2,
                        residue_index=pos2,
                        component_index=pos2,
                        atom_name=atom2 if atom2 and not _is_int(atom2) else None,
                        atom_index=int(atom2) if _is_int(atom2) else None,
                    ),
                )
            )

        return ComplexInput(
            name=name,
            proteins=proteins,
            rnas=rnas,
            dnas=dnas,
            ligands=ligands,
            ions=ions,
            bonds=bonds,
        )

    @staticmethod
    def dump(ci: ComplexInput) -> str:
        seqs: List[Dict[str, Any]] = []
        # Entity order matches sequences appended
        for p in ci.proteins:
            seqs.append({"proteinChain": {"sequence": p.sequence, "count": len(p.ids)}})
        for p in ci.dnas:
            seqs.append({"dnaSequence": {"sequence": p.sequence, "count": len(p.ids)}})
        for p in ci.rnas:
            seqs.append({"rnaSequence": {"sequence": p.sequence, "count": len(p.ids)}})
        for l in ci.ligands:
            if l.ccd_codes is not None:
                seqs.append(
                    {
                        "ligand": {
                            "ligand": "CCD_" + "_".join(l.ccd_codes),
                            "count": len(l.ids),
                        }
                    }
                )
            else:
                seqs.append({"ligand": {"ligand": l.smiles, "count": len(l.ids)}})
        for i in ci.ions:
            seqs.append({"ion": {"ion": i.code, "count": len(i.ids)}})

        # Build mapping from entity index to chain ids to support bonds
        entity_chain_ids: List[List[str]] = []
        for p in ci.proteins:
            entity_chain_ids.append(p.ids)
        for p in ci.dnas:
            entity_chain_ids.append(p.ids)
        for p in ci.rnas:
            entity_chain_ids.append(p.ids)
        for l in ci.ligands:
            entity_chain_ids.append(l.ids)
        for i in ci.ions:
            entity_chain_ids.append(i.ids)

        def find_entity_and_copy(chain_id: str) -> Tuple[int, int]:
            for e_idx, id_list in enumerate(entity_chain_ids):
                if chain_id in id_list:
                    return e_idx + 1, id_list.index(chain_id) + 1
            raise ValueError(f"Unknown chain id in bond: {chain_id}")

        covalent_bonds: List[Dict[str, Any]] = []
        for b in ci.bonds:
            e1, c1 = find_entity_and_copy(b.atom1.chain_id)
            e2, c2 = find_entity_and_copy(b.atom2.chain_id)
            covalent_bonds.append(
                {
                    "entity1": str(e1),
                    "copy1": c1,
                    "position1": str(
                        b.atom1.residue_index or b.atom1.component_index or 1
                    ),
                    "atom1": b.atom1.atom_name
                    if b.atom1.atom_name is not None
                    else (
                        b.atom1.atom_index if b.atom1.atom_index is not None else None
                    ),
                    "entity2": str(e2),
                    "copy2": c2,
                    "position2": str(
                        b.atom2.residue_index or b.atom2.component_index or 1
                    ),
                    "atom2": b.atom2.atom_name
                    if b.atom2.atom_name is not None
                    else (
                        b.atom2.atom_index if b.atom2.atom_index is not None else None
                    ),
                }
            )

        out = [
            {
                "name": ci.name or "spir-job",
                "sequences": seqs,
                "covalent_bonds": covalent_bonds if covalent_bonds else None,
            }
        ]
        if out[0]["covalent_bonds"] is None:
            del out[0]["covalent_bonds"]
        return json.dumps(out, indent=2)


def _is_int(v: Union[str, int, None]) -> bool:
    if isinstance(v, int):
        return True
    if isinstance(v, str) and v.isdigit():
        return True
    return False

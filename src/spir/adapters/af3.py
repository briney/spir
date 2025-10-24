from __future__ import annotations

import json
from typing import Any, Dict, List, Union

from pydantic import BaseModel

from ..ir.model import AtomRef, Bond, ComplexInput, Ion, Ligand, PolymerChain


class AF3ReaderWriter:
    @staticmethod
    def load_str(data: str) -> ComplexInput:
        obj = json.loads(data)
        if not isinstance(obj, dict) or obj.get("dialect") != "alphafold3":
            raise ValueError("Not an AlphaFold3 JSON (dialect=alphafold3)")
        sequences = obj.get("sequences", [])
        proteins: List[PolymerChain] = []
        rnas: List[PolymerChain] = []
        dnas: List[PolymerChain] = []
        ligands: List[Ligand] = []
        ions: List[Ion] = []
        bonds: List[Bond] = []
        for entry in sequences:
            if "protein" in entry:
                p = entry["protein"]
                ids = p["id"] if isinstance(p["id"], list) else [p["id"]]
                proteins.append(
                    PolymerChain(
                        type="protein",
                        ids=ids,
                        sequence=p["sequence"],
                    )
                )
            elif "rna" in entry:
                p = entry["rna"]
                ids = p["id"] if isinstance(p["id"], list) else [p["id"]]
                rnas.append(PolymerChain(type="rna", ids=ids, sequence=p["sequence"]))
            elif "dna" in entry:
                p = entry["dna"]
                ids = p["id"] if isinstance(p["id"], list) else [p["id"]]
                dnas.append(PolymerChain(type="dna", ids=ids, sequence=p["sequence"]))
            elif "ligand" in entry:
                p = entry["ligand"]
                ids = p["id"] if isinstance(p["id"], list) else [p["id"]]
                if "ccdCodes" in p:
                    ligands.append(Ligand(ids=ids, ccd_codes=list(p["ccdCodes"])))
                elif "smiles" in p:
                    ligands.append(Ligand(ids=ids, smiles=p["smiles"]))
                else:
                    raise ValueError("AF3 ligand requires ccdCodes or smiles")
            else:
                raise ValueError("Unknown sequence entry in AF3 JSON")

        for b in obj.get("bondedAtomPairs", []) or []:
            (e1, r1, a1), (e2, r2, a2) = b
            bonds.append(
                Bond(
                    atom1=AtomRef(
                        chain_id=str(e1), residue_index=int(r1), atom_name=str(a1)
                    ),
                    atom2=AtomRef(
                        chain_id=str(e2), residue_index=int(r2), atom_name=str(a2)
                    ),
                )
            )

        seeds = obj.get("modelSeeds", [])
        name = obj.get("name")
        user_ccd = obj.get("userCCD") or obj.get("userCCDPath")
        return ComplexInput(
            name=name,
            seeds=[int(s) for s in seeds],
            proteins=proteins,
            rnas=rnas,
            dnas=dnas,
            ligands=ligands,
            ions=ions,
            bonds=bonds,
            user_ccd=user_ccd,
        )

    @staticmethod
    def dump(ci: ComplexInput, version: int = 4) -> str:
        sequences: List[Dict[str, Any]] = []
        for p in ci.proteins:
            ident: Union[str, List[str]] = p.ids if len(p.ids) > 1 else p.ids[0]
            sequences.append({"protein": {"id": ident, "sequence": p.sequence}})
        for p in ci.rnas:
            ident = p.ids if len(p.ids) > 1 else p.ids[0]
            sequences.append({"rna": {"id": ident, "sequence": p.sequence}})
        for p in ci.dnas:
            ident = p.ids if len(p.ids) > 1 else p.ids[0]
            sequences.append({"dna": {"id": ident, "sequence": p.sequence}})
        for l in ci.ligands:
            ident = l.ids if len(l.ids) > 1 else l.ids[0]
            if l.ccd_codes is not None:
                sequences.append({"ligand": {"id": ident, "ccdCodes": l.ccd_codes}})
            else:
                sequences.append({"ligand": {"id": ident, "smiles": l.smiles}})

        bonded: List[List[List[Union[str, int]]]] = []
        for b in ci.bonds:
            a1 = [b.atom1.chain_id, b.atom1.residue_index or 1, b.atom1.atom_name or ""]
            a2 = [b.atom2.chain_id, b.atom2.residue_index or 1, b.atom2.atom_name or ""]
            bonded.append([a1, a2])

        out: Dict[str, Any] = {
            "name": ci.name or "spir-job",
            "modelSeeds": ci.seeds or [0],
            "sequences": sequences,
            "bondedAtomPairs": bonded if bonded else None,
            "dialect": "alphafold3",
            "version": version,
        }
        # Remove None fields
        if out["bondedAtomPairs"] is None:
            del out["bondedAtomPairs"]
        if ci.user_ccd:
            # heuristic: assume string path or mmCIF text
            key = "userCCDPath" if "\n" not in ci.user_ccd else "userCCD"
            out[key] = ci.user_ccd
        return json.dumps(out, indent=2)

from __future__ import annotations

from typing import Any, Dict, List, Union

import yaml

from ..ir.model import AtomRef, Bond, ComplexInput, Ion, Ligand, PolymerChain


class BoltzReaderWriter:
    @staticmethod
    def load_str(data: str) -> ComplexInput:
        obj = yaml.safe_load(data)
        if not isinstance(obj, dict) or "sequences" not in obj:
            raise ValueError("Not a Boltz YAML input")
        proteins: List[PolymerChain] = []
        rnas: List[PolymerChain] = []
        dnas: List[PolymerChain] = []
        ligands: List[Ligand] = []
        ions: List[Ion] = []
        bonds: List[Bond] = []

        for entry in obj.get("sequences", []) or []:
            if "protein" in entry:
                p = entry["protein"]
                ids = p.get("id")
                ids = ids if isinstance(ids, list) else [ids]
                proteins.append(
                    PolymerChain(type="protein", ids=ids, sequence=p["sequence"])
                )
            elif "rna" in entry:
                p = entry["rna"]
                ids = p.get("id")
                ids = ids if isinstance(ids, list) else [ids]
                rnas.append(PolymerChain(type="rna", ids=ids, sequence=p["sequence"]))
            elif "dna" in entry:
                p = entry["dna"]
                ids = p.get("id")
                ids = ids if isinstance(ids, list) else [ids]
                dnas.append(PolymerChain(type="dna", ids=ids, sequence=p["sequence"]))
            elif "ligand" in entry:
                p = entry["ligand"]
                ids = p.get("id")
                ids = ids if isinstance(ids, list) else [ids]
                if "ccd" in p:
                    ligands.append(Ligand(ids=ids, ccd_codes=[p["ccd"]]))
                elif "smiles" in p:
                    ligands.append(Ligand(ids=ids, smiles=p["smiles"]))
                else:
                    raise ValueError("Boltz ligand requires ccd or smiles")
            else:
                raise ValueError("Unknown sequence entry in Boltz YAML")

        for c in obj.get("constraints", []) or []:
            if "bond" in c:
                b = c["bond"]
                (cid1, pos1, atom1) = b["atom1"]
                (cid2, pos2, atom2) = b["atom2"]
                bonds.append(
                    Bond(
                        atom1=AtomRef(
                            chain_id=str(cid1),
                            residue_index=int(pos1),
                            atom_name=str(atom1),
                        ),
                        atom2=AtomRef(
                            chain_id=str(cid2),
                            residue_index=int(pos2),
                            atom_name=str(atom2),
                        ),
                    )
                )

        return ComplexInput(
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
        for p in ci.proteins:
            sequences.append(
                {
                    "protein": {
                        "id": p.ids if len(p.ids) > 1 else p.ids[0],
                        "sequence": p.sequence,
                    }
                }
            )
        for p in ci.rnas:
            sequences.append(
                {
                    "rna": {
                        "id": p.ids if len(p.ids) > 1 else p.ids[0],
                        "sequence": p.sequence,
                    }
                }
            )
        for p in ci.dnas:
            sequences.append(
                {
                    "dna": {
                        "id": p.ids if len(p.ids) > 1 else p.ids[0],
                        "sequence": p.sequence,
                    }
                }
            )
        for l in ci.ligands:
            if l.ccd_codes is not None and len(l.ccd_codes) == 1:
                sequences.append(
                    {
                        "ligand": {
                            "id": l.ids if len(l.ids) > 1 else l.ids[0],
                            "ccd": l.ccd_codes[0],
                        }
                    }
                )
            elif l.smiles is not None:
                sequences.append(
                    {
                        "ligand": {
                            "id": l.ids if len(l.ids) > 1 else l.ids[0],
                            "smiles": l.smiles,
                        }
                    }
                )
            else:
                # Multi-CCD glycans will still be modeled as ligand with single id; bonds must express connectivity
                sequences.append(
                    {
                        "ligand": {
                            "id": l.ids if len(l.ids) > 1 else l.ids[0],
                            "ccd": (l.ccd_codes or [""])[0],
                        }
                    }
                )

        constraints: List[Dict[str, Any]] = []
        for b in ci.bonds:
            if b.atom1.atom_name and b.atom2.atom_name:
                constraints.append(
                    {
                        "bond": {
                            "atom1": [
                                b.atom1.chain_id,
                                b.atom1.residue_index or 1,
                                b.atom1.atom_name,
                            ],
                            "atom2": [
                                b.atom2.chain_id,
                                b.atom2.residue_index or 1,
                                b.atom2.atom_name,
                            ],
                        }
                    }
                )

        out = {
            "version": 1,
            "sequences": sequences,
        }
        if constraints:
            out["constraints"] = constraints
        return yaml.safe_dump(out, sort_keys=False)

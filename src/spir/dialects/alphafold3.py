from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from spir.io.json import read_json, write_json
from spir.ir.models import (
    AtomRef,
    CovalentBond,
    DocumentIR,
    Glycan,
    GlycanAttachment,
    GlycanEdge,
    GlycanNode,
    JobIR,
    Ligand,
    LigandReprType,
    Modification,
    PolymerChain,
    PolymerType,
)
from spir.validate import ValidationResult


class AlphaFold3Dialect:
    name = "alphafold3"

    def parse(self, path: str) -> DocumentIR:
        payload = read_json(path)
        job = _parse_job(payload)
        glycans = _detect_glycans(job)
        if glycans:
            job = job.model_copy(update={"glycans": glycans})
        return DocumentIR(jobs=[job])

    def validate(self, path: str) -> ValidationResult:
        result = ValidationResult()
        try:
            payload = read_json(path)
        except Exception as e:
            result.add_error(f"Failed to parse JSON: {e}")
            return result

        # Check required fields
        if not isinstance(payload, dict):
            result.add_error("Input must be a JSON object")
            return result

        # Check dialect and version
        dialect = payload.get("dialect")
        if dialect is None:
            result.add_warning("Missing 'dialect' field (expected 'alphafold3')")
        elif dialect != "alphafold3":
            result.add_warning(f"Unexpected dialect '{dialect}' (expected 'alphafold3')")

        version = payload.get("version")
        if version is None:
            result.add_warning("Missing 'version' field (expected 4)")
        elif version != 4:
            result.add_warning(f"Unexpected version {version} (expected 4)")

        # Check modelSeeds
        seeds = payload.get("modelSeeds")
        if seeds is None:
            result.add_error("Missing required field 'modelSeeds'")
        elif not isinstance(seeds, list):
            result.add_error("'modelSeeds' must be a list of integers")
        elif len(seeds) == 0:
            result.add_error("'modelSeeds' must contain at least one seed")

        # Check sequences
        sequences = payload.get("sequences")
        if sequences is None:
            result.add_error("Missing required field 'sequences'")
        elif not isinstance(sequences, list):
            result.add_error("'sequences' must be a list")
        else:
            entity_ids = set()
            for idx, entry in enumerate(sequences):
                loc = f"sequences[{idx}]"
                _validate_sequence_entry(entry, loc, entity_ids, result)

            # Validate bondedAtomPairs references
            bonded_pairs = payload.get("bondedAtomPairs") or []
            for bond_idx, pair in enumerate(bonded_pairs):
                loc = f"bondedAtomPairs[{bond_idx}]"
                if not isinstance(pair, list) or len(pair) != 2:
                    result.add_error(f"Bond must be a list of two atom references", loc)
                    continue
                for atom_idx, atom in enumerate(pair):
                    if not isinstance(atom, list) or len(atom) != 3:
                        result.add_error(
                            f"Atom reference must be [entity_id, position, atom_name]",
                            f"{loc}[{atom_idx}]",
                        )
                        continue
                    entity_id = atom[0]
                    if entity_id not in entity_ids:
                        result.add_error(
                            f"Entity '{entity_id}' not found in sequences",
                            f"{loc}[{atom_idx}]",
                        )

        # Try full parse to catch additional issues
        if result.is_valid:
            try:
                self.parse(path)
            except Exception as e:
                result.add_error(f"Validation passed but parsing failed: {e}")

        return result

    def render(self, doc: DocumentIR, out_path: str) -> None:
        if len(doc.jobs) != 1:
            raise ValueError("AlphaFold3 (non-server) expects exactly one job per JSON file.")
        job = doc.jobs[0]

        seeds = job.seeds if job.seeds else [1]

        sequences = _render_sequences(job)
        bonded_pairs = _render_bonded_pairs(job)

        payload = {
            "name": job.name,
            "modelSeeds": seeds,
            "sequences": sequences,
            "bondedAtomPairs": bonded_pairs or None,
            "dialect": "alphafold3",
            "version": 4,
        }
        write_json(out_path, payload)


def _parse_job(payload: dict) -> JobIR:
    name = payload.get("name", "job")
    seeds = payload.get("modelSeeds") or []
    polymers: List[PolymerChain] = []
    ligands: List[Ligand] = []
    bonds: List[CovalentBond] = []

    for entry in payload.get("sequences", []):
        if "protein" in entry:
            p = entry["protein"]
            mods = p.get("modifications") or []
            msa_path = p.get("unpairedMsaPath")
            polymers.append(
                PolymerChain(
                    id=p["id"],
                    type=PolymerType.protein,
                    sequence=p["sequence"],
                    modifications=[
                        Modification(position=m["ptmPosition"], ccd=m["ptmType"]) for m in mods
                    ],
                    msa_path=msa_path,
                )
            )
        elif "dna" in entry:
            p = entry["dna"]
            polymers.append(
                PolymerChain(id=p["id"], type=PolymerType.dna, sequence=p["sequence"])
            )
        elif "rna" in entry:
            p = entry["rna"]
            msa_path = p.get("unpairedMsaPath")
            polymers.append(
                PolymerChain(
                    id=p["id"],
                    type=PolymerType.rna,
                    sequence=p["sequence"],
                    msa_path=msa_path,
                )
            )
        elif "ligand" in entry:
            l = entry["ligand"]
            if "ccdCodes" in l:
                ligands.append(
                    Ligand(id=l["id"], repr_type=LigandReprType.ccd, ccd_codes=l["ccdCodes"])
                )
            elif "smiles" in l:
                ligands.append(
                    Ligand(id=l["id"], repr_type=LigandReprType.smiles, smiles=l["smiles"])
                )
            elif "file" in l:
                ligands.append(
                    Ligand(id=l["id"], repr_type=LigandReprType.file, file_path=l["file"])
                )

    for pair in payload.get("bondedAtomPairs") or []:
        a = pair[0]
        b = pair[1]
        bonds.append(
            CovalentBond(
                a=AtomRef(entity_id=a[0], position=int(a[1]), atom=a[2]),
                b=AtomRef(entity_id=b[0], position=int(b[1]), atom=b[2]),
            )
        )

    return JobIR(
        name=name,
        seeds=seeds,
        polymers=polymers,
        ligands=ligands,
        covalent_bonds=bonds,
    )


def _detect_glycans(job: JobIR) -> List[Glycan]:
    polymer_ids = {p.id for p in job.polymers}
    glycans: List[Glycan] = []
    for lig in job.ligands:
        if lig.repr_type != LigandReprType.ccd or len(lig.ccd_codes) < 2:
            continue
        glycan_id = f"{lig.id}_glycan"
        nodes = [
            GlycanNode(node_id=f"{glycan_id}.n{i}", ccd=ccd)
            for i, ccd in enumerate(lig.ccd_codes)
        ]
        edges: List[GlycanEdge] = []
        attachments: List[GlycanAttachment] = []
        for bond in job.covalent_bonds:
            a = bond.a
            b = bond.b
            if a.entity_id == lig.id and b.entity_id == lig.id:
                edge = _edge_from_bond(glycan_id, a, b, len(nodes))
                if edge is not None:
                    edges.append(edge)
            elif a.entity_id in polymer_ids and b.entity_id == lig.id:
                attachments.append(_attachment_from_bond(glycan_id, a, b))
            elif b.entity_id in polymer_ids and a.entity_id == lig.id:
                attachments.append(_attachment_from_bond(glycan_id, b, a))
        if edges or attachments:
            glycans.append(
                Glycan(glycan_id=glycan_id, nodes=nodes, edges=edges, attachments=attachments)
            )
    return glycans


def _edge_from_bond(glycan_id: str, a: AtomRef, b: AtomRef, node_count: int) -> Optional[GlycanEdge]:
    if not isinstance(a.atom, str) or not isinstance(b.atom, str):
        return None
    if a.atom.startswith("O") and b.atom == "C1":
        parent_pos, child_pos = a.position, b.position
        parent_atom, child_atom = a.atom, b.atom
    elif b.atom.startswith("O") and a.atom == "C1":
        parent_pos, child_pos = b.position, a.position
        parent_atom, child_atom = b.atom, a.atom
    else:
        return None
    if parent_pos < 1 or child_pos < 1 or parent_pos > node_count or child_pos > node_count:
        return None
    return GlycanEdge(
        parent=f"{glycan_id}.n{parent_pos - 1}",
        child=f"{glycan_id}.n{child_pos - 1}",
        parent_atom=parent_atom,
        child_atom=child_atom,
    )


def _attachment_from_bond(glycan_id: str, polymer: AtomRef, ligand: AtomRef) -> GlycanAttachment:
    return GlycanAttachment(
        polymer_id=polymer.entity_id,
        polymer_residue_index=polymer.position,
        polymer_atom=polymer.atom if isinstance(polymer.atom, str) else None,
        root_node=f"{glycan_id}.n{ligand.position - 1}",
        root_atom=ligand.atom if isinstance(ligand.atom, str) else None,
    )


def _render_sequences(job: JobIR) -> List[dict]:
    sequences: List[dict] = []
    for p in job.polymers:
        if p.type.value == "protein":
            mods = [
                {"ptmType": m.ccd, "ptmPosition": m.position} for m in p.modifications
            ]
            protein_entry = {
                "id": p.id,
                "sequence": p.sequence,
                "modifications": mods or None,
            }
            if p.msa_path:
                protein_entry["unpairedMsaPath"] = p.msa_path
                protein_entry["pairedMsa"] = ""
            sequences.append({"protein": protein_entry})
        elif p.type.value == "dna":
            sequences.append({"dna": {"id": p.id, "sequence": p.sequence}})
        elif p.type.value == "rna":
            rna_entry = {"id": p.id, "sequence": p.sequence}
            if p.msa_path:
                rna_entry["unpairedMsaPath"] = p.msa_path
            sequences.append({"rna": rna_entry})

    for lig in job.ligands:
        if lig.repr_type.value == "ccd":
            sequences.append({"ligand": {"id": lig.id, "ccdCodes": lig.ccd_codes}})
        elif lig.repr_type.value == "smiles":
            sequences.append({"ligand": {"id": lig.id, "smiles": lig.smiles}})
        else:
            raise ValueError(
                "AF3 non-server JSON does not accept FILE ligands directly (use CCD or SMILES)."
            )

    for ion in job.ions:
        sequences.append({"ligand": {"id": ion.id, "ccdCodes": [ion.ccd]}})

    glycan_ligands, _ = _expand_glycans(job)
    sequences.extend({"ligand": {"id": l.id, "ccdCodes": l.ccd_codes}} for l in glycan_ligands)

    return sequences


def _bond_key(bond: CovalentBond) -> Tuple[Tuple[str, int, str], Tuple[str, int, str]]:
    a = (bond.a.entity_id, bond.a.position, str(bond.a.atom))
    b = (bond.b.entity_id, bond.b.position, str(bond.b.atom))
    return tuple(sorted([a, b]))  # type: ignore[return-value]


def _render_bonded_pairs(job: JobIR) -> List[list]:
    bonds = list(job.covalent_bonds)
    glycan_ligands, glycan_bonds = _expand_glycans(job)
    if glycan_ligands:
        existing = {_bond_key(b) for b in bonds}
        for b in glycan_bonds:
            key = _bond_key(b)
            if key not in existing:
                bonds.append(b)
                existing.add(key)

    return [
        [[b.a.entity_id, b.a.position, b.a.atom], [b.b.entity_id, b.b.position, b.b.atom]]
        for b in bonds
    ]


def _expand_glycans(job: JobIR) -> Tuple[List[Ligand], List[CovalentBond]]:
    if not job.glycans:
        return [], []
    used_ids = {p.id for p in job.polymers} | {l.id for l in job.ligands} | {i.id for i in job.ions}
    polymer_by_id = {p.id: p for p in job.polymers}
    ligands: List[Ligand] = []
    bonds: List[CovalentBond] = []
    for g in job.glycans:
        root_node = g.attachments[0].root_node if g.attachments else g.nodes[0].node_id
        ligand_id = _unique_id(g.glycan_id, used_ids)
        used_ids.add(ligand_id)

        nodes = {n.node_id: n.ccd for n in g.nodes}
        children = defaultdict(list)
        for e in g.edges:
            children[e.parent].append(e.child)

        order = _topo_tree_order(root_node, children)
        idx_map = {node_id: i + 1 for i, node_id in enumerate(order)}
        ligands.append(
            Ligand(id=ligand_id, repr_type=LigandReprType.ccd, ccd_codes=[nodes[n] for n in order])
        )

        for e in g.edges:
            bonds.append(
                CovalentBond(
                    a=AtomRef(
                        entity_id=ligand_id,
                        position=idx_map[e.parent],
                        atom=e.parent_atom or "O4",
                    ),
                    b=AtomRef(
                        entity_id=ligand_id,
                        position=idx_map[e.child],
                        atom=e.child_atom or "C1",
                    ),
                )
            )

        for att in g.attachments:
            polymer = polymer_by_id.get(att.polymer_id)
            polymer_atom = att.polymer_atom
            if polymer_atom is None and polymer is not None:
                polymer_atom = _default_polymer_atom(polymer, att.polymer_residue_index)
            bonds.append(
                CovalentBond(
                    a=AtomRef(
                        entity_id=att.polymer_id,
                        position=att.polymer_residue_index,
                        atom=polymer_atom or "ND2",
                    ),
                    b=AtomRef(
                        entity_id=ligand_id,
                        position=idx_map[att.root_node],
                        atom=att.root_atom or "C1",
                    ),
                )
            )

    return ligands, bonds


def _topo_tree_order(root: str, children_map: Dict[str, List[str]]) -> List[str]:
    out: List[str] = []
    stack = [root]
    while stack:
        node = stack.pop()
        out.append(node)
        kids = children_map.get(node, [])
        for k in reversed(kids):
            stack.append(k)
    return out


def _unique_id(base: str, used: set) -> str:
    if base not in used:
        return base
    i = 1
    while f"{base}_{i}" in used:
        i += 1
    return f"{base}_{i}"


def _default_polymer_atom(polymer: PolymerChain, residue_index: int) -> str:
    try:
        residue = polymer.sequence[residue_index - 1]
    except Exception:
        return "ND2"
    if polymer.type.value != "protein":
        return "ND2"
    if residue == "N":
        return "ND2"
    if residue == "S":
        return "OG"
    if residue == "T":
        return "OG1"
    return "ND2"


def _validate_sequence_entry(
    entry: dict, loc: str, entity_ids: set, result: ValidationResult
) -> None:
    """Validate a single sequence entry and collect entity IDs."""
    if not isinstance(entry, dict):
        result.add_error("Sequence entry must be an object", loc)
        return

    entry_types = ["protein", "dna", "rna", "ligand"]
    found_type = None
    for t in entry_types:
        if t in entry:
            found_type = t
            break

    if found_type is None:
        result.add_error(
            f"Sequence entry must contain one of: {', '.join(entry_types)}", loc
        )
        return

    data = entry[found_type]
    if not isinstance(data, dict):
        result.add_error(f"'{found_type}' must be an object", loc)
        return

    # Check for required id
    entity_id = data.get("id")
    if entity_id is None:
        result.add_error(f"Missing required 'id' field", f"{loc}.{found_type}")
    else:
        if entity_id in entity_ids:
            result.add_error(f"Duplicate entity ID '{entity_id}'", f"{loc}.{found_type}")
        entity_ids.add(entity_id)

    # Check for required sequence (polymers) or ligand definition
    if found_type in ("protein", "dna", "rna"):
        sequence = data.get("sequence")
        if sequence is None:
            result.add_error(f"Missing required 'sequence' field", f"{loc}.{found_type}")
        elif not isinstance(sequence, str):
            result.add_error(f"'sequence' must be a string", f"{loc}.{found_type}")
        elif len(sequence) == 0:
            result.add_error(f"'sequence' cannot be empty", f"{loc}.{found_type}")

        # Validate modifications
        mods = data.get("modifications") or []
        for mod_idx, mod in enumerate(mods):
            mod_loc = f"{loc}.{found_type}.modifications[{mod_idx}]"
            if not isinstance(mod, dict):
                result.add_error("Modification must be an object", mod_loc)
                continue
            ptm_pos = mod.get("ptmPosition")
            if ptm_pos is None:
                result.add_error("Missing 'ptmPosition'", mod_loc)
            elif sequence and isinstance(ptm_pos, int) and ptm_pos > len(sequence):
                result.add_error(
                    f"ptmPosition {ptm_pos} exceeds sequence length {len(sequence)}",
                    mod_loc,
                )
            if mod.get("ptmType") is None:
                result.add_error("Missing 'ptmType'", mod_loc)
    else:
        # Ligand: must have ccdCodes, smiles, or file
        has_repr = any(k in data for k in ("ccdCodes", "smiles", "file"))
        if not has_repr:
            result.add_error(
                "Ligand must have 'ccdCodes', 'smiles', or 'file'",
                f"{loc}.{found_type}",
            )

from __future__ import annotations

from typing import Dict, List, Tuple

from spir.io.yaml import read_yaml, write_yaml
from spir.ir.models import (
    AtomRef,
    ContactConstraint,
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
    PocketConstraint,
    PolymerChain,
    PolymerType,
)
from spir.validate import ValidationResult


class Boltz2Dialect:
    name = "boltz2"

    def parse(self, path: str) -> DocumentIR:
        payload = read_yaml(path)
        job = _parse_job(payload)
        return DocumentIR(jobs=[job])

    def validate(self, path: str) -> ValidationResult:
        result = ValidationResult()
        try:
            payload = read_yaml(path)
        except Exception as e:
            result.add_error(f"Failed to parse YAML: {e}")
            return result

        if not isinstance(payload, dict):
            result.add_error("Input must be a YAML mapping")
            return result

        # Check version
        version = payload.get("version")
        if version is None:
            result.add_warning("Missing 'version' field (expected 1)")
        elif version != 1:
            result.add_warning(f"Unexpected version {version} (expected 1)")

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
                _validate_boltz_sequence_entry(entry, loc, entity_ids, result)

            # Validate constraints (bonds, contacts, pockets)
            constraints = payload.get("constraints") or []
            for con_idx, con in enumerate(constraints):
                loc = f"constraints[{con_idx}]"
                _validate_boltz_constraint(con, loc, entity_ids, result)

        # Try full parse to catch additional issues
        if result.is_valid:
            try:
                self.parse(path)
            except Exception as e:
                result.add_error(f"Validation passed but parsing failed: {e}")

        return result

    def render(self, doc: DocumentIR, out_path: str) -> None:
        if len(doc.jobs) != 1:
            raise ValueError("Boltz-2 expects exactly one job per YAML file.")
        job = doc.jobs[0]
        payload = _render_job(job)
        write_yaml(out_path, payload)


def _parse_job(payload: dict) -> JobIR:
    name = payload.get("name", "job")
    polymers: List[PolymerChain] = []
    ligands: List[Ligand] = []
    bonds: List[CovalentBond] = []
    ir_constraints: List[object] = []
    constraints = payload.get("constraints") or []

    for entry in payload.get("sequences", []):
        if "protein" in entry:
            data = entry["protein"]
            ids = data.get("id")
            ids_list = ids if isinstance(ids, list) else [ids]
            mods = data.get("modifications") or []
            msa_path = data.get("msa")
            # "empty" is a special Boltz keyword meaning no MSA; treat as None
            if msa_path == "empty":
                msa_path = None
            for chain_id in ids_list:
                polymers.append(
                    PolymerChain(
                        id=chain_id,
                        type=PolymerType.protein,
                        sequence=data["sequence"],
                        modifications=[
                            Modification(position=m["position"], ccd=m["ccd"]) for m in mods
                        ],
                        msa_path=msa_path,
                    )
                )
        elif "dna" in entry:
            data = entry["dna"]
            ids = data.get("id")
            ids_list = ids if isinstance(ids, list) else [ids]
            mods = data.get("modifications") or []
            for chain_id in ids_list:
                polymers.append(
                    PolymerChain(
                        id=chain_id,
                        type=PolymerType.dna,
                        sequence=data["sequence"],
                        modifications=[
                            Modification(position=m["position"], ccd=m["ccd"]) for m in mods
                        ],
                    )
                )
        elif "rna" in entry:
            data = entry["rna"]
            ids = data.get("id")
            ids_list = ids if isinstance(ids, list) else [ids]
            mods = data.get("modifications") or []
            for chain_id in ids_list:
                polymers.append(
                    PolymerChain(
                        id=chain_id,
                        type=PolymerType.rna,
                        sequence=data["sequence"],
                        modifications=[
                            Modification(position=m["position"], ccd=m["ccd"]) for m in mods
                        ],
                    )
                )
        elif "ligand" in entry:
            data = entry["ligand"]
            ids = data.get("id")
            ids_list = ids if isinstance(ids, list) else [ids]
            for chain_id in ids_list:
                if "ccd" in data:
                    ligands.append(
                        Ligand(
                            id=chain_id,
                            repr_type=LigandReprType.ccd,
                            ccd_codes=[data["ccd"]],
                        )
                    )
                elif "smiles" in data:
                    ligands.append(
                        Ligand(
                            id=chain_id,
                            repr_type=LigandReprType.smiles,
                            smiles=data["smiles"],
                        )
                    )

    for item in constraints:
        if "bond" in item:
            bond = item["bond"]
            a = bond["atom1"]
            b = bond["atom2"]
            bonds.append(
                CovalentBond(
                    a=AtomRef(entity_id=a[0], position=int(a[1]), atom=a[2]),
                    b=AtomRef(entity_id=b[0], position=int(b[1]), atom=b[2]),
                )
            )
        elif "contact" in item:
            contact = item["contact"]
            ir_constraints.append(
                ContactConstraint(
                    token1=_parse_token(contact["token1"]),
                    token2=_parse_token(contact["token2"]),
                    max_distance_angstrom=float(contact["max_distance"]),
                )
            )
        elif "pocket" in item:
            pocket = item["pocket"]
            ir_constraints.append(
                PocketConstraint(
                    binder_entity_id=pocket["binder"],
                    contacts=[_parse_token(t) for t in pocket.get("contacts", [])],
                    max_distance_angstrom=float(pocket["max_distance"]),
                )
            )

    job = JobIR(
        name=name,
        polymers=polymers,
        ligands=ligands,
        covalent_bonds=bonds,
        constraints=ir_constraints,
    )
    glycans, filtered_ligands = _detect_glycans(job)
    if glycans:
        job = job.model_copy(update={"glycans": glycans, "ligands": filtered_ligands})
    return job


def _render_job(job: JobIR) -> dict:
    sequences: List[dict] = []

    for p in job.polymers:
        key = p.type.value
        entry = {
            "id": p.id,
            "sequence": p.sequence,
            "modifications": [
                {"position": m.position, "ccd": m.ccd} for m in p.modifications
            ]
            or None,
        }
        # Boltz only supports MSA for proteins
        if p.msa_path and p.type == PolymerType.protein:
            entry["msa"] = p.msa_path
        sequences.append({key: entry})

    for lig in job.ligands:
        if lig.repr_type.value == "ccd":
            sequences.append({"ligand": {"id": lig.id, "ccd": lig.ccd_codes[0]}})
        elif lig.repr_type.value == "smiles":
            sequences.append({"ligand": {"id": lig.id, "smiles": lig.smiles}})

    for ion in job.ions:
        sequences.append({"ligand": {"id": ion.id, "ccd": ion.ccd}})

    glycan_ligands, glycan_bonds = _expand_glycans(job)
    for lig in glycan_ligands:
        sequences.append({"ligand": {"id": lig.id, "ccd": lig.ccd_codes[0]}})

    constraints: List[dict] = []
    bonds = list(job.covalent_bonds)
    bonds.extend(_dedupe_bonds(job.covalent_bonds, glycan_bonds))

    for b in bonds:
        constraints.append(
            {
                "bond": {
                    "atom1": [b.a.entity_id, b.a.position, b.a.atom],
                    "atom2": [b.b.entity_id, b.b.position, b.b.atom],
                }
            }
        )

    for constraint in job.constraints:
        if isinstance(constraint, ContactConstraint):
            constraints.append(
                {
                    "contact": {
                        "token1": _token_from_atomref(constraint.token1),
                        "token2": _token_from_atomref(constraint.token2),
                        "max_distance": constraint.max_distance_angstrom,
                    }
                }
            )
        elif isinstance(constraint, PocketConstraint):
            constraints.append(
                {
                    "pocket": {
                        "binder": constraint.binder_entity_id,
                        "contacts": [_token_from_atomref(a) for a in constraint.contacts],
                        "max_distance": constraint.max_distance_angstrom,
                    }
                }
            )

    payload = {
        "version": 1,
        "sequences": sequences,
        "constraints": constraints or None,
    }
    return payload


def _token_from_atomref(a: AtomRef) -> List:
    if isinstance(a.atom, str):
        return [a.entity_id, a.position, a.atom]
    return [a.entity_id, a.position, a.atom]


def _parse_token(token: List) -> AtomRef:
    if len(token) == 3:
        return AtomRef(entity_id=token[0], position=int(token[1]), atom=token[2])
    if len(token) == 2:
        entity_id = token[0]
        val = token[1]
        if isinstance(val, int):
            return AtomRef(entity_id=entity_id, position=val, atom="CA")
        if isinstance(val, str) and val.isdigit():
            return AtomRef(entity_id=entity_id, position=int(val), atom="CA")
        return AtomRef(entity_id=entity_id, position=1, atom=val)
    raise ValueError(f"Unsupported token format: {token}")


def _detect_glycans(job: JobIR) -> Tuple[List[Glycan], List[Ligand]]:
    ligand_by_id = {
        lig.id: lig
        for lig in job.ligands
        if lig.repr_type == LigandReprType.ccd and len(lig.ccd_codes) == 1
    }
    if not ligand_by_id or not job.covalent_bonds:
        return [], job.ligands

    polymer_ids = {p.id for p in job.polymers}
    edges: List[Tuple[str, str, str, str]] = []
    attachments_by_ligand: Dict[str, List[GlycanAttachment]] = {}

    for bond in job.covalent_bonds:
        a = bond.a
        b = bond.b
        a_is_lig = a.entity_id in ligand_by_id
        b_is_lig = b.entity_id in ligand_by_id
        if a_is_lig and b_is_lig and isinstance(a.atom, str) and isinstance(b.atom, str):
            if a.atom.startswith("O") and b.atom == "C1":
                edges.append((a.entity_id, b.entity_id, a.atom, b.atom))
            elif b.atom.startswith("O") and a.atom == "C1":
                edges.append((b.entity_id, a.entity_id, b.atom, a.atom))
        elif a_is_lig and b.entity_id in polymer_ids and isinstance(a.atom, str):
            if a.atom == "C1":
                attachments_by_ligand.setdefault(a.entity_id, []).append(
                    GlycanAttachment(
                        polymer_id=b.entity_id,
                        polymer_residue_index=b.position,
                        polymer_atom=b.atom if isinstance(b.atom, str) else None,
                        root_node="",
                        root_atom=a.atom,
                    )
                )
        elif b_is_lig and a.entity_id in polymer_ids and isinstance(b.atom, str):
            if b.atom == "C1":
                attachments_by_ligand.setdefault(b.entity_id, []).append(
                    GlycanAttachment(
                        polymer_id=a.entity_id,
                        polymer_residue_index=a.position,
                        polymer_atom=a.atom if isinstance(a.atom, str) else None,
                        root_node="",
                        root_atom=b.atom,
                    )
                )

    glycan_candidates = set()
    for parent_id, child_id, _, _ in edges:
        glycan_candidates.add(parent_id)
        glycan_candidates.add(child_id)
    glycan_candidates.update(attachments_by_ligand.keys())
    if not glycan_candidates:
        return [], job.ligands

    adjacency: Dict[str, List[str]] = {lid: [] for lid in glycan_candidates}
    for parent_id, child_id, _, _ in edges:
        if parent_id in adjacency and child_id in adjacency:
            adjacency[parent_id].append(child_id)
            adjacency[child_id].append(parent_id)

    visited = set()
    glycans: List[Glycan] = []
    glycan_index = 1
    for start in sorted(glycan_candidates):
        if start in visited:
            continue
        stack = [start]
        component = []
        visited.add(start)
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbor in adjacency.get(node, []):
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)

        glycan_id = f"glycan{glycan_index}"
        glycan_index += 1
        node_ids = {lig_id: f"{glycan_id}.n{i}" for i, lig_id in enumerate(component)}
        nodes = [
            GlycanNode(node_id=node_ids[lig_id], ccd=ligand_by_id[lig_id].ccd_codes[0])
            for lig_id in component
        ]
        glycan_edges = [
            GlycanEdge(
                parent=node_ids[parent_id],
                child=node_ids[child_id],
                parent_atom=parent_atom,
                child_atom=child_atom,
            )
            for parent_id, child_id, parent_atom, child_atom in edges
            if parent_id in node_ids and child_id in node_ids
        ]
        attachments: List[GlycanAttachment] = []
        for lig_id in component:
            for att in attachments_by_ligand.get(lig_id, []):
                attachments.append(att.model_copy(update={"root_node": node_ids[lig_id]}))
        glycans.append(
            Glycan(glycan_id=glycan_id, nodes=nodes, edges=glycan_edges, attachments=attachments)
        )

    filtered_ligands = [lig for lig in job.ligands if lig.id not in glycan_candidates]
    return glycans, filtered_ligands


def _expand_glycans(job: JobIR) -> Tuple[List[Ligand], List[CovalentBond]]:
    if not job.glycans:
        return [], []
    used_ids = {p.id for p in job.polymers} | {l.id for l in job.ligands} | {i.id for i in job.ions}
    polymer_by_id = {p.id: p for p in job.polymers}
    ligands: List[Ligand] = []
    bonds: List[CovalentBond] = []
    for g in job.glycans:
        node_to_ligand: Dict[str, str] = {}
        for idx, node in enumerate(g.nodes, start=1):
            lig_id = _unique_id(f"{g.glycan_id}_{idx}", used_ids)
            used_ids.add(lig_id)
            node_to_ligand[node.node_id] = lig_id
            ligands.append(
                Ligand(id=lig_id, repr_type=LigandReprType.ccd, ccd_codes=[node.ccd])
            )
        for e in g.edges:
            bonds.append(
                CovalentBond(
                    a=AtomRef(
                        entity_id=node_to_ligand[e.parent],
                        position=1,
                        atom=e.parent_atom or "O4",
                    ),
                    b=AtomRef(
                        entity_id=node_to_ligand[e.child],
                        position=1,
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
                        entity_id=node_to_ligand[att.root_node],
                        position=1,
                        atom=att.root_atom or "C1",
                    ),
                )
            )
    return ligands, bonds


def _dedupe_bonds(existing: List[CovalentBond], extra: List[CovalentBond]) -> List[CovalentBond]:
    seen = {_bond_key(b) for b in existing}
    out: List[CovalentBond] = []
    for b in extra:
        key = _bond_key(b)
        if key not in seen:
            out.append(b)
            seen.add(key)
    return out


def _bond_key(bond: CovalentBond) -> Tuple[Tuple[str, int, str], Tuple[str, int, str]]:
    a = (bond.a.entity_id, bond.a.position, str(bond.a.atom))
    b = (bond.b.entity_id, bond.b.position, str(bond.b.atom))
    return tuple(sorted([a, b]))  # type: ignore[return-value]


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


def _validate_boltz_sequence_entry(
    entry: dict, loc: str, entity_ids: set, result: ValidationResult
) -> None:
    """Validate a single Boltz sequence entry and collect entity IDs."""
    if not isinstance(entry, dict):
        result.add_error("Sequence entry must be a mapping", loc)
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
        result.add_error(f"'{found_type}' must be a mapping", loc)
        return

    # Check for required id (can be string or list)
    entity_id = data.get("id")
    if entity_id is None:
        result.add_error(f"Missing required 'id' field", f"{loc}.{found_type}")
    else:
        ids_list = entity_id if isinstance(entity_id, list) else [entity_id]
        for eid in ids_list:
            if eid in entity_ids:
                result.add_error(f"Duplicate entity ID '{eid}'", f"{loc}.{found_type}")
            entity_ids.add(eid)

    # Check for required sequence (polymers) or ligand definition
    if found_type in ("protein", "dna", "rna"):
        sequence = data.get("sequence")
        if sequence is None:
            result.add_error(f"Missing required 'sequence' field", f"{loc}.{found_type}")
        elif not isinstance(sequence, str):
            result.add_error(f"'sequence' must be a string", f"{loc}.{found_type}")
        elif len(sequence) == 0:
            result.add_error(f"'sequence' cannot be empty", f"{loc}.{found_type}")
    else:
        # Ligand: must have ccd or smiles
        has_repr = any(k in data for k in ("ccd", "smiles"))
        if not has_repr:
            result.add_error(
                "Ligand must have 'ccd' or 'smiles'",
                f"{loc}.{found_type}",
            )


def _validate_boltz_constraint(
    constraint: dict, loc: str, entity_ids: set, result: ValidationResult
) -> None:
    """Validate a single Boltz constraint."""
    if not isinstance(constraint, dict):
        result.add_error("Constraint must be a mapping", loc)
        return

    constraint_types = ["bond", "contact", "pocket"]
    found_type = None
    for t in constraint_types:
        if t in constraint:
            found_type = t
            break

    if found_type is None:
        result.add_error(
            f"Constraint must contain one of: {', '.join(constraint_types)}", loc
        )
        return

    data = constraint[found_type]
    if not isinstance(data, dict):
        result.add_error(f"'{found_type}' must be a mapping", loc)
        return

    if found_type == "bond":
        for atom_key in ("atom1", "atom2"):
            atom = data.get(atom_key)
            if atom is None:
                result.add_error(f"Missing '{atom_key}'", f"{loc}.{found_type}")
            elif isinstance(atom, list) and len(atom) >= 1:
                entity_id = atom[0]
                if entity_id not in entity_ids:
                    result.add_error(
                        f"Entity '{entity_id}' not found in sequences",
                        f"{loc}.{found_type}.{atom_key}",
                    )
    elif found_type == "contact":
        for token_key in ("token1", "token2"):
            token = data.get(token_key)
            if token is None:
                result.add_error(f"Missing '{token_key}'", f"{loc}.{found_type}")
            elif isinstance(token, list) and len(token) >= 1:
                entity_id = token[0]
                if entity_id not in entity_ids:
                    result.add_error(
                        f"Entity '{entity_id}' not found in sequences",
                        f"{loc}.{found_type}.{token_key}",
                    )
        if data.get("max_distance") is None:
            result.add_error("Missing 'max_distance'", f"{loc}.{found_type}")
    elif found_type == "pocket":
        binder = data.get("binder")
        if binder is None:
            result.add_error("Missing 'binder'", f"{loc}.{found_type}")
        elif binder not in entity_ids:
            result.add_error(
                f"Binder entity '{binder}' not found in sequences",
                f"{loc}.{found_type}",
            )
        if data.get("max_distance") is None:
            result.add_error("Missing 'max_distance'", f"{loc}.{found_type}")

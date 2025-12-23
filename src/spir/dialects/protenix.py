from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Tuple

from spir.io.json import read_json, write_json
from spir.ir.models import (
    AtomRef,
    CovalentBond,
    DocumentIR,
    Glycan,
    Ion,
    JobIR,
    Ligand,
    LigandReprType,
    Modification,
    PolymerChain,
    PolymerType,
)
from spir.validate import ValidationResult


class ProtenixDialect:
    name = "protenix"

    def parse(self, path: str) -> DocumentIR:
        payload = read_json(path)
        jobs_payload = payload
        if isinstance(payload, dict) and "jobs" in payload:
            jobs_payload = payload["jobs"]
        if not isinstance(jobs_payload, list):
            raise ValueError("Protenix input must be a list of jobs.")
        jobs = [_parse_job(job) for job in jobs_payload]
        return DocumentIR(jobs=jobs)

    def validate(self, path: str) -> ValidationResult:
        result = ValidationResult()
        try:
            payload = read_json(path)
        except Exception as e:
            result.add_error(f"Failed to parse JSON: {e}")
            return result

        # Protenix format is a list of jobs (or dict with "jobs" key)
        jobs_payload = payload
        if isinstance(payload, dict) and "jobs" in payload:
            jobs_payload = payload["jobs"]

        if not isinstance(jobs_payload, list):
            result.add_error(
                "Protenix input must be a list of jobs "
                "(or an object with a 'jobs' array)"
            )
            return result

        if len(jobs_payload) == 0:
            result.add_error("Jobs list cannot be empty")
            return result

        for job_idx, job in enumerate(jobs_payload):
            loc = f"jobs[{job_idx}]"
            _validate_protenix_job(job, loc, result)

        # Try full parse to catch additional issues
        if result.is_valid:
            try:
                self.parse(path)
            except Exception as e:
                result.add_error(f"Validation passed but parsing failed: {e}")

        return result

    def render(self, doc: DocumentIR, out_path: str) -> None:
        jobs_payload = [_render_job(job) for job in doc.jobs]
        write_json(out_path, jobs_payload)


def _parse_job(payload: dict) -> JobIR:
    name = payload.get("name", "job")
    polymers: List[PolymerChain] = []
    ligands: List[Ligand] = []
    ions: List[Ion] = []
    bonds: List[CovalentBond] = []

    sequences = payload.get("sequences", [])
    entry_ids: List[List[str]] = []
    entry_types: List[str] = []

    entity_counter = 0
    for entry in sequences:
        entity_counter += 1
        if "proteinChain" in entry:
            p = entry["proteinChain"]
            count = int(p.get("count", 1))
            ids = []
            for _ in range(count):
                chain_id = f"P{entity_counter}_{len(ids) + 1}"
                ids.append(chain_id)
                mods = p.get("modifications") or []
                polymers.append(
                    PolymerChain(
                        id=chain_id,
                        type=PolymerType.protein,
                        sequence=p["sequence"],
                        modifications=[
                            Modification(position=m["ptmPosition"], ccd=m["ptmType"])
                            for m in mods
                        ],
                    )
                )
            entry_ids.append(ids)
            entry_types.append("polymer")
        elif "dnaSequence" in entry:
            p = entry["dnaSequence"]
            count = int(p.get("count", 1))
            ids = []
            for _ in range(count):
                chain_id = f"D{entity_counter}_{len(ids) + 1}"
                ids.append(chain_id)
                mods = p.get("modifications") or []
                polymers.append(
                    PolymerChain(
                        id=chain_id,
                        type=PolymerType.dna,
                        sequence=p["sequence"],
                        modifications=[
                            Modification(position=m["basePosition"], ccd=m["modificationType"])
                            for m in mods
                        ],
                    )
                )
            entry_ids.append(ids)
            entry_types.append("polymer")
        elif "rnaSequence" in entry:
            p = entry["rnaSequence"]
            count = int(p.get("count", 1))
            ids = []
            for _ in range(count):
                chain_id = f"R{entity_counter}_{len(ids) + 1}"
                ids.append(chain_id)
                mods = p.get("modifications") or []
                polymers.append(
                    PolymerChain(
                        id=chain_id,
                        type=PolymerType.rna,
                        sequence=p["sequence"],
                        modifications=[
                            Modification(position=m["basePosition"], ccd=m["modificationType"])
                            for m in mods
                        ],
                    )
                )
            entry_ids.append(ids)
            entry_types.append("polymer")
        elif "ligand" in entry:
            l = entry["ligand"]
            count = int(l.get("count", 1))
            ids = []
            for _ in range(count):
                lig_id = f"L{entity_counter}_{len(ids) + 1}"
                ids.append(lig_id)
                ligands.append(_parse_ligand(lig_id, l["ligand"]))
            entry_ids.append(ids)
            entry_types.append("ligand")
        elif "ion" in entry:
            i = entry["ion"]
            count = int(i.get("count", 1))
            ids = []
            for _ in range(count):
                ion_id = f"I{entity_counter}_{len(ids) + 1}"
                ids.append(ion_id)
                ions.append(Ion(id=ion_id, ccd=i["ion"]))
            entry_ids.append(ids)
            entry_types.append("ion")

    for bond in payload.get("covalent_bonds") or []:
        entity1 = int(bond.get("entity1") or bond.get("left_entity"))
        entity2 = int(bond.get("entity2") or bond.get("right_entity"))
        copy1 = int(bond.get("copy1") or bond.get("left_copy") or 1)
        copy2 = int(bond.get("copy2") or bond.get("right_copy") or 1)
        position1 = int(bond.get("position1") or bond.get("left_position"))
        position2 = int(bond.get("position2") or bond.get("right_position"))
        atom1 = bond.get("atom1") or bond.get("left_atom")
        atom2 = bond.get("atom2") or bond.get("right_atom")

        entity1_ids = entry_ids[entity1 - 1] if entity1 - 1 < len(entry_ids) else []
        entity2_ids = entry_ids[entity2 - 1] if entity2 - 1 < len(entry_ids) else []
        if not entity1_ids or not entity2_ids:
            continue
        if copy1 - 1 >= len(entity1_ids) or copy2 - 1 >= len(entity2_ids):
            continue
        bonds.append(
            CovalentBond(
                a=AtomRef(
                    entity_id=entity1_ids[copy1 - 1],
                    position=position1,
                    atom=_parse_atom(atom1),
                ),
                b=AtomRef(
                    entity_id=entity2_ids[copy2 - 1],
                    position=position2,
                    atom=_parse_atom(atom2),
                ),
            )
        )

    return JobIR(
        name=name,
        polymers=polymers,
        ligands=ligands,
        ions=ions,
        covalent_bonds=bonds,
    )


def _parse_ligand(lig_id: str, ligand: str) -> Ligand:
    if ligand.startswith("CCD_"):
        codes = ligand[4:].split("_")
        return Ligand(id=lig_id, repr_type=LigandReprType.ccd, ccd_codes=codes)
    if ligand.startswith("FILE_"):
        return Ligand(id=lig_id, repr_type=LigandReprType.file, file_path=ligand[5:])
    if ligand.isalnum() and ligand.upper() == ligand:
        return Ligand(id=lig_id, repr_type=LigandReprType.ccd, ccd_codes=[ligand])
    return Ligand(id=lig_id, repr_type=LigandReprType.smiles, smiles=ligand)


def _parse_atom(atom) -> object:
    if atom is None:
        return ""
    if isinstance(atom, int):
        return atom
    if isinstance(atom, str) and atom.isdigit():
        return int(atom)
    return atom


def _render_job(job: JobIR) -> dict:
    sequences: List[dict] = []
    entity_ids: List[str] = []
    entity_map: Dict[str, int] = {}

    for p in job.polymers:
        if p.type.value == "protein":
            key = "proteinChain"
            mods = [
                {"ptmType": _prefix_ccd(m.ccd), "ptmPosition": m.position}
                for m in p.modifications
            ]
        elif p.type.value == "dna":
            key = "dnaSequence"
            mods = [
                {"modificationType": _prefix_ccd(m.ccd), "basePosition": m.position}
                for m in p.modifications
            ]
        else:
            key = "rnaSequence"
            mods = [
                {"modificationType": _prefix_ccd(m.ccd), "basePosition": m.position}
                for m in p.modifications
            ]
        sequences.append(
            {
                key: {
                    "sequence": p.sequence,
                    "modifications": mods or None,
                    "count": 1,
                }
            }
        )
        entity_ids.append(p.id)
    for lig in job.ligands:
        sequences.append({"ligand": {"ligand": _ligand_string(lig), "count": 1}})
        entity_ids.append(lig.id)
    for ion in job.ions:
        sequences.append({"ion": {"ion": ion.ccd, "count": 1}})
        entity_ids.append(ion.id)

    glycan_ligands, glycan_bonds = _expand_glycans(job, entity_ids)
    for lig in glycan_ligands:
        sequences.append({"ligand": {"ligand": _ligand_string(lig), "count": 1}})
        entity_ids.append(lig.id)

    for idx, entity_id in enumerate(entity_ids, start=1):
        entity_map[entity_id] = idx

    bonds = list(job.covalent_bonds)
    bonds.extend(_dedupe_bonds(job.covalent_bonds, glycan_bonds))
    covalent_bonds = [_render_bond(b, entity_map) for b in bonds if _render_bond(b, entity_map)]

    return {
        "name": job.name,
        "sequences": sequences,
        "covalent_bonds": covalent_bonds or None,
    }


def _render_bond(bond: CovalentBond, entity_map: Dict[str, int]) -> dict | None:
    if bond.a.entity_id not in entity_map or bond.b.entity_id not in entity_map:
        return None
    return {
        "entity1": str(entity_map[bond.a.entity_id]),
        "copy1": 1,
        "position1": str(bond.a.position),
        "atom1": bond.a.atom,
        "entity2": str(entity_map[bond.b.entity_id]),
        "copy2": 1,
        "position2": str(bond.b.position),
        "atom2": bond.b.atom,
    }


def _ligand_string(lig: Ligand) -> str:
    if lig.repr_type.value == "ccd":
        return _prefix_ccd("_".join(lig.ccd_codes))
    if lig.repr_type.value == "file":
        return f"FILE_{lig.file_path}"
    if lig.repr_type.value == "smiles":
        return lig.smiles or ""
    return ""


def _prefix_ccd(code: str) -> str:
    return code if code.startswith("CCD_") else f"CCD_{code}"


def _expand_glycans(job: JobIR, entity_ids: List[str]) -> Tuple[List[Ligand], List[CovalentBond]]:
    if not job.glycans:
        return [], []
    used_ids = set(entity_ids)
    polymer_by_id = {p.id: p for p in job.polymers}
    ligands: List[Ligand] = []
    bonds: List[CovalentBond] = []
    for g in job.glycans:
        root_node = g.attachments[0].root_node if g.attachments else g.nodes[0].node_id
        nodes = {n.node_id: n.ccd for n in g.nodes}
        children = defaultdict(list)
        for e in g.edges:
            children[e.parent].append(e.child)
        order = _topo_tree_order(root_node, children)
        idx_map = {node_id: i + 1 for i, node_id in enumerate(order)}
        lig_id = _unique_id(g.glycan_id, used_ids)
        used_ids.add(lig_id)
        ligands.append(
            Ligand(id=lig_id, repr_type=LigandReprType.ccd, ccd_codes=[nodes[n] for n in order])
        )
        for e in g.edges:
            bonds.append(
                CovalentBond(
                    a=AtomRef(
                        entity_id=lig_id,
                        position=idx_map[e.parent],
                        atom=e.parent_atom or "O4",
                    ),
                    b=AtomRef(
                        entity_id=lig_id,
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
                        entity_id=lig_id,
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


def _validate_protenix_job(job: dict, loc: str, result: ValidationResult) -> None:
    """Validate a single Protenix job."""
    if not isinstance(job, dict):
        result.add_error("Job must be an object", loc)
        return

    # Check sequences
    sequences = job.get("sequences")
    if sequences is None:
        result.add_error("Missing required field 'sequences'", loc)
        return
    if not isinstance(sequences, list):
        result.add_error("'sequences' must be a list", loc)
        return

    entity_count = 0
    for seq_idx, entry in enumerate(sequences):
        seq_loc = f"{loc}.sequences[{seq_idx}]"
        entity_count += 1
        _validate_protenix_sequence_entry(entry, seq_loc, result)

    # Validate covalent_bonds
    bonds = job.get("covalent_bonds") or []
    for bond_idx, bond in enumerate(bonds):
        bond_loc = f"{loc}.covalent_bonds[{bond_idx}]"
        _validate_protenix_bond(bond, bond_loc, entity_count, result)


def _validate_protenix_sequence_entry(
    entry: dict, loc: str, result: ValidationResult
) -> None:
    """Validate a single Protenix sequence entry."""
    if not isinstance(entry, dict):
        result.add_error("Sequence entry must be an object", loc)
        return

    entry_types = ["proteinChain", "dnaSequence", "rnaSequence", "ligand", "ion"]
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

    # Check for required sequence (polymers) or ligand/ion definition
    if found_type in ("proteinChain", "dnaSequence", "rnaSequence"):
        sequence = data.get("sequence")
        if sequence is None:
            result.add_error(f"Missing required 'sequence' field", f"{loc}.{found_type}")
        elif not isinstance(sequence, str):
            result.add_error(f"'sequence' must be a string", f"{loc}.{found_type}")
        elif len(sequence) == 0:
            result.add_error(f"'sequence' cannot be empty", f"{loc}.{found_type}")

        # Validate modifications
        if found_type == "proteinChain":
            mods = data.get("modifications") or []
            for mod_idx, mod in enumerate(mods):
                mod_loc = f"{loc}.{found_type}.modifications[{mod_idx}]"
                if not isinstance(mod, dict):
                    result.add_error("Modification must be an object", mod_loc)
                    continue
                if mod.get("ptmPosition") is None:
                    result.add_error("Missing 'ptmPosition'", mod_loc)
                if mod.get("ptmType") is None:
                    result.add_error("Missing 'ptmType'", mod_loc)
        else:
            mods = data.get("modifications") or []
            for mod_idx, mod in enumerate(mods):
                mod_loc = f"{loc}.{found_type}.modifications[{mod_idx}]"
                if not isinstance(mod, dict):
                    result.add_error("Modification must be an object", mod_loc)
                    continue
                if mod.get("basePosition") is None:
                    result.add_error("Missing 'basePosition'", mod_loc)
                if mod.get("modificationType") is None:
                    result.add_error("Missing 'modificationType'", mod_loc)
    elif found_type == "ligand":
        if data.get("ligand") is None:
            result.add_error("Missing 'ligand' identifier", f"{loc}.{found_type}")
    elif found_type == "ion":
        if data.get("ion") is None:
            result.add_error("Missing 'ion' CCD code", f"{loc}.{found_type}")


def _validate_protenix_bond(
    bond: dict, loc: str, entity_count: int, result: ValidationResult
) -> None:
    """Validate a single Protenix covalent bond."""
    if not isinstance(bond, dict):
        result.add_error("Bond must be an object", loc)
        return

    # Check entity references
    for side in ("1", "2"):
        entity_key = f"entity{side}"
        alt_entity_key = f"{'left' if side == '1' else 'right'}_entity"
        entity = bond.get(entity_key) or bond.get(alt_entity_key)
        if entity is None:
            result.add_error(f"Missing entity reference ('{entity_key}')", loc)
        else:
            try:
                entity_num = int(entity)
                if entity_num < 1 or entity_num > entity_count:
                    result.add_error(
                        f"Entity {entity_num} out of range (1-{entity_count})", loc
                    )
            except (TypeError, ValueError):
                result.add_error(f"Entity reference must be an integer", loc)

        pos_key = f"position{side}"
        alt_pos_key = f"{'left' if side == '1' else 'right'}_position"
        position = bond.get(pos_key) or bond.get(alt_pos_key)
        if position is None:
            result.add_error(f"Missing position ('{pos_key}')", loc)

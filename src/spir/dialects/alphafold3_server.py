from __future__ import annotations

from typing import Dict, List

from spir.io.json import read_json, write_json
from spir.ir.glycans.parse_af3_server import parse_af3_server_glycan_string
from spir.ir.glycans.render_af3_server import render_af3_server_glycan_string
from spir.ir.models import (
    DocumentIR,
    Glycan,
    GlycanAttachment,
    Ion,
    JobIR,
    Ligand,
    LigandReprType,
    Modification,
    PolymerChain,
    PolymerType,
)
from spir.validate import ValidationResult


class AlphaFold3ServerDialect:
    name = "alphafold3server"

    def parse(self, path: str) -> DocumentIR:
        payload = read_json(path)
        jobs_payload = payload
        if isinstance(payload, dict) and "jobs" in payload:
            jobs_payload = payload["jobs"]
        if not isinstance(jobs_payload, list):
            raise ValueError("AlphaFold3 Server input must be a list of jobs.")
        jobs = [_parse_job(job) for job in jobs_payload]
        return DocumentIR(jobs=jobs)

    def validate(self, path: str) -> ValidationResult:
        result = ValidationResult()
        try:
            payload = read_json(path)
        except Exception as e:
            result.add_error(f"Failed to parse JSON: {e}")
            return result

        # AF3 Server format is a list of jobs (or dict with "jobs" key)
        jobs_payload = payload
        if isinstance(payload, dict) and "jobs" in payload:
            jobs_payload = payload["jobs"]

        if not isinstance(jobs_payload, list):
            result.add_error(
                "AlphaFold3 Server input must be a list of jobs "
                "(or an object with a 'jobs' array)"
            )
            return result

        if len(jobs_payload) == 0:
            result.add_error("Jobs list cannot be empty")
            return result

        for job_idx, job in enumerate(jobs_payload):
            loc = f"jobs[{job_idx}]"
            _validate_server_job(job, loc, result)

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
    seeds = payload.get("modelSeeds") or []
    polymers: List[PolymerChain] = []
    ligands: List[Ligand] = []
    ions: List[Ion] = []
    glycans: List[Glycan] = []

    chain_counter = 0
    ligand_counter = 0
    ion_counter = 0

    for entry in payload.get("sequences", []):
        if "proteinChain" in entry:
            p = entry["proteinChain"]
            count = int(p.get("count", 1))
            mods = p.get("modifications") or []
            glycan_entries = p.get("glycans") or []
            # Non-standard extension: allow msa_path for conversion to other formats
            msa_path = p.get("msa_path")
            for _ in range(count):
                chain_counter += 1
                chain_id = _chain_id(chain_counter)
                polymers.append(
                    PolymerChain(
                        id=chain_id,
                        type=PolymerType.protein,
                        sequence=p["sequence"],
                        modifications=[
                            Modification(position=m["ptmPosition"], ccd=m["ptmType"]) for m in mods
                        ],
                        msa_path=msa_path,
                    )
                )
                for g_idx, g_entry in enumerate(glycan_entries):
                    glycan_id = f"{chain_id}_glycan{g_idx + 1}"
                    glycan = parse_af3_server_glycan_string(glycan_id, g_entry["residues"])
                    glycan = glycan.model_copy(
                        update={
                            "attachments": [
                                GlycanAttachment(
                                    polymer_id=chain_id,
                                    polymer_residue_index=int(g_entry["position"]),
                                    polymer_atom=None,
                                    root_node=f"{glycan_id}.n0",
                                    root_atom="C1",
                                )
                            ]
                        }
                    )
                    glycans.append(glycan)
        elif "dnaSequence" in entry:
            p = entry["dnaSequence"]
            count = int(p.get("count", 1))
            mods = p.get("modifications") or []
            # Non-standard extension: allow msa_path for conversion to other formats
            msa_path = p.get("msa_path")
            for _ in range(count):
                chain_counter += 1
                chain_id = _chain_id(chain_counter)
                polymers.append(
                    PolymerChain(
                        id=chain_id,
                        type=PolymerType.dna,
                        sequence=p["sequence"],
                        modifications=[
                            Modification(position=m["basePosition"], ccd=m["modificationType"])
                            for m in mods
                        ],
                        msa_path=msa_path,
                    )
                )
        elif "rnaSequence" in entry:
            p = entry["rnaSequence"]
            count = int(p.get("count", 1))
            mods = p.get("modifications") or []
            # Non-standard extension: allow msa_path for conversion to other formats
            msa_path = p.get("msa_path")
            for _ in range(count):
                chain_counter += 1
                chain_id = _chain_id(chain_counter)
                polymers.append(
                    PolymerChain(
                        id=chain_id,
                        type=PolymerType.rna,
                        sequence=p["sequence"],
                        modifications=[
                            Modification(position=m["basePosition"], ccd=m["modificationType"])
                            for m in mods
                        ],
                        msa_path=msa_path,
                    )
                )
        elif "ligand" in entry:
            l = entry["ligand"]
            count = int(l.get("count", 1))
            for _ in range(count):
                ligand_counter += 1
                ligands.append(
                    Ligand(
                        id=f"L{ligand_counter}",
                        repr_type=LigandReprType.ccd,
                        ccd_codes=[l["ligand"]],
                    )
                )
        elif "ion" in entry:
            i = entry["ion"]
            count = int(i.get("count", 1))
            for _ in range(count):
                ion_counter += 1
                ions.append(Ion(id=f"I{ion_counter}", ccd=i["ion"]))

    return JobIR(
        name=name,
        seeds=seeds,
        polymers=polymers,
        ligands=ligands,
        glycans=glycans,
        ions=ions,
    )


def _render_job(job: JobIR) -> dict:
    sequences: List[dict] = []
    glycans_by_polymer: Dict[str, List[Glycan]] = {}
    for g in job.glycans:
        for att in g.attachments:
            glycans_by_polymer.setdefault(att.polymer_id, []).append(g)

    for p in job.polymers:
        if p.type.value == "protein":
            glycan_entries = []
            for g in glycans_by_polymer.get(p.id, []):
                for att in g.attachments:
                    if att.polymer_id != p.id:
                        continue
                    residues = render_af3_server_glycan_string(g, att.root_node)
                    glycan_entries.append(
                        {
                            "residues": residues,
                            "position": att.polymer_residue_index,
                        }
                    )
            sequences.append(
                {
                    "proteinChain": {
                        "sequence": p.sequence,
                        "modifications": [
                            {
                                "ptmType": _prefix_ccd(m.ccd),
                                "ptmPosition": m.position,
                            }
                            for m in p.modifications
                        ]
                        or None,
                        "glycans": glycan_entries or None,
                        "count": 1,
                    }
                }
            )
        elif p.type.value == "dna":
            sequences.append(
                {
                    "dnaSequence": {
                        "sequence": p.sequence,
                        "modifications": [
                            {
                                "modificationType": _prefix_ccd(m.ccd),
                                "basePosition": m.position,
                            }
                            for m in p.modifications
                        ]
                        or None,
                        "count": 1,
                    }
                }
            )
        elif p.type.value == "rna":
            sequences.append(
                {
                    "rnaSequence": {
                        "sequence": p.sequence,
                        "modifications": [
                            {
                                "modificationType": _prefix_ccd(m.ccd),
                                "basePosition": m.position,
                            }
                            for m in p.modifications
                        ]
                        or None,
                        "count": 1,
                    }
                }
            )

    for lig in job.ligands:
        if lig.repr_type.value != "ccd" or not lig.ccd_codes:
            continue
        sequences.append({"ligand": {"ligand": _prefix_ccd(lig.ccd_codes[0]), "count": 1}})

    for ion in job.ions:
        sequences.append({"ion": {"ion": _prefix_ccd(ion.ccd), "count": 1}})

    return {
        "name": job.name,
        "modelSeeds": job.seeds,
        "sequences": sequences,
        "dialect": "alphafoldserver",
        "version": 1,
    }


def _chain_id(idx: int) -> str:
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    out = ""
    n = idx
    while n > 0:
        n -= 1
        out = alphabet[n % 26] + out
        n //= 26
    return out


def _prefix_ccd(code: str) -> str:
    return code if code.startswith("CCD_") else f"CCD_{code}"


def _validate_server_job(job: dict, loc: str, result: ValidationResult) -> None:
    """Validate a single AlphaFold Server job."""
    if not isinstance(job, dict):
        result.add_error("Job must be an object", loc)
        return

    # Check sequences
    sequences = job.get("sequences")
    if sequences is None:
        result.add_error("Missing required field 'sequences'", loc)
    elif not isinstance(sequences, list):
        result.add_error("'sequences' must be a list", loc)
    else:
        for seq_idx, entry in enumerate(sequences):
            seq_loc = f"{loc}.sequences[{seq_idx}]"
            _validate_server_sequence_entry(entry, seq_loc, result)


def _validate_server_sequence_entry(entry: dict, loc: str, result: ValidationResult) -> None:
    """Validate a single AlphaFold Server sequence entry."""
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
        result.add_error(f"Sequence entry must contain one of: {', '.join(entry_types)}", loc)
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

            # Validate glycans
            glycans = data.get("glycans") or []
            for glyc_idx, glyc in enumerate(glycans):
                glyc_loc = f"{loc}.{found_type}.glycans[{glyc_idx}]"
                if not isinstance(glyc, dict):
                    result.add_error("Glycan must be an object", glyc_loc)
                    continue
                if glyc.get("residues") is None:
                    result.add_error("Missing 'residues'", glyc_loc)
                if glyc.get("position") is None:
                    result.add_error("Missing 'position'", glyc_loc)
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
            result.add_error("Missing 'ligand' CCD code", f"{loc}.{found_type}")
    elif found_type == "ion":
        if data.get("ion") is None:
            result.add_error("Missing 'ion' CCD code", f"{loc}.{found_type}")

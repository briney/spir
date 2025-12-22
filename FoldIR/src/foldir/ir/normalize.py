from __future__ import annotations

from typing import Iterable, List, Optional

from foldir.ir.glycans.resolve_linkages import DefaultSugarLinkageResolver, fill_missing_linkages
from foldir.ir.ids import ensure_unique_entity_ids, ensure_unique_glycan_ids
from foldir.ir.models import DocumentIR, Glycan, GlycanAttachment, GlycanNode, JobIR, Ligand, Modification


def normalize_ccd(code: str) -> str:
    return code[4:] if code.startswith("CCD_") else code


def _normalize_mods(mods: Iterable[Modification]) -> List[Modification]:
    out: List[Modification] = []
    for m in mods:
        out.append(m.model_copy(update={"ccd": normalize_ccd(m.ccd)}))
    return out


def _normalize_ligands(ligands: Iterable[Ligand]) -> List[Ligand]:
    out: List[Ligand] = []
    for lig in ligands:
        if lig.repr_type.value == "ccd":
            out.append(lig.model_copy(update={"ccd_codes": [normalize_ccd(c) for c in lig.ccd_codes]}))
        else:
            out.append(lig)
    return out


def _normalize_glycan_nodes(nodes: Iterable[GlycanNode]) -> List[GlycanNode]:
    return [n.model_copy(update={"ccd": normalize_ccd(n.ccd)}) for n in nodes]


def normalize_glycan(g: Glycan) -> Glycan:
    return g.model_copy(update={"nodes": _normalize_glycan_nodes(g.nodes)})


def _default_polymer_atom(polymer, residue_index: int, opts: object) -> str:
    default_asn = getattr(opts, "default_asn_atom", "ND2")
    default_ser = getattr(opts, "default_ser_atom", "OG")
    default_thr = getattr(opts, "default_thr_atom", "OG1")
    try:
        residue = polymer.sequence[residue_index - 1]
    except Exception:
        return default_asn
    if polymer.type.value != "protein":
        return default_asn
    if residue == "N":
        return default_asn
    if residue == "S":
        return default_ser
    if residue == "T":
        return default_thr
    return default_asn


def _fill_glycan_defaults(job: JobIR, opts: object) -> JobIR:
    resolver = DefaultSugarLinkageResolver(
        default_parent_atom=getattr(opts, "default_glycan_parent_atom", "O4"),
        default_child_atom=getattr(opts, "default_glycan_child_atom", "C1"),
    )
    polymer_by_id = {p.id: p for p in job.polymers}
    new_glycans = []
    for g in job.glycans:
        g_filled = fill_missing_linkages(g, resolver)
        attachments = []
        for att in g_filled.attachments:
            polymer_atom = att.polymer_atom
            if polymer_atom is None:
                polymer = polymer_by_id.get(att.polymer_id)
                polymer_atom = _default_polymer_atom(polymer, att.polymer_residue_index, opts)
            root_atom = att.root_atom or getattr(opts, "default_glycan_child_atom", "C1")
            attachments.append(
                GlycanAttachment(
                    polymer_id=att.polymer_id,
                    polymer_residue_index=att.polymer_residue_index,
                    polymer_atom=polymer_atom,
                    root_node=att.root_node,
                    root_atom=root_atom,
                )
            )
        g_filled = g_filled.model_copy(update={"attachments": attachments})
        new_glycans.append(g_filled)
    return job.model_copy(update={"glycans": new_glycans})


def normalize_job(job: JobIR, opts: Optional[object] = None) -> JobIR:
    polymers = [p.model_copy(update={"modifications": _normalize_mods(p.modifications)}) for p in job.polymers]
    ligands = _normalize_ligands(job.ligands)
    ions = [i.model_copy(update={"ccd": normalize_ccd(i.ccd)}) for i in job.ions]
    glycans = [normalize_glycan(g) for g in job.glycans]

    polymers, ligands, ions = ensure_unique_entity_ids(polymers, ligands, ions)
    glycans = ensure_unique_glycan_ids(glycans)

    new_job = job.model_copy(
        update={
            "polymers": polymers,
            "ligands": ligands,
            "ions": ions,
            "glycans": glycans,
        }
    )
    if opts is not None:
        new_job = _fill_glycan_defaults(new_job, opts)
    return new_job


def normalize_document(doc: DocumentIR, opts: Optional[object] = None) -> DocumentIR:
    jobs = [normalize_job(j, opts=opts) for j in doc.jobs]
    return doc.model_copy(update={"jobs": jobs})

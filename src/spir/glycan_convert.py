from __future__ import annotations

import json
from typing import List, Literal, Tuple

from .adapters.af3 import AF3ReaderWriter
from .adapters.boltz import BoltzReaderWriter
from .adapters.chai import ChaiReaderWriter
from .adapters.protenix import ProtenixReaderWriter
from .ir.glycan import (
    GlycanGraph,
    parse_af3_server_glycan,
    render_chai_glycan,
)
from .ir.model import AtomRef, Bond, ComplexInput, Ligand, PolymerChain


def _parse_af3_server_components(glycan_residues_spec: str) -> List[str]:
    """Extract component CCD codes from an AF3 Server glycan residues string.

    This mirrors the lightweight extraction used in the AF3 Server adapter:
    take tokens that look like CCD codes (alnum/underscore/hyphen), ignoring
    parentheses and link annotations.
    """
    comp: List[str] = []
    buf: List[str] = []
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
    # Filter out empty tokens
    comp = [c for c in comp if c]
    if not comp:
        raise ValueError("No CCD components could be parsed from glycan string")
    return comp


def _full_text_for_format(components: List[str], target: str) -> str:
    ci = ComplexInput(
        name="spir-glycan",
        ligands=[Ligand(ids=["G"], ccd_codes=list(components))],
    )
    if target in {"af3", "alphafold3"}:
        return AF3ReaderWriter.dump(ci)
    if target == "boltz":
        return BoltzReaderWriter.dump(ci)
    if target == "protenix":
        return ProtenixReaderWriter.dump(ci)
    if target == "chai":
        fasta_text, restraints_text = ChaiReaderWriter.dump(ci)
        # If no restraints, include an empty header to aid authoring
        if restraints_text is None:
            header = (
                "chainA,res_idxA,chainB,res_idxB,connection_type,confidence,"
                "min_distance_angstrom,max_distance_angstrom,comment,restraint_id\n"
            )
            return f"{fasta_text}\n# restraints.csv (empty)\n{header}"
        return f"{fasta_text}\n---- restraints.csv ----\n{restraints_text}"
    if target == "af3-server":
        job = [
            {
                "name": "spir-job",
                "modelSeeds": [],
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": "",  # illustrative placeholder
                            "count": 1,
                            "glycans": [{"residues": "_PLACEHOLDER_", "position": 1}],
                        }
                    }
                ],
                "dialect": "alphafoldserver",
                "version": 1,
            }
        ]
        # Fill actual residues string in the illustrative slot for clarity
        job[0]["sequences"][0]["proteinChain"]["glycans"][0]["residues"] = (
            "_FROM_INPUT_"
        )
        return json.dumps(job, indent=2)
    raise ValueError(f"Unsupported format: {target}")


def _graph_to_ci_for_af3(graph: GlycanGraph) -> ComplexInput:
    # Protein mock anchor
    protein = PolymerChain(type="protein", ids=["A"], sequence="N")
    ligand = Ligand(ids=["G"], ccd_codes=graph.components_as_ccd_list())
    bonds: List[Bond] = []
    # Intra-glycan bonds
    for b in graph.bonds:
        bonds.append(
            Bond(
                atom1=AtomRef(
                    chain_id="G",
                    residue_index=b.parent_index,
                    component_index=b.parent_index,
                    atom_name=b.parent_atom,
                ),
                atom2=AtomRef(
                    chain_id="G",
                    residue_index=b.child_index,
                    component_index=b.child_index,
                    atom_name=b.child_atom,
                ),
            )
        )
    # Anchor: Asn ND2 ↔ root C1 (root is component 1)
    bonds.append(
        Bond(
            atom1=AtomRef(chain_id="A", residue_index=1, atom_name="ND2"),
            atom2=AtomRef(
                chain_id="G", residue_index=1, component_index=1, atom_name="C1"
            ),
        )
    )
    return ComplexInput(
        name="spir-glycan", proteins=[protein], ligands=[ligand], bonds=bonds
    )


def _graph_to_ci_for_protenix(graph: GlycanGraph) -> ComplexInput:
    # Same structure as AF3, but we will dump with Protenix writer
    protein = PolymerChain(type="protein", ids=["A"], sequence="N")
    ligand = Ligand(ids=["G"], ccd_codes=graph.components_as_ccd_list())
    bonds: List[Bond] = []
    for b in graph.bonds:
        bonds.append(
            Bond(
                atom1=AtomRef(
                    chain_id="G",
                    residue_index=b.parent_index,
                    component_index=b.parent_index,
                    atom_name=b.parent_atom,
                ),
                atom2=AtomRef(
                    chain_id="G",
                    residue_index=b.child_index,
                    component_index=b.child_index,
                    atom_name=b.child_atom,
                ),
            )
        )
    bonds.append(
        Bond(
            atom1=AtomRef(chain_id="A", residue_index=1, atom_name="ND2"),
            atom2=AtomRef(
                chain_id="G", residue_index=1, component_index=1, atom_name="C1"
            ),
        )
    )
    return ComplexInput(
        name="spir-glycan", proteins=[protein], ligands=[ligand], bonds=bonds
    )


def _graph_to_ci_for_boltz(graph: GlycanGraph) -> ComplexInput:
    # Represent each sugar as a separate ligand entity (G1..Gn), residue_index is 1 for each
    n = len(graph.components)
    ligands: List[Ligand] = []
    for i in range(n):
        ligands.append(Ligand(ids=[f"G{i+1}"], ccd_codes=[graph.components[i].ccd]))
    protein = PolymerChain(type="protein", ids=["A"], sequence="N")
    bonds: List[Bond] = []
    for b in graph.bonds:
        bonds.append(
            Bond(
                atom1=AtomRef(
                    chain_id=f"G{b.parent_index}",
                    residue_index=1,
                    atom_name=b.parent_atom,
                ),
                atom2=AtomRef(
                    chain_id=f"G{b.child_index}",
                    residue_index=1,
                    atom_name=b.child_atom,
                ),
            )
        )
    bonds.append(
        Bond(
            atom1=AtomRef(chain_id="A", residue_index=1, atom_name="ND2"),
            atom2=AtomRef(chain_id="G1", residue_index=1, atom_name="C1"),
        )
    )
    return ComplexInput(
        name="spir-glycan", proteins=[protein], ligands=ligands, bonds=bonds
    )


def _graph_to_fasta_and_restraints_for_chai(graph: GlycanGraph) -> Tuple[str, str]:
    # Build FASTA with a protein and a glycan entry
    glycan_spec = render_chai_glycan(graph)
    fasta_lines = [
        ">protein|A",
        "N",
        ">glycan|G",
        glycan_spec,
    ]
    fasta_text = "\n".join(fasta_lines) + "\n"
    # One covalent row anchoring A:N1@N to glycan root @C1
    header = (
        "chainA,res_idxA,chainB,res_idxB,connection_type,confidence,"
        "min_distance_angstrom,max_distance_angstrom,comment,restraint_id\n"
    )
    row = "A,N1@N,G,@C1,covalent,1.0,0.0,0.0,protein-glycan,bond1\n"
    restraints_text = header + row
    return fasta_text, restraints_text


def convert_af3_server_glycan(
    glycan: str,
    to: Literal[
        "alphafold3", "af3", "af3-server", "boltz", "chai", "protenix"
    ] = "alphafold3",
    detail: Literal["minimal", "full"] = "minimal",
) -> str:
    """Convert an AF3 Server glycan string into a per-format view/snippet.

    - When detail="minimal": returns a small, copy-pasteable snippet showing how
      the glycan would be represented in the target format.
    - When detail="full": returns a complete example file text where applicable.
    """
    target = "alphafold3" if to == "af3" else to
    components = _parse_af3_server_components(glycan)

    if detail == "full":
        # Build a fully connected representation using deterministic linkage inference
        graph = parse_af3_server_glycan(glycan)
        if target == "alphafold3":
            ci = _graph_to_ci_for_af3(graph)
            return AF3ReaderWriter.dump(ci)
        if target == "protenix":
            ci = _graph_to_ci_for_protenix(graph)
            return ProtenixReaderWriter.dump(ci)
        if target == "boltz":
            ci = _graph_to_ci_for_boltz(graph)
            return BoltzReaderWriter.dump(ci)
        if target == "chai":
            fasta_text, restraints_text = _graph_to_fasta_and_restraints_for_chai(graph)
            return f"{fasta_text}\n---- restraints.csv ----\n{restraints_text}"
        if target == "af3-server":
            # Echo back original glycan spec embedded in a minimal AF3-Server job with mock protein
            # This keeps focus on conversion to other formats; server input is just reference.
            return _full_text_for_format(components, target)
        raise ValueError(f"Unsupported format: {target}")

    # Minimal snippets
    if target == "alphafold3":
        frag = {"ligand": {"id": "G", "ccdCodes": components}}
        return json.dumps(frag, indent=2)
    if target == "af3-server":
        return glycan.strip()
    if target == "boltz":
        lines = [
            f"# components: [{', '.join(components)}]",
            "ligand:",
            "  id: G",
            f"  ccd: {components[0]}",
        ]
        return "\n".join(lines) + "\n"
    if target == "chai":
        # Produce a minimal FASTA glycan entry with full nested string
        graph = parse_af3_server_glycan(glycan)
        spec = render_chai_glycan(graph)
        return f">glycan|G\n{spec}\n"
    if target == "protenix":
        frag = {"ligand": {"ligand": "CCD_" + "_".join(components), "count": 1}}
        return json.dumps(frag)
    raise ValueError(f"Unsupported format: {target}")

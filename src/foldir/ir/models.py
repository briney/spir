from __future__ import annotations

from enum import Enum
from typing import Dict, List, Literal, Optional, Union

from pydantic import BaseModel, Field, model_validator


class PolymerType(str, Enum):
    protein = "protein"
    dna = "dna"
    rna = "rna"


class LigandReprType(str, Enum):
    ccd = "ccd"
    smiles = "smiles"
    file = "file"


class Modification(BaseModel):
    position: int = Field(ge=1)
    ccd: str


class PolymerChain(BaseModel):
    id: str
    type: PolymerType
    sequence: str
    modifications: List[Modification] = Field(default_factory=list)

    @model_validator(mode="after")
    def _validate_mod_positions(self) -> "PolymerChain":
        for m in self.modifications:
            if m.position > len(self.sequence):
                raise ValueError(
                    f"Modification at {m.position} exceeds sequence length {len(self.sequence)}"
                )
        return self


class Ligand(BaseModel):
    id: str
    repr_type: LigandReprType
    ccd_codes: List[str] = Field(default_factory=list)
    smiles: Optional[str] = None
    file_path: Optional[str] = None


class Ion(BaseModel):
    id: str
    ccd: str


AtomName = Union[str, int]


class AtomRef(BaseModel):
    """
    Generic atom address in IR:
    - entity_id: chain/molecule ID (AF3/Boltz style)
    - position: residue index for polymers, or component index for multi-CCD ligands, or 1 for single ligands
    - atom: atom name (CCD) or atom index (for SMILES/FILE contexts)
    """

    entity_id: str
    position: int = Field(ge=1)
    atom: AtomName


class CovalentBond(BaseModel):
    a: AtomRef
    b: AtomRef


class ConstraintType(str, Enum):
    contact = "contact"
    pocket = "pocket"


class ContactConstraint(BaseModel):
    type: Literal["contact"] = "contact"
    token1: AtomRef
    token2: AtomRef
    max_distance_angstrom: float = Field(gt=0)


class PocketConstraint(BaseModel):
    type: Literal["pocket"] = "pocket"
    binder_entity_id: str
    contacts: List[AtomRef]
    max_distance_angstrom: float = Field(gt=0)


Constraint = Union[ContactConstraint, PocketConstraint]


class GlycanNode(BaseModel):
    node_id: str
    ccd: str


class GlycanEdge(BaseModel):
    parent: str
    child: str
    parent_atom: Optional[str] = None
    child_atom: Optional[str] = "C1"


class GlycanAttachment(BaseModel):
    polymer_id: str
    polymer_residue_index: int = Field(ge=1)
    polymer_atom: Optional[str] = None
    root_node: str
    root_atom: Optional[str] = "C1"


class Glycan(BaseModel):
    glycan_id: str
    nodes: List[GlycanNode]
    edges: List[GlycanEdge]
    attachments: List[GlycanAttachment] = Field(default_factory=list)


class JobIR(BaseModel):
    name: str
    seeds: List[int] = Field(default_factory=list)

    polymers: List[PolymerChain] = Field(default_factory=list)
    ligands: List[Ligand] = Field(default_factory=list)
    ions: List[Ion] = Field(default_factory=list)

    glycans: List[Glycan] = Field(default_factory=list)
    covalent_bonds: List[CovalentBond] = Field(default_factory=list)
    constraints: List[Constraint] = Field(default_factory=list)


class DocumentIR(BaseModel):
    """Some formats are list-of-jobs."""

    jobs: List[JobIR]

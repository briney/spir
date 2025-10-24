from __future__ import annotations

from enum import Enum
from typing import List, Literal, Optional

from pydantic import BaseModel, Field, model_validator


class EntityType(str, Enum):
    protein = "protein"
    rna = "rna"
    dna = "dna"
    ligand = "ligand"
    ion = "ion"


class Modification(BaseModel):
    ccd: str = Field(..., description="CCD code, e.g., HY3, 6OG")
    position: int = Field(..., ge=1, description="1-based residue/base position")


class AtomRef(BaseModel):
    chain_id: str = Field(..., description="Chain identifier, e.g., 'A'")
    copy_index: Optional[int] = Field(
        default=None, ge=1, description="1-based copy index when applicable"
    )
    residue_index: Optional[int] = Field(
        default=None, ge=1, description="1-based residue index for polymers"
    )
    component_index: Optional[int] = Field(
        default=None, ge=1, description="1-based component index for multi-CCD ligands"
    )
    atom_name: Optional[str] = Field(
        default=None, description="CCD atom name when available"
    )
    atom_index: Optional[int] = Field(
        default=None, ge=0, description="Atom index for SMILES/FILE ligands (0-based)"
    )

    @model_validator(mode="after")
    def _validate_address(self) -> "AtomRef":
        # Require at least one of (residue_index, component_index) for addressable tokens
        if self.residue_index is None and self.component_index is None:
            # For ions or unpositioned ligands, allow bare chain-level reference only if atom specified
            if self.atom_name is None and self.atom_index is None:
                raise ValueError(
                    "AtomRef must include residue_index or component_index, or specify an atom by name/index"
                )
        # Disallow both atom_name and atom_index simultaneously
        if self.atom_name is not None and self.atom_index is not None:
            raise ValueError("Specify either atom_name or atom_index, not both")
        return self


class Bond(BaseModel):
    atom1: AtomRef
    atom2: AtomRef


class PolymerChain(BaseModel):
    type: Literal["protein", "rna", "dna"]
    ids: List[str] = Field(..., min_length=1, description="Chain IDs for copies")
    sequence: str
    modifications: List[Modification] = Field(default_factory=list)
    msa_path: Optional[str] = None
    cyclic: Optional[bool] = False

    @model_validator(mode="after")
    def _validate_sequence(self) -> "PolymerChain":
        if not self.sequence or not isinstance(self.sequence, str):
            raise ValueError("sequence must be a non-empty string")
        return self


class Ligand(BaseModel):
    ids: List[str] = Field(..., min_length=1)
    ccd_codes: Optional[List[str]] = None
    smiles: Optional[str] = None
    # Optional: preserve AF3 Server glycan residues string for round-tripping
    af3_residues: Optional[str] = None

    @model_validator(mode="after")
    def _validate_ligand(self) -> "Ligand":
        if (self.ccd_codes is None) == (self.smiles is None):
            raise ValueError("Ligand requires exactly one of ccd_codes or smiles")
        if self.ccd_codes is not None and len(self.ccd_codes) == 0:
            raise ValueError("ccd_codes cannot be empty if provided")
        return self


class Ion(BaseModel):
    ids: List[str] = Field(..., min_length=1)
    code: str = Field(..., description="Ion CCD code, e.g., MG")


class ComplexInput(BaseModel):
    name: Optional[str] = None
    seeds: List[int] = Field(default_factory=list)
    proteins: List[PolymerChain] = Field(default_factory=list)
    rnas: List[PolymerChain] = Field(default_factory=list)
    dnas: List[PolymerChain] = Field(default_factory=list)
    ligands: List[Ligand] = Field(default_factory=list)
    ions: List[Ion] = Field(default_factory=list)
    bonds: List[Bond] = Field(default_factory=list)
    user_ccd: Optional[str] = None

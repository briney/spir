## Design goals

`FoldIR` should give you:

1. A **single, explicit Intermediate Representation (IR)** that can losslessly represent the union of:

   * polymers (protein/DNA/RNA), modifications, (optionally) MSAs/templates
   * ligands/ions (CCD, SMILES, file-backed)
   * **covalent bonds**
   * optional “soft” constraints (contact/pocket) where supported
   * **glycans as a graph** (because formats disagree sharply here)
2. A set of **dialect adapters**: `parse(target_format) -> IR` and `render(IR) -> target_format`.
3. A packaging layout that is “pip-installable” and modern (`pyproject.toml`, `src/` layout, typed code, tests).

Key format facts that drive the IR design:

* **AF3 Server**: glycans are compact **tree strings** in `proteinChain.glycans[]` with only a protein residue position; **no linkage atom control**.  
* **AF3 (non-Server)**: glycans are **multi-CCD ligands + explicit `bondedAtomPairs`**; converter in AF3 does **not** translate AF3 Server glycans.  
* **Chai-1**: glycan connectivity is encoded inline as `NAG(4-1 NAG...)` and protein↔glycan attachment is a **covalent** entry in a restraints CSV.   
* **Boltz-2**: no glycan type; represent rings as ligands and connect with `constraints: - bond:`.  
* **Protenix**: no `glycans` field; glycans are ligands (including multi-CCD via `CCD_...` concatenation) plus explicit `covalent_bonds` with entity indices.  

---

## IR overview

### High-level data flow

```
                ┌─────────────────────────┐
Input files  ──►│ dialect.parse(...)      │
                └──────────┬──────────────┘
                           ▼
                  ┌────────────────┐
                  │ Raw IR         │  (may have missing IDs, ambiguous glycans,
                  └───────┬────────┘   CCD_ prefixes, etc.)
                          ▼
                  ┌────────────────┐
                  │ normalize(IR)  │  (canonical IDs, CCD code normalization,
                  └───────┬────────┘   glycan graph building/annotation)
                          ▼
                  ┌────────────────┐
                  │ target.render  │  (apply target rules/limits)
                  └───────┬────────┘
                          ▼
                     Output files
```

### Core IR decisions

1. **Everything is explicit and 1-based where formats are 1-based**
   Most of these inputs use 1-based indexing for residues (AF3, AF3 Server, Boltz, Protenix, Chai). Keep IR positions **1-based** and only convert if a target truly differs.

2. **Glycans are represented as a graph of monosaccharide nodes**
   This is the only way to translate between:

   * AF3 Server’s tree string (no linkage atoms)
   * Chai’s inline linkage notation (`4-1`)
   * AF3/Boltz/Protenix explicit atom bonds

3. **Covalent bonds use a shared `AtomRef` model**
   Targets disagree on addressing (chain IDs vs entity indices vs component indices). IR should stay stable and resolvable.

---

## Suggested package layout (`src/` + `pyproject.toml`)

```
FoldIR/
  pyproject.toml
  README.md
  src/foldir/
    __init__.py
    cli.py
    convert.py

    ir/
      __init__.py
      models.py
      normalize.py
      ids.py
      glycans/
        __init__.py
        model.py
        parse_af3_server.py
        parse_chai.py
        render_af3_server.py
        render_chai.py
        resolve_linkages.py

    dialects/
      __init__.py
      base.py
      alphafold3.py
      alphafold3_server.py
      boltz2.py
      chai1.py
      protenix.py

    io/
      __init__.py
      json.py
      yaml.py
      fasta.py
      csv.py

  tests/
    test_glycans_roundtrip.py
    test_convert_smoke.py
```

### `pyproject.toml` skeleton (Hatch)

```toml
[build-system]
requires = ["hatchling>=1.24"]
build-backend = "hatchling.build"

[project]
name = "FoldIR"
version = "0.1.0"
description = "Intermediate representation and converters for protein folding model inputs"
readme = "README.md"
requires-python = ">=3.10"
license = { text = "MIT" }
dependencies = [
  "pydantic>=2.7",
  "PyYAML>=6.0",
  "typer>=0.12",
]

[project.optional-dependencies]
dev = ["pytest>=8", "ruff>=0.6", "mypy>=1.10"]
# Optional: for richer FASTA parsing if you want it
bio = ["biopython>=1.83"]

[project.scripts]
foldir = "foldir.cli:app"

[tool.hatch.build.targets.wheel]
packages = ["src/foldir"]

[tool.ruff]
line-length = 100

[tool.mypy]
python_version = "3.10"
warn_unused_ignores = true
disallow_untyped_defs = true
```

---

## Core IR models (Pydantic)

Put this in `src/foldir/ir/models.py`.

```python
from __future__ import annotations

from enum import Enum
from typing import Literal, Optional, Union, List, Dict
from pydantic import BaseModel, Field, model_validator


class PolymerType(str, Enum):
    protein = "protein"
    dna = "dna"
    rna = "rna"


class LigandReprType(str, Enum):
    ccd = "ccd"       # CCD code(s)
    smiles = "smiles" # SMILES string
    file = "file"     # structure file path


class Modification(BaseModel):
    position: int = Field(ge=1)
    ccd: str  # normalized, no leading "CCD_"


class PolymerChain(BaseModel):
    id: str
    type: PolymerType
    sequence: str
    modifications: List[Modification] = Field(default_factory=list)

    @model_validator(mode="after")
    def _validate_mod_positions(self):
        for m in self.modifications:
            if m.position > len(self.sequence):
                raise ValueError(f"Modification at {m.position} exceeds sequence length {len(self.sequence)}")
        return self


class Ligand(BaseModel):
    id: str
    repr_type: LigandReprType
    # If repr_type == ccd, allow multiple codes for multi-CCD ligands (AF3 non-server)
    ccd_codes: List[str] = Field(default_factory=list)
    smiles: Optional[str] = None
    file_path: Optional[str] = None


class Ion(BaseModel):
    id: str
    ccd: str  # e.g., "MG", "NA"


AtomName = Union[str, int]  # str for CCD atom names, int for atom-index addressing (Protenix SMILES/FILE)


class AtomRef(BaseModel):
    """
    Generic atom address in IR:
    - entity_id: chain/molecule ID (AF3/Boltz style)
    - position: residue index for polymers, or "component index" for multi-CCD ligands, or 1 for single ligands
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
    token1: AtomRef  # allow residue-level by using atom="CA" or leaving conventions to renderer
    token2: AtomRef
    max_distance_angstrom: float = Field(gt=0)


class PocketConstraint(BaseModel):
    type: Literal["pocket"] = "pocket"
    binder_entity_id: str
    contacts: List[AtomRef]
    max_distance_angstrom: float = Field(gt=0)


Constraint = Union[ContactConstraint, PocketConstraint]


# ------------------ Glycans (first-class) ------------------

class GlycanNode(BaseModel):
    node_id: str
    ccd: str  # e.g. NAG, MAN, BMA ...


class GlycanEdge(BaseModel):
    parent: str
    child: str
    parent_atom: Optional[str] = None  # e.g. "O4", "O6"
    child_atom: Optional[str] = "C1"   # usually "C1"


class GlycanAttachment(BaseModel):
    polymer_id: str
    polymer_residue_index: int = Field(ge=1)
    polymer_atom: Optional[str] = None  # e.g. ND2/OG/OG1/N (Chai often uses "N")
    root_node: str
    root_atom: Optional[str] = "C1"


class Glycan(BaseModel):
    glycan_id: str
    nodes: List[GlycanNode]
    edges: List[GlycanEdge]
    attachments: List[GlycanAttachment] = Field(default_factory=list)


class JobIR(BaseModel):
    name: str
    # Only meaningful for AlphaFold inputs; others may ignore
    seeds: List[int] = Field(default_factory=list)

    polymers: List[PolymerChain] = Field(default_factory=list)
    ligands: List[Ligand] = Field(default_factory=list)
    ions: List[Ion] = Field(default_factory=list)

    glycans: List[Glycan] = Field(default_factory=list)
    covalent_bonds: List[CovalentBond] = Field(default_factory=list)
    constraints: List[Constraint] = Field(default_factory=list)


class DocumentIR(BaseModel):
    """Some formats are list-of-jobs (AF3 Server, Protenix)."""
    jobs: List[JobIR]
```

---

## Normalization layer

Normalization belongs in `src/foldir/ir/normalize.py`. It should:

1. Normalize CCD codes:

   * Strip leading `CCD_` when present (AF3 Server ligands/mods; Protenix mods/examples)  
2. Ensure entity IDs are unique and stable:

   * AF3 non-server requires explicit IDs; AF3 Server doesn’t provide them (only `count`).  
3. Resolve glycan ambiguity where possible:

   * AF3 Server glycan strings provide topology but **no linkage atoms**. 
   * Decide whether to store missing linkage atoms as `None` and defer to renderers (recommended).

Example normalization helper:

```python
def normalize_ccd(code: str) -> str:
    return code[4:] if code.startswith("CCD_") else code
```

---

## Glycan parsers and renderers

This is the “hard part” and deserves dedicated modules.

### 1) AF3 Server glycan string parser

AF3 Server glycan is a rooted tree string like `NAG(NAG)(BMA)` with up to 2 children per node; no linkage positions.  

`src/foldir/ir/glycans/parse_af3_server.py`:

```python
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple
import re

from foldir.ir.models import Glycan, GlycanNode, GlycanEdge

CCD_RE = re.compile(r"[A-Za-z0-9]{3}")

class ParseError(ValueError):
    pass


def parse_af3_server_glycan_string(glycan_id: str, s: str) -> Glycan:
    """
    Parses AF3 Server compact glycan tree string into a Glycan graph.
    Linkage atoms are unknown in this format -> edges have parent_atom=None.
    """
    i = 0
    nodes: List[GlycanNode] = []
    edges: List[GlycanEdge] = []
    node_counter = 0

    def parse_node() -> str:
        nonlocal i, node_counter
        m = CCD_RE.match(s, i)
        if not m:
            raise ParseError(f"Expected CCD code at offset {i}: ...{s[i:i+10]!r}")
        ccd = m.group(0)
        i = m.end()

        node_id = f"{glycan_id}.n{node_counter}"
        node_counter += 1
        nodes.append(GlycanNode(node_id=node_id, ccd=ccd))

        # zero, one, or two children encoded as (...)(...)
        while i < len(s) and s[i] == "(":
            i += 1
            child_id = parse_node()
            if i >= len(s) or s[i] != ")":
                raise ParseError(f"Missing ')' at offset {i}")
            i += 1
            edges.append(GlycanEdge(parent=node_id, child=child_id, parent_atom=None, child_atom="C1"))

        return node_id

    root = parse_node()
    if i != len(s):
        raise ParseError(f"Trailing junk at offset {i}: {s[i:]!r}")

    return Glycan(glycan_id=glycan_id, nodes=nodes, edges=edges, attachments=[])
```

### 2) Chai glycan string parser

Chai encodes linkage positions inline: `NAG(4-1 NAG(4-1 BMA(3-1 MAN)(6-1 MAN)))`.  

`src/foldir/ir/glycans/parse_chai.py`:

```python
from __future__ import annotations
import re
from typing import List, Optional

from foldir.ir.models import Glycan, GlycanNode, GlycanEdge

CCD_RE = re.compile(r"[A-Za-z0-9]{3}")
INT_RE = re.compile(r"\d+")

class ParseError(ValueError):
    pass


def parse_chai_glycan_string(glycan_id: str, s: str) -> Glycan:
    i = 0
    nodes: List[GlycanNode] = []
    edges: List[GlycanEdge] = []
    node_counter = 0

    def read_ccd() -> str:
        nonlocal i
        m = CCD_RE.match(s, i)
        if not m:
            raise ParseError(f"Expected CCD at {i}")
        i = m.end()
        return m.group(0)

    def read_int() -> int:
        nonlocal i
        m = INT_RE.match(s, i)
        if not m:
            raise ParseError(f"Expected int at {i}")
        i = m.end()
        return int(m.group(0))

    def skip_ws():
        nonlocal i
        while i < len(s) and s[i].isspace():
            i += 1

    def parse_node() -> str:
        nonlocal i, node_counter
        ccd = read_ccd()
        node_id = f"{glycan_id}.n{node_counter}"
        node_counter += 1
        nodes.append(GlycanNode(node_id=node_id, ccd=ccd))

        while i < len(s) and s[i] == "(":
            i += 1
            parent_pos = read_int()
            if i >= len(s) or s[i] != "-":
                raise ParseError(f"Expected '-' after parent_pos at {i}")
            i += 1
            child_pos = read_int()
            skip_ws()
            child_id = parse_node()
            if i >= len(s) or s[i] != ")":
                raise ParseError(f"Missing ')' at {i}")
            i += 1

            edges.append(
                GlycanEdge(
                    parent=node_id,
                    child=child_id,
                    parent_atom=f"O{parent_pos}",
                    child_atom=f"C{child_pos}",
                )
            )

        return node_id

    root = parse_node()
    if i != len(s):
        raise ParseError(f"Trailing junk: {s[i:]!r}")

    return Glycan(glycan_id=glycan_id, nodes=nodes, edges=edges, attachments=[])
```

### 3) Render glycan to AF3 Server string

Renderer must validate AF3 Server limitations: <=2 children per node, <=8 residues, plus root restrictions by residue type if you choose to enforce them. 

`src/foldir/ir/glycans/render_af3_server.py`:

```python
from __future__ import annotations
from collections import defaultdict
from foldir.ir.models import Glycan

class RenderError(ValueError):
    pass


def render_af3_server_glycan_string(g: Glycan, root_node_id: str) -> str:
    nodes = {n.node_id: n for n in g.nodes}
    children = defaultdict(list)
    parents = {}
    for e in g.edges:
        children[e.parent].append(e.child)
        if e.child in parents:
            raise RenderError("Not a tree (child has multiple parents).")
        parents[e.child] = e.parent

    if len(g.nodes) > 8:
        raise RenderError("AF3 Server supports up to 8 glycan residues.")
    for p, kids in children.items():
        if len(kids) > 2:
            raise RenderError("AF3 Server glycan nodes may have at most 2 children.")

    def rec(nid: str) -> str:
        s = nodes[nid].ccd
        for kid in children.get(nid, []):
            s += f"({rec(kid)})"
        return s

    return rec(root_node_id)
```

### 4) Render glycan to Chai string

Requires linkage positions. If edges lack `parent_atom/child_atom`, either:

* error, or
* fill defaults via a linkage resolver (recommended; see below). 

`src/foldir/ir/glycans/render_chai.py`:

```python
from __future__ import annotations
from collections import defaultdict
from foldir.ir.models import Glycan

class RenderError(ValueError):
    pass

def _atom_to_pos(atom: str) -> int:
    # "O4" -> 4, "C1" -> 1
    if len(atom) < 2 or not atom[1:].isdigit():
        raise RenderError(f"Cannot convert atom to pos: {atom}")
    return int(atom[1:])

def render_chai_glycan_string(g: Glycan, root_node_id: str) -> str:
    nodes = {n.node_id: n for n in g.nodes}
    children = defaultdict(list)  # parent -> list[(edge, child)]
    parents = {}
    for e in g.edges:
        if e.parent_atom is None or e.child_atom is None:
            raise RenderError("Chai requires explicit linkage positions/atoms.")
        children[e.parent].append((e, e.child))
        if e.child in parents:
            raise RenderError("Not a tree (child has multiple parents).")
        parents[e.child] = e.parent

    def rec(nid: str) -> str:
        s = nodes[nid].ccd
        for edge, kid in children.get(nid, []):
            ppos = _atom_to_pos(edge.parent_atom)
            cpos = _atom_to_pos(edge.child_atom)
            s += f"({ppos}-{cpos} {rec(kid)})"
        return s

    return rec(root_node_id)
```

---

## Linkage resolution strategy (critical for AF3 Server → explicit formats)

AF3 Server doesn’t let you specify glycosidic linkage atoms. So if you convert an AF3 Server glycan to AF3 non-server / Chai / Boltz / Protenix, you must decide linkages.

Do this via a pluggable resolver:

`src/foldir/ir/glycans/resolve_linkages.py`:

```python
from __future__ import annotations
from typing import Protocol, Tuple
from foldir.ir.models import Glycan, GlycanEdge

class LinkageResolver(Protocol):
    def resolve(self, parent_ccd: str, child_ccd: str) -> tuple[str, str]:
        """
        Return (parent_atom, child_atom), e.g. ("O4","C1").
        """

class DefaultSugarLinkageResolver:
    """
    Conservative defaults:
    - use child_atom = C1 (common)
    - use parent_atom = O4 (common) unless otherwise configured
    """
    def __init__(self, default_parent_atom: str = "O4", default_child_atom: str = "C1"):
        self.default_parent_atom = default_parent_atom
        self.default_child_atom = default_child_atom

    def resolve(self, parent_ccd: str, child_ccd: str) -> tuple[str, str]:
        return (self.default_parent_atom, self.default_child_atom)


def fill_missing_linkages(g: Glycan, resolver: LinkageResolver) -> Glycan:
    nodes = {n.node_id: n.ccd for n in g.nodes}
    new_edges = []
    for e in g.edges:
        if e.parent_atom is None or e.child_atom is None:
            p_atom, c_atom = resolver.resolve(nodes[e.parent], nodes[e.child])
            new_edges.append(e.model_copy(update={"parent_atom": p_atom, "child_atom": c_atom}))
        else:
            new_edges.append(e)
    return g.model_copy(update={"edges": new_edges})
```

In practice, you’ll likely want:

* a per-project default (`O4->C1`)
* optional overrides per edge / per monosaccharide pair
* attachment atom defaults per residue type when missing (Asn vs Ser/Thr)

This mirrors the fact that AF3 Server chooses linkages heuristically. 

---

## Dialect adapters

Define a common interface in `src/foldir/dialects/base.py`:

```python
from __future__ import annotations
from typing import Protocol, Any
from foldir.ir.models import DocumentIR

class Dialect(Protocol):
    name: str
    def parse(self, path: str) -> DocumentIR: ...
    def render(self, doc: DocumentIR, out_path: str) -> None: ...
```

Then implement one module per target.

### AlphaFold3 (non-server) dialect notes

* Single-job JSON object with `dialect: alphafold3` and `version` (example shows 4). 
* Requires `modelSeeds` non-empty. 
* Glycans: represent as ligand `ccdCodes:[...]` plus explicit `bondedAtomPairs`.  

Renderer core (sketch) in `src/foldir/dialects/alphafold3.py`:

```python
import json
from foldir.ir.models import DocumentIR

def render_alphafold3(doc: DocumentIR, out_path: str) -> None:
    if len(doc.jobs) != 1:
        raise ValueError("AlphaFold3 (non-server) expects exactly one job per JSON file.")
    job = doc.jobs[0]

    if not job.seeds:
        # deterministic default: choose 1, or hash(job.name) mod 2**32, etc.
        job = job.model_copy(update={"seeds": [1]})

    sequences = []
    for p in job.polymers:
        if p.type.value == "protein":
            sequences.append({"protein": {"id": p.id, "sequence": p.sequence,
                                          "modifications": [{"ptmType": m.ccd, "ptmPosition": m.position}
                                                            for m in p.modifications] or None}})
        elif p.type.value == "dna":
            sequences.append({"dna": {"id": p.id, "sequence": p.sequence}})
        elif p.type.value == "rna":
            sequences.append({"rna": {"id": p.id, "sequence": p.sequence}})

    for lig in job.ligands:
        if lig.repr_type.value == "ccd":
            sequences.append({"ligand": {"id": lig.id, "ccdCodes": lig.ccd_codes}})
        elif lig.repr_type.value == "smiles":
            sequences.append({"ligand": {"id": lig.id, "smiles": lig.smiles}})
        else:
            raise ValueError("AF3 non-server JSON does not accept FILE ligands directly (use userCCD/CCD or SMILES).")

    for ion in job.ions:
        # AF3 treats ions as ligands with CCD codes, e.g. ["MG"] :contentReference[oaicite:27]{index=27}
        sequences.append({"ligand": {"id": ion.id, "ccdCodes": [ion.ccd]}})

    bonded_pairs = [
        [[b.a.entity_id, b.a.position, b.a.atom], [b.b.entity_id, b.b.position, b.b.atom]]
        for b in job.covalent_bonds
    ]

    payload = {
        "name": job.name,
        "modelSeeds": job.seeds,
        "sequences": sequences,
        "bondedAtomPairs": bonded_pairs or None,
        "dialect": "alphafold3",
        "version": 4,
    }

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
```

Glycans should be injected into `job.ligands` + `job.covalent_bonds` during normalization/render prep (see “glycan expansion” below).

### AlphaFold3 Server dialect notes

* Top-level JSON is a **list** of job dicts. 
* Per-protein glycans are `proteinChain.glycans[] = {residues: "...", position: int}`. 
* Ligands/ions are limited in the server UI, but you may still want to parse/emit them. 

### Boltz-2 dialect notes

* YAML schema with `sequences`, `constraints`, `templates`, `properties`. 
* Bond constraints use `atom1: [CHAIN_ID, RES_IDX, ATOM_NAME]`. 
* Glycans are represented as ligands (CCD) + bond constraints (recommended for translation). 

### Chai-1 dialect notes

* You need at least:

  * FASTA entries for polymers and glycans
  * restraints CSV containing covalent bonds and/or contact/pocket
* Glycan intra-connectivity is encoded in glycan FASTA record string; protein↔glycan via a covalent row in restraints.  

### Protenix dialect notes

* Top-level JSON is **list of jobs** like AF3 Server. 
* Glycans: use ligand strings like `CCD_NAG_BMA_BGC` + `covalent_bonds` referencing **entity indices**.  

---

## Glycan expansion into target-specific representations

A practical pattern is:

* Keep `JobIR.glycans` as the canonical representation
* When rendering to explicit-bond formats (AF3 non-server, Boltz, Protenix, Chai), **expand** glycans into:

  * target ligand entities (one or many)
  * target bond objects (bondedAtomPairs / constraints.bond / covalent_bonds / restraints row)
* When parsing those explicit-bond formats, optionally **attempt to detect glycans** and populate `JobIR.glycans` (heuristics), but always preserve the explicit bonds too.

### Expansion to AF3 non-server (multi-CCD ligand + bondedAtomPairs)

AF3 wants glycan rings inside one ligand entity as `ccdCodes: [...]` and internal edges are bonds where `position` indexes the component in `ccdCodes` (1..N).  

Pseudo-implementation:

```python
from __future__ import annotations
from collections import defaultdict, deque

from foldir.ir.models import Glycan, Ligand, LigandReprType, CovalentBond, AtomRef

def _topo_tree_order(root: str, children_map: dict[str, list[str]]) -> list[str]:
    # deterministic DFS pre-order
    out = []
    stack = [root]
    while stack:
        n = stack.pop()
        out.append(n)
        kids = children_map.get(n, [])
        # reverse for stable left-to-right
        for k in reversed(kids):
            stack.append(k)
    return out

def expand_glycan_to_af3(g: Glycan, ligand_id: str, root_node: str) -> tuple[Ligand, list[CovalentBond], dict[str, int]]:
    nodes = {n.node_id: n.ccd for n in g.nodes}
    children = defaultdict(list)
    for e in g.edges:
        children[e.parent].append(e.child)

    order = _topo_tree_order(root_node, children)
    idx_map = {node_id: i + 1 for i, node_id in enumerate(order)}  # 1-based component index
    ligand = Ligand(id=ligand_id, repr_type=LigandReprType.ccd, ccd_codes=[nodes[n] for n in order])

    bonds: list[CovalentBond] = []
    for e in g.edges:
        # edge direction is parent -> child, positions from idx_map
        bonds.append(
            CovalentBond(
                a=AtomRef(entity_id=ligand_id, position=idx_map[e.parent], atom=e.parent_atom or "O4"),
                b=AtomRef(entity_id=ligand_id, position=idx_map[e.child],  atom=e.child_atom or "C1"),
            )
        )
    for att in g.attachments:
        bonds.append(
            CovalentBond(
                a=AtomRef(entity_id=att.polymer_id, position=att.polymer_residue_index, atom=att.polymer_atom or "ND2"),
                b=AtomRef(entity_id=ligand_id, position=idx_map[att.root_node], atom=att.root_atom or "C1"),
            )
        )

    return ligand, bonds, idx_map
```

**Important**: defaulting `ND2`/`O4` may be wrong for a specific system; treat these as configurable policies, not silent magic. The IR should keep `None` until a user or resolver fills it.

### Expansion to Boltz-2 (many ligand entries + constraints.bond)

Boltz bonds refer to `[chain_id, residue_index, atom_name]`, and for ligands residue index is typically `1`.  

Recommended pattern: one ligand entity per monosaccharide node; then each edge becomes a constraint bond.

### Expansion to Protenix (multi-CCD ligand string + covalent_bonds)

Protenix uses entity indices in `sequences` (1-based), and for multi-CCD ligands `position` is the component number.  

So the mapping step is:

* determine the ligand entity index in the output sequences list
* map glycan node → component index in the concatenated list

### Expansion to Chai-1 (glycan FASTA + covalent restraint row)

* Glycan string requires linkage positions (`4-1` etc). 
* Protein↔glycan attachment is a covalent row in restraints CSV.
* Restraints want residue-letter redundancy like `N436@ND2` (Chai checks consistency). 

Your Chai renderer should therefore:

1. Order FASTA entries deterministically (because Chai assigns chain letters in order). 
2. Create glycan FASTA record using `render_chai_glycan_string(...)`.
3. Write restraints CSV row:

   * `chainA` = protein chain letter
   * `res_idxA` = `<AA><pos>@<atom>`
   * `chainB` = glycan chain letter
   * `res_idxB` = `@C1` (or `@<root_atom>`)

---

## “Converter” orchestration API

In `src/foldir/convert.py`:

```python
from __future__ import annotations
from dataclasses import dataclass
from foldir.ir.models import DocumentIR
from foldir.ir.normalize import normalize_document
from foldir.dialects import get_dialect

@dataclass(frozen=True)
class ConvertOptions:
    # policies for ambiguous cases
    default_seed: int = 1
    default_glycan_parent_atom: str = "O4"
    default_glycan_child_atom: str = "C1"
    default_asn_atom: str = "ND2"
    default_ser_atom: str = "OG"
    default_thr_atom: str = "OG1"

def convert(in_path: str, in_dialect: str, out_path: str, out_dialect: str, opts: ConvertOptions) -> None:
    src = get_dialect(in_dialect)
    dst = get_dialect(out_dialect)

    doc = src.parse(in_path)
    doc = normalize_document(doc, opts=opts)
    dst.render(doc, out_path)
```

---

## CLI (Typer)

`src/foldir/cli.py`:

```python
import typer
from foldir.convert import convert, ConvertOptions

app = typer.Typer(no_args_is_help=True)

@app.command()
def convert_cmd(
    in_path: str = typer.Argument(...),
    in_dialect: str = typer.Option(..., "--from"),
    out_path: str = typer.Argument(...),
    out_dialect: str = typer.Option(..., "--to"),
):
    opts = ConvertOptions()
    convert(in_path, in_dialect, out_path, out_dialect, opts)

if __name__ == "__main__":
    app()
```

Usage:

```bash
foldir convert --from alphafoldserver input.json --to alphafold3 output
foldir convert --from alphafold3 input.json --to boltz2 output
foldir convert --from chai1 ./case/ --to protenix output
```

(For Chai you’ll typically accept a directory containing `*.fasta` + `*.csv`.)

---

## Format-specific mapping checklist (what your adapters must implement)

### AlphaFold3 Server → IR

* Parse list-of-jobs top-level. 
* Expand each `proteinChain` count into distinct polymer IDs or store copies.
* Parse `proteinChain.glycans[].residues` into `Glycan` via `parse_af3_server_glycan_string`.
* Store attachment: `position` is residue index; linkage atoms unknown. 

### AlphaFold3 (non-Server) → IR

* Read single JSON object with `sequences`.
* Read `bondedAtomPairs` and store as `CovalentBond`s. 
* Detect glycan candidates:

  * multi-CCD ligand with internal O*→C1 bonds
  * protein↔ligand bond to a “C1” atom
* Convert those into `JobIR.glycans[]` where possible; otherwise leave as generic ligands + bonds.

### Boltz-2 → IR

* Parse YAML sequences list.
* Map `constraints: - bond:` to covalent bonds. 
* Keep pocket/contact constraints if present. 

### Chai-1 → IR

* Parse FASTA entries; identify `>glycan|...` and parse glycan strings via `parse_chai_glycan_string`. 
* Parse restraints CSV:

  * extract covalent rows and convert to attachments/bonds
  * pocket/contact rows map to IR constraints. 

### Protenix → IR

* Parse list-of-jobs JSON. 
* Parse ligands:

  * `CCD_...` codes (possibly multi-CCD concatenation)
  * SMILES
  * `FILE_...` paths 
* Parse `covalent_bonds`:

  * map entity indices → IR entity IDs using output ordering
  * positions are residue index or component index 

---

## Testing strategy (minimal but high-value)

1. **Golden glycan roundtrips**

   * AF3 Server string → IR → AF3 Server string (exact match)
   * Chai string → IR → Chai string (exact match)
2. **Cross-model glycan conversions**

   * AF3 Server → IR → AF3 non-server should generate:

     * one multi-CCD ligand
     * expected number of bonds (edges + attachment)
   * Chai → IR → Boltz should generate:

     * one ligand per node
     * constraints.bond count matches edges + attachment
3. **Schema validation tests**

   * Ensure rendered outputs conform to required top-level shapes:

     * AF3 Server and Protenix are lists-of-jobs  
     * AF3 non-server is single-job object with `dialect/version` 
     * Boltz is YAML with `sequences` and optional `constraints` 

---

## Practical caveats you should encode as explicit errors/warnings

1. **AF3 Server glycan ambiguity**

   * If converting AF3 Server glycans into explicit-bond formats, linkage atoms are not provided. Your conversion must either:

     * require user-supplied linkage policy, or
     * apply documented defaults and emit a warning.

2. **AF3 Server limitations**

   * Max 8 residues and max 2 children per node in glycan strings. 
     When converting *to* AF3 Server, reject glycans that exceed these constraints.

3. **Atom naming mismatches**

   * Chai/AF3/Boltz/Protenix require correct atom names for bonds (CCD atom naming).
     The IR should preserve whatever the user specified and avoid “guessing” silently.

4. **Seeds**

   * AF3 non-server requires a seed list, whereas AF3 Server allows empty seeds (random).
     Your normalizer should deterministically fill missing seeds when rendering to AF3 non-server.

---

## Summary of the IR fields needed specifically for glycans

Your `Glycan` model (nodes/edges/attachment) covers the “union” of glycan requirements described across all five tools. 

* Nodes: `ccd` per monosaccharide instance
* Edges:

  * for Chai/Boltz/AF3/Protenix: `parent_atom` and `child_atom`
  * for AF3 Server: may be unknown (`None`)
* Attachment:

  * polymer chain id + residue index
  * polymer atom name (may be unknown)
  * root atom name (usually `C1`)

That structure is sufficient to render:

* AF3 Server tree strings (ignore linkage atoms)
* Chai inline linkage strings + covalent restraint row
* AF3 non-server multi-CCD ligand + `bondedAtomPairs`
* Boltz ligands + `constraints: bond`
* Protenix multi-CCD ligand string + `covalent_bonds`

---

If you implement `DocumentIR` + normalization + five dialect adapters as above, you’ll have a maintainable `FoldIR` codebase that can translate glycans and other covalent connectivity across AlphaFold3 Server, AlphaFold3 (non-server), Chai-1, Boltz-2, and Protenix, while keeping ambiguity and format limitations explicit rather than hidden.

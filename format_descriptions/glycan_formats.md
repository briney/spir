# Glycan formats by model (AlphaFold3 Server, AlphaFold3, Chai-1, Boltz-2, Protenix)

This document summarizes how each model expects **glycans** to be represented in its input format, since glycan handling is one of the biggest sources of friction when translating between formats.

Where relevant, examples below reuse the Man3 / Man6 strings shown in `AF3-server_glycan-formatting.md`.

## Quick comparison (what you need to model in an IR)

| Model | Where glycans are represented | “Glycan object” representation | Explicit glycosidic linkage atoms/positions? | Protein↔glycan attachment |
|---|---|---|---|---|
| **AlphaFold3 Server** | `proteinChain.glycans[]` | Compact **tree string** using CCD residue codes (e.g. `NAG(NAG(MAN...))`) | **No** (Server chooses bond atoms heuristically) | `proteinChain.glycans[].position` (1-based residue index; attachment chemistry not user-specifiable) |
| **AlphaFold3 (non-Server)** | `sequences[]` + `bondedAtomPairs[]` | **Multi-CCD ligand** (`ccdCodes: [...]`) + explicit bond list | **Yes** (`bondedAtomPairs`) | Explicit bond in `bondedAtomPairs` between protein atom and glycan atom |
| **Chai-1** | FASTA glycan record + restraints CSV | CCD codes + **inline linkage positions** (e.g. `NAG(4-1 NAG)`) | **Yes** (e.g. `4-1` means O4→C1) | `connection_type=covalent` line in restraints CSV |
| **Boltz-2** | YAML `sequences[]` + `constraints[]` | No special glycan type; treat sugars as ligands (CCD or SMILES) | **Yes** (via `constraints: - bond:`) | Explicit `bond` constraint between protein atom and ligand atom |
| **Protenix** | `sequences[]` + `covalent_bonds[]` | No `glycans` field; use ligands (CCD/SMILES/FILE). Multi-CCD ligands supported via `CCD_...` concatenation | **Yes** (`covalent_bonds`) | Explicit bond(s) in `covalent_bonds` |

## Suggested intermediate representation (IR) fields for glycan translation

To translate across all five targets, your IR generally needs to capture a **graph**:

- **Nodes**: monosaccharide (or “chemical component”) instances
  - `node_id` (stable ID within glycan)
  - `ccd_code` (e.g. `NAG`, `MAN`, `BMA`, `FUC`, ...)
- **Edges**: covalent linkages between nodes
  - `parent_node_id`, `child_node_id`
  - `parent_atom` (e.g. `O4`, `O6`, `O3`) and `child_atom` (commonly `C1`) when required by the target format
  - If linkage atoms are **unknown/unspecified** (e.g. AlphaFold3 Server), store this explicitly and allow the target to choose defaults.
- **Attachment(s) to polymer(s)**:
  - `protein_entity_id` / `chain_id` / `entity_index` (depending on model)
  - `protein_residue_index` (1-based in every format shown here)
  - `protein_atom_name` (e.g. `ND2`, `OG`, `OG1`, or whatever the target needs)
  - `glycan_root_node_id` and `glycan_root_atom_name` (often `C1`)

The remainder of this document explains how each model encodes (parts of) this information.

---

## AlphaFold3 Server (`dialect: alphafoldserver`)

### Where glycans are specified (Server JSON)

Within each `proteinChain` entity, glycans are declared as an optional list:

- `proteinChain.glycans[]` entries contain:
  - `residues` (**string**): glycan encoded as a compact tree string (see below)
  - `position` (**int**, 1-based): residue index in the protein sequence the glycan is attached to

### Glycan “residues” string format

AlphaFold3 Server uses a compact syntax that describes a **rooted tree** of monosaccharide residues:

- **Residue identifiers** are **3-letter PDB CCD codes** (Chemical Components Dictionary).
  - Stereoisomers use different CCD codes (e.g. mannose can be `MAN` vs `BMA`).
- **Tree structure** is encoded with parentheses:
  - `NAG` is a single residue.
  - `NAG(BMA)` means “NAG has one child BMA”.
  - `NAG(FUC)(NAG)` means “NAG has two children: FUC and NAG”.
- **Branching limits**:
  - Each residue can have **0–2 children**.
  - Up to **8 total glycan residues** are supported.
- **Allowed root residues depend on the protein residue type**:
  - Attached to **N (Asn)**: `BGC`, `BMA`, `GLC`, `MAN`, `NAG`
  - Attached to **S (Ser)** / **T (Thr)**: `BGC`, `BMA`, `FUC`, `GLC`, `MAN`, `NAG`
- **No explicit linkage atoms/positions**:
  - You cannot specify which atoms form the glycosidic bond.
  - The Server chooses bond atoms heuristically based on frequent occurrences in the PDB.
  - The Server expects glycan-glycan connections to be chemically valid (e.g. `GLC(NAG)(MAN)` is cited as invalid).

### Examples (glycan strings)

These are directly usable in `proteinChain.glycans[].residues`:

- **Single residue**: `NAG`
- **One child**: `NAG(BMA)`
- **Linear chain**: `NAG(BMA(BGC))`
- **Two children**: `NAG(FUC)(NAG)`
- **Man3 (linear example)**: `NAG(NAG(MAN(MAN(MAN))))`
- **Man6 (branched example)**: `NAG(NAG(MAN(MAN(MAN)(MAN(MAN)(MAN))))))`

### Example JSON snippet (protein glycosylation)

```json
{
  "name": "Example AF3 Server glycosylation",
  "modelSeeds": [],
  "sequences": [
    {
      "proteinChain": {
        "sequence": "PREACHINGS",
        "count": 1,
        "glycans": [
          {
            "residues": "NAG(NAG(MAN(MAN(MAN))))",
            "position": 8
          }
        ]
      }
    }
  ],
  "dialect": "alphafoldserver",
  "version": 1
}
```

---

## AlphaFold3 (non-Server) (`dialect: alphafold3`, `version: 4`)

### Key idea: “multi-CCD ligand + explicit bonds”

AlphaFold3 (non-Server) does not have an `alphafoldserver`-style glycan string. Instead:

- A glycan is represented as a **ligand entity** whose `ccdCodes` list contains **one CCD code per glycan component**.
- Covalent connectivity (protein↔glycan and glycan↔glycan) is expressed with **explicit bonds** in `bondedAtomPairs`.

### Ligand definition (multi-CCD)

- `ligand.id` is a single uppercase letter (or list of letters for multiple copies).
- `ligand.ccdCodes` is a list of CCD codes. If the list has length > 1, you typically also provide `bondedAtomPairs` to connect components.

### Bond representation (`bondedAtomPairs`)

Each bonded atom is addressed by three fields:

- **Entity ID** (`str`): the `id` of the entity (e.g. `"A"` for a protein chain, `"G"` for a ligand).
- **Residue index within the entity** (`int`, 1-based):
  - For proteins/RNA/DNA: the residue index in the sequence.
  - For multi-CCD ligands: the **component index** within `ccdCodes` (1..N).
  - For single-residue ligands: always `1`.
- **Atom name** (`str`): the CCD atom name (e.g. `C1`, `O4`, `O6`, `ND2`, ...).

Important constraints called out in the AF3 format description:

- Bonds are covalent; other bond types are not supported.
- SMILES ligands do **not** support `bondedAtomPairs` (no stable atom naming); use CCD codes or a user-provided CCD if you need bonded connectivity.

### Example: Man3 as a multi-CCD ligand

This example attaches the glycan to residue 8 of the protein chain (sequence `"PREACHINGS"` has `N` at position 8). Bond atoms shown are illustrative and should be validated against the CCD atom naming for the specific components you use.

```json
{
  "name": "Example AF3 Man3 (multi-CCD ligand)",
  "modelSeeds": [1],
  "sequences": [
    {
      "protein": {
        "id": "A",
        "sequence": "PREACHINGS"
      }
    },
    {
      "ligand": {
        "id": "G",
        "ccdCodes": ["NAG", "NAG", "MAN", "MAN", "MAN"],
        "description": "Man3-like linear chain: NAG-NAG-MAN-MAN-MAN"
      }
    }
  ],
  "bondedAtomPairs": [
    [["A", 8, "ND2"], ["G", 1, "C1"]],
    [["G", 1, "O4"], ["G", 2, "C1"]],
    [["G", 2, "O4"], ["G", 3, "C1"]],
    [["G", 3, "O6"], ["G", 4, "C1"]],
    [["G", 4, "O6"], ["G", 5, "C1"]]
  ],
  "dialect": "alphafold3",
  "version": 4
}
```

### Example: Man6 (branched) as a multi-CCD ligand

This matches the **topology** of the AlphaFold3 Server Man6 example string `NAG(NAG(MAN(MAN(MAN)(MAN(MAN)(MAN))))))` (2 NAG + 6 MAN = 8 components), but uses explicit bonds to encode branching.

```json
{
  "name": "Example AF3 Man6 (multi-CCD ligand, branched)",
  "modelSeeds": [1],
  "sequences": [
    {
      "protein": {
        "id": "A",
        "sequence": "PREACHINGS"
      }
    },
    {
      "ligand": {
        "id": "G",
        "ccdCodes": ["NAG", "NAG", "MAN", "MAN", "MAN", "MAN", "MAN", "MAN"],
        "description": "Man6-like branched tree (see bondedAtomPairs)"
      }
    }
  ],
  "bondedAtomPairs": [
    [["A", 8, "ND2"], ["G", 1, "C1"]],

    [["G", 1, "O4"], ["G", 2, "C1"]],
    [["G", 2, "O4"], ["G", 3, "C1"]],
    [["G", 3, "O6"], ["G", 4, "C1"]],

    [["G", 4, "O3"], ["G", 5, "C1"]],
    [["G", 4, "O6"], ["G", 6, "C1"]],

    [["G", 6, "O3"], ["G", 7, "C1"]],
    [["G", 6, "O6"], ["G", 8, "C1"]]
  ],
  "dialect": "alphafold3",
  "version": 4
}
```

---

## Chai-1

### Where glycans are specified (FASTA + restraints)

Chai-1 represents glycans in two places:

1. **FASTA input** includes a glycan record, e.g.:
   - `>glycan|example-name`
   - followed by a single-line glycan string (see below)
2. **Restraints CSV** specifies covalent attachments (including protein↔glycan) using a row with `connection_type=covalent`.

### Glycan string format (inline linkage positions)

Chai-1 uses CCD codes and an inline bond syntax:

- Start with the **root** CCD code (e.g. `NAG`).
- Attach a child by writing:
  - `(PARENT_POS-CHILD_POS CHILD_CCD)`
  - Example: `NAG(4-1 NAG)` means bond **O4** (parent) → **C1** (child).
- Build outward left-to-right; parentheses always attach to the residue immediately preceding them.
- Branch by adding multiple parenthetical attachments:
  - `BMA(3-1 MAN)(6-1 MAN)` attaches two children to `BMA`.

### Protein↔glycan attachment (restraints CSV)

In a restraints `.csv`, add a covalent bond row. In the example below:

- `chainA=A` is the first FASTA entry (protein)
- `chainB=B` is the second FASTA entry (glycan)
- `res_idxA` includes both a residue+index and an atom name (e.g. `N8@ND2`)
- `res_idxB` can refer to an atom in the root glycan as `@C1`

```text
chainA|res_idxA|chainB|res_idxB|connection_type|confidence|min_distance_angstrom|max_distance_angstrom|comment|restraint_id
A|N8@ND2|B|@C1|covalent|1.0|0.0|0.0|protein-glycan|bond1
```

### Examples (Chai-1 glycan FASTA strings)

- **Single residue**:
  - `NAG`
- **Two residues**:
  - `NAG(4-1 NAG)`
- **Man3 (linear example, topology matching the AF3 Server Man3 string)**:
  - `NAG(4-1 NAG(4-1 MAN(6-1 MAN(6-1 MAN))))`
- **Man6 (branched example, topology matching the AF3 Server Man6 string)**:
  - `NAG(4-1 NAG(4-1 MAN(6-1 MAN(3-1 MAN)(6-1 MAN(3-1 MAN)(6-1 MAN)))))`

Chai-1 also notes that sugar rings include hydroxyl groups that leave when bonds form; it attempts to drop glycan leaving atoms automatically for glycan rings.

---

## Boltz-2

### Key idea: “no glycan type; use ligands + bond constraints”

Boltz-2’s input schema has:

- `sequences`: proteins/DNA/RNA and `ligand` entries
- `constraints`: optional covalent `bond` constraints between two atoms

There is no special “glycan string” field. A glycan can be represented by:

- **Option A (single ligand)**: one `ligand` with a full-glycan **SMILES** (most compact, but SMILES authoring is non-trivial).
- **Option B (recommended for translation)**: one `ligand` per monosaccharide ring (CCD code), plus `constraints: - bond:` entries that connect rings and attach the root ring to the protein.

Boltz-2 states that `bond` constraints are currently supported only for **CCD ligands** and **canonical residues**. Atom names should be verified against the component’s CCD mmCIF.

### Example: Man3 (linear) as multiple CCD ligands + bonds

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: PREACHINGS
      msa: empty
  - ligand:
      id: G1
      ccd: NAG
  - ligand:
      id: G2
      ccd: NAG
  - ligand:
      id: G3
      ccd: MAN
  - ligand:
      id: G4
      ccd: MAN
  - ligand:
      id: G5
      ccd: MAN
constraints:
  - bond:
      atom1: [A, 8, ND2]
      atom2: [G1, 1, C1]
  - bond:
      atom1: [G1, 1, O4]
      atom2: [G2, 1, C1]
  - bond:
      atom1: [G2, 1, O4]
      atom2: [G3, 1, C1]
  - bond:
      atom1: [G3, 1, O6]
      atom2: [G4, 1, C1]
  - bond:
      atom1: [G4, 1, O6]
      atom2: [G5, 1, C1]
```

### Example: Man6 (branched) as multiple CCD ligands + bonds

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: PREACHINGS
      msa: empty
  - ligand: { id: G1, ccd: NAG }
  - ligand: { id: G2, ccd: NAG }
  - ligand: { id: G3, ccd: MAN }
  - ligand: { id: G4, ccd: MAN }
  - ligand: { id: G5, ccd: MAN }
  - ligand: { id: G6, ccd: MAN }
  - ligand: { id: G7, ccd: MAN }
  - ligand: { id: G8, ccd: MAN }
constraints:
  - bond: { atom1: [A, 8, ND2], atom2: [G1, 1, C1] }

  - bond: { atom1: [G1, 1, O4], atom2: [G2, 1, C1] }
  - bond: { atom1: [G2, 1, O4], atom2: [G3, 1, C1] }
  - bond: { atom1: [G3, 1, O6], atom2: [G4, 1, C1] }

  - bond: { atom1: [G4, 1, O3], atom2: [G5, 1, C1] }
  - bond: { atom1: [G4, 1, O6], atom2: [G6, 1, C1] }

  - bond: { atom1: [G6, 1, O3], atom2: [G7, 1, C1] }
  - bond: { atom1: [G6, 1, O6], atom2: [G8, 1, C1] }
```

---

## Protenix

### Key idea: “no glycans field; use ligands + covalent_bonds”

Protenix explicitly states:

- There is **no supported `glycans` field**.
- Glycans can be represented as:
  - multiple ligands with defined bonding, or
  - a single ligand described by a full-molecule SMILES / structure file
- A ligand can be specified as a CCD code prefixed with `CCD_`.
  - For glycans, Protenix supports a **multi-CCD ligand string** that concatenates multiple CCD codes, e.g. `CCD_NAG_BMA_BGC`.

### Covalent bonds (`covalent_bonds`)

`covalent_bonds` entries identify atoms by:

- `entity1` / `entity2`: 1-based entity index in the `sequences` list
- `copy1` / `copy2`: 1-based copy index (optional)
- `position1` / `position2`:
  - For polymers: residue index in the sequence (1-based)
  - For multi-CCD ligands: **component index** within the concatenated CCD list (1-based)
  - For single-CCD / SMILES / FILE ligands: always `1`
- `atom1` / `atom2`:
  - For polymers or CCD-defined ligands: **CCD atom names**
  - For SMILES/FILE ligands: atoms can be specified by **atom index** (0-based), per the Protenix doc

### Example: Man3 as a multi-CCD ligand + bonds

```json
[
  {
    "name": "Example Protenix Man3",
    "sequences": [
      {
        "proteinChain": {
          "sequence": "PREACHINGS",
          "count": 1
        }
      },
      {
        "ligand": {
          "ligand": "CCD_NAG_NAG_MAN_MAN_MAN",
          "count": 1
        }
      }
    ],
    "covalent_bonds": [
      { "entity1": "1", "copy1": 1, "position1": "8", "atom1": "ND2", "entity2": "2", "copy2": 1, "position2": "1", "atom2": "C1" },
      { "entity1": "2", "copy1": 1, "position1": "1", "atom1": "O4", "entity2": "2", "copy2": 1, "position2": "2", "atom2": "C1" },
      { "entity1": "2", "copy1": 1, "position1": "2", "atom1": "O4", "entity2": "2", "copy2": 1, "position2": "3", "atom2": "C1" },
      { "entity1": "2", "copy1": 1, "position1": "3", "atom1": "O6", "entity2": "2", "copy2": 1, "position2": "4", "atom2": "C1" },
      { "entity1": "2", "copy1": 1, "position1": "4", "atom1": "O6", "entity2": "2", "copy2": 1, "position2": "5", "atom2": "C1" }
    ]
  }
]
```

### Example: Man6 (branched) as a multi-CCD ligand + bonds

```json
[
  {
    "name": "Example Protenix Man6 (branched)",
    "sequences": [
      {
        "proteinChain": {
          "sequence": "PREACHINGS",
          "count": 1
        }
      },
      {
        "ligand": {
          "ligand": "CCD_NAG_NAG_MAN_MAN_MAN_MAN_MAN_MAN",
          "count": 1
        }
      }
    ],
    "covalent_bonds": [
      { "entity1": "1", "copy1": 1, "position1": "8", "atom1": "ND2", "entity2": "2", "copy2": 1, "position2": "1", "atom2": "C1" },

      { "entity1": "2", "copy1": 1, "position1": "1", "atom1": "O4", "entity2": "2", "copy2": 1, "position2": "2", "atom2": "C1" },
      { "entity1": "2", "copy1": 1, "position1": "2", "atom1": "O4", "entity2": "2", "copy2": 1, "position2": "3", "atom2": "C1" },
      { "entity1": "2", "copy1": 1, "position1": "3", "atom1": "O6", "entity2": "2", "copy2": 1, "position2": "4", "atom2": "C1" },

      { "entity1": "2", "copy1": 1, "position1": "4", "atom1": "O3", "entity2": "2", "copy2": 1, "position2": "5", "atom2": "C1" },
      { "entity1": "2", "copy1": 1, "position1": "4", "atom1": "O6", "entity2": "2", "copy2": 1, "position2": "6", "atom2": "C1" },

      { "entity1": "2", "copy1": 1, "position1": "6", "atom1": "O3", "entity2": "2", "copy2": 1, "position2": "7", "atom2": "C1" },
      { "entity1": "2", "copy1": 1, "position1": "6", "atom1": "O6", "entity2": "2", "copy2": 1, "position2": "8", "atom2": "C1" }
    ]
  }
]
```

---

## Practical notes for formats requiring explicit atoms (AF3 / Chai / Boltz / Protenix)

- **Atom names must match the CCD atom naming** for the component (e.g. `C1`, `O4`, `O3`, `O6` are common for sugars, and `ND2`/`OG`/`OG1` for typical glycosylation sites).
- If your linkage atoms differ from the examples above, adjust the bonds accordingly.
- AlphaFold3 Server is the only format here that **does not let you specify** linkage atoms/positions; it will infer them.

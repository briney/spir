Below is a description of how each model expects glycans to be **specified in inputs**, and contains two examples of AlphaFold‑Server strings converted into the proper format for each model:

* **AlphaFold 3 (local model; not the web server)**
* **Boltz‑2**
* **Chai‑1**
* **Protenix**

> **Important note on link chemistry.**
> AF3 **Server** accepts a parenthesis string like `NAG(NAG(MAN(MAN(MAN))))` and infers both **topology** and **glycosidic link positions** automatically. All four models below **do not** infer link positions from that bare string. You must express the bonds explicitly (either as atom pairs or via a glycan grammar that includes the numeric linkage, e.g., `4‑1`, `6‑1`). Where a bond is needed, the convention is: **parent O[n] ↔ child C1** (so a “4‑1” linkage means *O4 of the parent sugar bonded to C1 of the child*). This is how AF3 (local) and Chai‑1 document bond semantics, and it’s also how to wire ligands in Boltz‑2 and Protenix. ([GitHub][1])

---

## How to read the two examples

* **Example A (linear):** `NAG(NAG(MAN(MAN(MAN))))`
  *One chain of five monosaccharides. For a concrete illustration of formatting, I’ve used canonical N‑glycan‑core linkages:*
  **NAG(4‑1)NAG(4‑1)MAN(2‑1)MAN(2‑1)MAN** → bonds at **O4→C1** for NAG→NAG and **O2→C1** for MAN→MAN.

* **Example B (branched):** `NAG(NAG(MAN(MAN(MAN)(MAN(MAN)(MAN)))))`
  *This is the Man3 core (two GlcNAc + a branched mannose arm). I’ve used the standard N‑glycan branching:* the middle MAN has branches **3‑1** and **6‑1** to two MAN residues, and each branch has a terminal **2‑1** MAN. In atom terms: from the branching MAN, **O3→C1** and **O6→C1** start the two arms; each arm extends by **O2→C1** to its terminal mannose.*

If you need a different linkage pattern, just change the O[n] on the **parent** side (the **child** is always **C1**).

---

## 1) AlphaFold 3 (local model; not the server)

**What the model expects.** You pass a JSON with entities (proteins, ligands, etc.) and an optional `bondedAtomPairs` list. Each atom is addressed as `[entity_id, residue_index (1‑based inside that ligand), atom_name]`. AF3 (local) explicitly shows a glycan bond within a multi‑CCD ligand as `["J", 1, "O6"]` ↔ `["J", 2, "C1"]`. Glycans are defined as ligands with **multiple CCD codes** plus explicit bonds (including bonds to the protein if glycosylated). ([GitHub][1])

### Example A (linear) — AF3 JSON fragment

```json
{
  "dialect": "alphafold3",
  "version": 4,
  "modelSeeds": [1],
  "name": "Linear glycan example",
  "sequences": [
    { "ligand": { "id": "G", "ccdCodes": ["NAG", "NAG", "MAN", "MAN", "MAN"] } }
  ],
  "bondedAtomPairs": [
    [["G", 1, "O4"], ["G", 2, "C1"]],  // NAG(4-1)NAG
    [["G", 2, "O4"], ["G", 3, "C1"]],  // NAG(4-1)MAN
    [["G", 3, "O2"], ["G", 4, "C1"]],  // MAN(2-1)MAN
    [["G", 4, "O2"], ["G", 5, "C1"]]   // MAN(2-1)MAN
  ]
}
```

### Example B (branched) — AF3 JSON fragment

```json
{
  "dialect": "alphafold3",
  "version": 4,
  "modelSeeds": [1],
  "name": "Branched glycan example",
  "sequences": [
    { "ligand": { "id": "G", "ccdCodes": ["NAG","NAG","MAN","MAN","MAN","MAN","MAN","MAN"] } }
  ],
  "bondedAtomPairs": [
    [["G", 1, "O4"], ["G", 2, "C1"]],  // NAG(4-1)NAG
    [["G", 2, "O4"], ["G", 3, "C1"]],  // NAG(4-1)MAN  (root MAN at index 3)

    [["G", 3, "O2"], ["G", 4, "C1"]],  // trunk: MAN(2-1)MAN
    [["G", 4, "O2"], ["G", 5, "C1"]],  // trunk: MAN(2-1)MAN

    [["G", 3, "O3"], ["G", 6, "C1"]],  // branch 1 from MAN@3: (3-1)
    [["G", 6, "O2"], ["G", 7, "C1"]],  // extend branch 1: (2-1)

    [["G", 3, "O6"], ["G", 8, "C1"]]   // branch 2 from MAN@3: (6-1)
  ]
}
```

*Why this is “correct” for AF3 (local).* AF3’s official `input.md` specifies multi‑CCD ligands and `bondedAtomPairs` with **O[n]** ↔ **C1** atom names and 1‑based residue indices inside the ligand; it even highlights glycan bonds using that addressing scheme. ([GitHub][1])

---

## 2) Boltz‑2

**What the model expects.** Boltz‑2 takes a YAML “schema” and lets you define **CCD ligands** and **bond constraints**. Bonds are specified as:

```yaml
constraints:
  - bond:
      atom1: [CHAIN_ID, RES_IDX, ATOM_NAME]
      atom2: [CHAIN_ID, RES_IDX, ATOM_NAME]
```

This is supported for CCD ligands and canonical residues; atom names come from the CCD. (SMILES ligands cannot be used for covalent bonds in Boltz.) ([GitHub][2])

> **Tip:** For a glycan made of several monosaccharides, define each sugar as a **ligand with `ccd: NAG/MAN/BMA`** (copy count = 1) and wire them with `constraints: bond:` pairs. (You can also group multiple CCD residues under one Boltz ligand entity if you’ve customized the loader, but the documented path is one CCD per ligand with bonds between them.)

### Example A (linear) — Boltz‑2 YAML fragment

```yaml
version: 1
sequences:
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
  - bond: {atom1: [G1, 1, O4], atom2: [G2, 1, C1]}  # NAG(4-1)NAG
  - bond: {atom1: [G2, 1, O4], atom2: [G3, 1, C1]}  # NAG(4-1)MAN
  - bond: {atom1: [G3, 1, O2], atom2: [G4, 1, C1]}  # MAN(2-1)MAN
  - bond: {atom1: [G4, 1, O2], atom2: [G5, 1, C1]}  # MAN(2-1)MAN
```

### Example B (branched) — Boltz‑2 YAML fragment

```yaml
version: 1
sequences:
  - ligand: {id: G1, ccd: NAG}
  - ligand: {id: G2, ccd: NAG}
  - ligand: {id: G3, ccd: MAN}  # branching MAN
  - ligand: {id: G4, ccd: MAN}  # trunk 1
  - ligand: {id: G5, ccd: MAN}  # trunk 2
  - ligand: {id: B1, ccd: MAN}  # branch 1 base
  - ligand: {id: B2, ccd: MAN}  # branch 1 tip
  - ligand: {id: C1, ccd: MAN}  # branch 2 base

constraints:
  - bond: {atom1: [G1, 1, O4], atom2: [G2, 1, C1]}  # NAG->NAG (4-1)
  - bond: {atom1: [G2, 1, O4], atom2: [G3, 1, C1]}  # NAG->MAN (4-1)

  - bond: {atom1: [G3, 1, O2], atom2: [G4, 1, C1]}  # trunk: (2-1)
  - bond: {atom1: [G4, 1, O2], atom2: [G5, 1, C1]}  # trunk: (2-1)

  - bond: {atom1: [G3, 1, O3], atom2: [B1, 1, C1]}  # branch 1 from MAN@G3: (3-1)
  - bond: {atom1: [B1, 1, O2], atom2: [B2, 1, C1]}  # extend branch 1: (2-1)

  - bond: {atom1: [G3, 1, O6], atom2: [C1, 1, C1]}  # branch 2 from MAN@G3: (6-1)
```

*Why this is “correct” for Boltz‑2.* The repo’s `prediction.md` and issues document **bond constraints** in YAML with `atom1/atom2` lists of `[CHAIN, RES_IDX, ATOM_NAME]`, and clarify that covalent bonds are supported for **CCD** ligands (not SMILES). ([GitHub][2])

---

## 3) Chai‑1

**What the model expects.** Chai‑1 supports **covalent bonds and glycans** via:

* a **FASTA line** with a **glycan string that includes numeric bonds** like `NAG(4-1 NAG)`; parentheses denote attachment to the immediately preceding sugar; `4‑1` means O4(parent)–C1(child).
* a separate **CSV “restraints” file** to connect the **root glycan residue** to a specific protein residue/atom if you’re modeling a glycosylated protein. The examples and docs show this format and explicitly call out the `4‑1` bond semantics. ([GitHub][3])

> **Only the glycan strings are shown below** (the restraint row that bonds the root sugar to the protein is separate). If you’re folding the glycan alone, you can omit the restraint.

### Example A (linear) — Chai‑1 glycan FASTA line

```
>glycan|linear
NAG(4-1 NAG(4-1 MAN(2-1 MAN(2-1 MAN))))
```

### Example B (branched) — Chai‑1 glycan FASTA line

```
>glycan|branched
NAG(4-1 NAG(4-1 MAN(2-1 MAN(2-1 MAN))(3-1 MAN(2-1 MAN))(6-1 MAN)))
```

This follows the grammar documented by Chai‑1 (explicit `n‑1` bonds; building “outward” left‑to‑right), and mirrors the corrected branched Man3 example given in their issue tracker. ([GitHub][3])

---

## 4) Protenix

**What the model expects.** Protenix uses a JSON **similar to AlphaFold‑Server** but adds:

* Support for multi‑CCD ligands by concatenating CCDs (e.g., `"CCD_NAG_BMA_BGC"`), and
* Explicit `covalent_bonds` specifying the two atoms to connect.
  For bonds, each side is addressed by `[entity_number, copy_index, position, atom_name]`, where **`position` is the serial index of the CCD within that ligand** (or residue index for polymers). You can also split sugars into separate ligands and connect them. ([Hugging Face][4])

### Example A (linear) — Protenix JSON fragment

```json
{
  "name": "Linear glycan example",
  "sequences": [
    { "ligand": { "ligand": "CCD_NAG_NAG_MAN_MAN_MAN", "count": 1 } }
  ],
  "covalent_bonds": [
    { "entity1": "1", "position1": "1", "atom1": "O4", "entity2": "1", "position2": "2", "atom2": "C1" },  // NAG(4-1)NAG
    { "entity1": "1", "position1": "2", "atom1": "O4", "entity2": "1", "position2": "3", "atom2": "C1" },  // NAG(4-1)MAN
    { "entity1": "1", "position1": "3", "atom1": "O2", "entity2": "1", "position2": "4", "atom2": "C1" },  // MAN(2-1)MAN
    { "entity1": "1", "position1": "4", "atom1": "O2", "entity2": "1", "position2": "5", "atom2": "C1" }   // MAN(2-1)MAN
  ]
}
```

### Example B (branched) — Protenix JSON fragment

Two ways are supported. I show (i) **one multi‑CCD ligand** with positions, or (ii) **split branches as separate ligands**.

**(i) One multi‑CCD ligand**

```json
{
  "name": "Branched glycan example (single ligand)",
  "sequences": [
    { "ligand": { "ligand": "CCD_NAG_NAG_MAN_MAN_MAN_MAN_MAN_MAN", "count": 1 } }
  ],
  "covalent_bonds": [
    { "entity1": "1", "position1": "1", "atom1": "O4", "entity2": "1", "position2": "2", "atom2": "C1" },  // NAG->NAG
    { "entity1": "1", "position1": "2", "atom1": "O4", "entity2": "1", "position2": "3", "atom2": "C1" },  // NAG->MAN (root @pos3)

    { "entity1": "1", "position1": "3", "atom1": "O2", "entity2": "1", "position2": "4", "atom2": "C1" },  // trunk: (2-1)
    { "entity1": "1", "position1": "4", "atom1": "O2", "entity2": "1", "position2": "5", "atom2": "C1" },  // trunk: (2-1)

    { "entity1": "1", "position1": "3", "atom1": "O3", "entity2": "1", "position2": "6", "atom2": "C1" },  // branch 1 from MAN@3: (3-1)
    { "entity1": "1", "position1": "6", "atom1": "O2", "entity2": "1", "position2": "7", "atom2": "C1" },  // extend branch 1: (2-1)

    { "entity1": "1", "position1": "3", "atom1": "O6", "entity2": "1", "position2": "8", "atom2": "C1" }   // branch 2 from MAN@3: (6-1)
  ]
}
```

**(ii) Split into multiple ligands** (easier to read/edit)

```json
{
  "name": "Branched glycan example (split ligands)",
  "sequences": [
    { "ligand": { "ligand": "CCD_NAG", "count": 1 } },  // entity 1, pos 1
    { "ligand": { "ligand": "CCD_NAG", "count": 1 } },  // entity 2, pos 1
    { "ligand": { "ligand": "CCD_MAN_MAN_MAN", "count": 1 } },        // entity 3, pos 1..3 (root/trunk)
    { "ligand": { "ligand": "CCD_MAN_MAN", "count": 1 } },            // entity 4, pos 1..2 (branch 1)
    { "ligand": { "ligand": "CCD_MAN", "count": 1 } }                 // entity 5, pos 1 (branch 2 leaf)
  ],
  "covalent_bonds": [
    { "entity1": "1", "position1": "1", "atom1": "O4", "entity2": "2", "position2": "1", "atom2": "C1" },  // NAG->NAG
    { "entity1": "2", "position1": "1", "atom1": "O4", "entity2": "3", "position2": "1", "atom2": "C1" },  // NAG->MAN (root)

    { "entity1": "3", "position1": "1", "atom1": "O2", "entity2": "3", "position2": "2", "atom2": "C1" },  // trunk (2-1)
    { "entity1": "3", "position1": "2", "atom1": "O2", "entity2": "3", "position2": "3", "atom2": "C1" },  // trunk (2-1)

    { "entity1": "3", "position1": "1", "atom1": "O3", "entity2": "4", "position2": "1", "atom2": "C1" },  // branch 1 start (3-1)
    { "entity1": "4", "position1": "1", "atom1": "O2", "entity2": "4", "position2": "2", "atom2": "C1" },  // branch 1 extend (2-1)

    { "entity1": "3", "position1": "1", "atom1": "O6", "entity2": "5", "position2": "1", "atom2": "C1" }   // branch 2 start (6-1)
  ]
}
```

*Why this is “correct” for Protenix.* The official format guide explains multi‑CCD ligands like `"CCD_NAG_BMA_BGC"` and the `covalent_bonds` addressing scheme with `entity`, `position`, and `atom` fields for bonding ligands to ligands or polymers. ([Hugging Face][4])

---

## Quick atom‑mapping cheat sheet (for NAG, MAN, BMA)

* **Child side** is always **`C1`**.
* **Parent side** atom equals **`O2` / `O3` / `O4` / `O6`** for **2‑1 / 3‑1 / 4‑1 / 6‑1** linkages respectively.
  This is the convention used in AF3 (local) `bondedAtomPairs` examples and Chai‑1’s glycan grammar (`n‑1` means parent **O[n]** ↔ child **C1**). ([GitHub][1])

---

## Sources

* **AlphaFold 3 (local)** input format, multi‑CCD ligands and `bondedAtomPairs` (with glycan example inside a ligand). ([GitHub][1])
* **Boltz (Boltz‑2) docs**: YAML schema and **bond constraints** (`atom1`/`atom2` = `[CHAIN_ID, RES_IDX, ATOM_NAME]`); CCD ligands only for covalent bonds. ([GitHub][2])
* **Chai‑1 glycans**: FASTA glycan grammar using **`(n‑1`)** bond notation and restraint CSV; branched Man3 example with corrected ordering. ([GitHub][3])
* **Protenix** JSON: multi‑CCD ligands via `"CCD_..."` concatenation; explicit `covalent_bonds` with `entity`/`position`/`atom` addressing for ligand–ligand or ligand–polymer bonds. ([Hugging Face][4])

---

### If you want these wired to a protein site

Add one more **bond** from the **root sugar’s C1** to the **protein atom** (e.g., **Asn ND2** for an N‑linked glycan) using each model’s bond syntax:

* AF3: `["A", <asn_index>, "ND2"]` ↔ `["G", <root_pos>, "C1"]`. ([GitHub][1])
* Boltz‑2: `atom1: [A, <asn_index>, ND2]`, `atom2: [G1, 1, C1]`. ([GitHub][2])
* Chai‑1: add one CSV restraint row with `connection_type=covalent` between `A,N<idx>@ND2` and the glycan root `@C1` as shown in their covalent‑bonds example. ([GitHub][5])
* Protenix: a `covalent_bonds` entry where `entity1` is the protein entity, `atom1="ND2"`, and `entity2` is the glycan ligand at `position` of the root sugar with `atom2="C1"`. ([Hugging Face][4])

If you’d like, I can also generate ready‑to‑run full input files (including a protein and the exact Asn site) for one of your sequences.

[1]: https://raw.githubusercontent.com/google-deepmind/alphafold3/main/docs/input.md "raw.githubusercontent.com"
[2]: https://raw.githubusercontent.com/jwohlwend/boltz/main/docs/prediction.md?utm_source=chatgpt.com "raw.githubusercontent.com"
[3]: https://github.com/chaidiscovery/chai-lab/blob/main/examples/covalent_bonds/README.md?utm_source=chatgpt.com "chai-lab/examples/covalent_bonds/README.md at main - GitHub"
[4]: https://huggingface.co/datasets/Yimingbear/protenix/blob/main/docs/infer_json_format.md "docs/infer_json_format.md · Yimingbear/protenix at main"
[5]: https://github.com/chaidiscovery/chai-lab/blob/main/examples/restraints/README.md?utm_source=chatgpt.com "chai-lab/examples/restraints/README.md at main - GitHub"

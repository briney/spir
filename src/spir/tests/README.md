## Test Suite Overview

This directory contains the automated tests for `spir`. Tests are grouped by type to help quickly identify failures during CI and triage issues.

### Unit tests

- **`test_chain_ids_and_detection.py`**
  - **Purpose**: Validates utility behavior for generating spreadsheet-style IDs used for chain labeling.
  - **What it checks**:
    - `spreadsheet_ids(1)` returns `["A"]`
    - `spreadsheet_ids(3)` returns `["A", "B", "C"]`
    - Boundary behavior around wrap‑around to double letters:
      - `spreadsheet_ids(28)[-2:] == ["AA", "AB"]`
      - `spreadsheet_ids(29)[-2:] == ["AB", "AC"]`
  - **Why it matters**: Ensures stable and predictable chain ID generation across single- and double-letter ranges, preventing downstream mismatches in adapters or file outputs that rely on these IDs.
  - **If this fails**: Likely causes include changes to the ID generation algorithm in `spir.utils.ids.spreadsheet_ids` (off‑by‑one, base‑26 logic), or regressions introduced by refactoring.

- **`test_convert_af3_server_glycan.py`**
  - **Purpose**: Verifies glycan conversion from AF3‑Server style strings into all supported formats with deterministic linkage inference and a mock protein anchor.
  - **What it checks**:
    - Minimal mode: AF3 contains full `ccdCodes`; Protenix uses concatenated `CCD_...`; Boltz and Chai show the first CCD as a snippet; AF3‑Server echo is preserved.
    - Full mode (linear and branched examples):
      - AF3: multi‑CCD ligand with complete `bondedAtomPairs` for `O[n]↔C1` intra‑glycan bonds, plus anchor `A:1@ND2 ↔ G:1@C1`.
      - Boltz: one CCD per ligand (`G1..Gn`) and `constraints.bond` entries for all edges and the anchor.
      - Chai: FASTA glycan string using `(n‑1 ...)` grammar and a single covalent restraint anchoring `A,1@ND2` to `G,@C1`.
      - Protenix: single multi‑CCD ligand with `covalent_bonds` addressing by positions, including the anchor.
  - **Why it matters**: Exercises the hardest part of format conversion—glycan topology and bond semantics—ensuring conformance with `format_descriptions/glycan_formats.md`.
  - **If this fails**: Inspect `parse_af3_server_glycan` and `_infer_parent_atom_for_link` in `spir.ir.glycan`, and the per‑format builders in `spir.glycan_convert` (`_graph_to_ci_for_af3/_boltz/_protenix` and `_graph_to_fasta_and_restraints_for_chai`).

### Integration / smoke tests

- **`test_roundtrip_smoke.py`**
  - **Purpose**: Smoke test for AlphaFold3 round‑trip through the new high‑level API (`spir.read()` → `IR.write_alphafold3()`).
  - **What it checks**:
    - Writes a minimal AF3 JSON to a temp file.
    - Calls `spir.read(temp_path)` to build an `IR`.
    - Calls `IR.write_alphafold3(temp_dir)` and validates the output contains `"alphafold3"`.
  - **Why it matters**: Ensures the public API integrates correctly with adapters and default naming behavior.
  - **If this fails**: Check `spir.io.read`, writer dispatch in `spir.ir.ir.IR`, and AF3 adapter serialization.

New high‑level API snippet:

```python
from spir import read

ir = read("example.yaml")
ir.write_alphafold3("out/")
```

## Conventions and maintenance

- **Add tests** adjacent to related code when feasible, and update this README with a one‑line purpose and brief “If this fails” triage note.
- **Name tests** to reflect their role: prefer `*_unit.py` for focused utilities and `*_smoke.py` or `*_e2e.py` for broader flows.
- **Running tests**:
  - All tests: `pytest -q`
  - Match by keyword: `pytest -k af3` or `pytest -k spreadsheet_ids`

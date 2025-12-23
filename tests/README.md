# SPIR Test Suite

This directory contains the test suite for SPIR (Structure Prediction Input Representation). Tests are organized by functionality and can be run using pytest:

```bash
# Run all tests
pytest

# Run tests with verbose output
pytest -v

# Run a specific test file
pytest tests/test_validate.py
```

## Test Categories

### CLI Tests (`test_cli.py`)

Tests for the command-line interface, ensuring the `spir` CLI works correctly.

| Test | Description | Purpose |
|------|-------------|---------|
| `test_convert_help` | Verifies `spir convert --help` displays correct options | Ensures CLI documentation is accessible |
| `test_convert_basic` | Tests basic conversion from AlphaFold Server to AlphaFold3 | Validates the convert command works end-to-end |
| `test_validate_help` | Verifies `spir validate --help` displays correct options | Ensures validation CLI documentation is accessible |
| `test_validate_valid_file` | Tests validation passes for a correct AlphaFold3 file | Confirms valid files are accepted |
| `test_validate_invalid_file` | Tests validation fails for a file missing required fields | Confirms invalid files are rejected with proper exit code |
| `test_validate_nonexistent_file` | Tests validation of a non-existent file | Ensures graceful handling of missing files |
| `test_validate_with_short_option` | Tests `-d` short option for `--dialect` | Validates CLI short options work correctly |
| `test_main_help` | Verifies `spir --help` shows both commands | Ensures top-level help lists all available commands |

---

### Validation Tests (`test_validate.py`)

Tests for the input file validation functionality across all supported dialects.

#### ValidationResult Unit Tests

| Test | Description | Purpose |
|------|-------------|---------|
| `test_empty_result_is_valid` | Empty ValidationResult should be valid | Baseline validation result behavior |
| `test_warnings_dont_invalidate` | Warnings alone should not fail validation | Distinguish warnings from errors |
| `test_errors_invalidate` | Errors should cause validation to fail | Ensure errors properly invalidate |
| `test_merge` | Merging results combines all issues | Support for modular validation |

#### AlphaFold3 Validation

| Test | Description | Purpose |
|------|-------------|---------|
| `test_valid_alphafold3` | Valid AF3 file passes validation | Confirms proper AF3 format is accepted |
| `test_missing_model_seeds` | Missing modelSeeds field is caught | Required field validation |
| `test_empty_model_seeds` | Empty modelSeeds array is rejected | At least one seed is required |
| `test_missing_sequence` | Protein without sequence field fails | Required field validation |
| `test_duplicate_entity_id` | Duplicate entity IDs are detected | Entity ID uniqueness enforcement |
| `test_invalid_bond_entity_ref` | Bonds referencing non-existent entities fail | Referential integrity checking |
| `test_wrong_dialect_warning` | Wrong dialect value triggers warning | Dialect mismatch detection (non-fatal) |
| `test_invalid_json` | Malformed JSON is caught | Parse error handling |

#### Boltz2 Validation

| Test | Description | Purpose |
|------|-------------|---------|
| `test_valid_boltz2` | Valid Boltz2 YAML passes validation | Confirms proper Boltz2 format is accepted |
| `test_missing_sequences` | Missing sequences field is caught | Required field validation |
| `test_invalid_bond_entity` | Bonds referencing non-existent entities fail | Referential integrity checking |

#### AlphaFold Server Validation

| Test | Description | Purpose |
|------|-------------|---------|
| `test_valid_server_format` | Valid AF Server format passes | Confirms proper server format is accepted |
| `test_not_a_list` | Single job (not a list) is rejected | AF Server requires job array |
| `test_empty_sequence` | Empty sequence string is rejected | Non-empty sequence required |

#### Protenix Validation

| Test | Description | Purpose |
|------|-------------|---------|
| `test_valid_protenix` | Valid Protenix format passes | Confirms proper Protenix format is accepted |
| `test_invalid_bond_entity_ref` | Out-of-range entity reference fails | Entity index bounds checking |

#### Chai1 Validation

| Test | Description | Purpose |
|------|-------------|---------|
| `test_valid_chai1` | Valid FASTA with protein and ligand passes | Confirms proper Chai1 format is accepted |
| `test_invalid_entity_type` | Unknown entity type in header fails | Entity type validation |
| `test_empty_sequence` | Empty sequence is rejected | Non-empty sequence required |
| `test_valid_with_restraints` | Valid FASTA + restraints CSV passes | Multi-file input validation |
| `test_invalid_connection_type` | Invalid connection_type in restraints fails | Constraint type validation |

---

### Conversion Tests

#### Smoke Tests (`test_convert_smoke.py`)

Quick sanity checks for basic conversion functionality.

| Test | Description | Purpose |
|------|-------------|---------|
| `test_af3_server_to_af3` | Convert AF Server to AF3 with glycans | Verifies glycans are converted to bondedAtomPairs |

#### Roundtrip Tests (`test_roundtrip.py`)

Tests that data survives conversion through intermediate formats and back.

| Test | Description | Purpose |
|------|-------------|---------|
| `test_roundtrip_af3_server_boltz_af3_server` | AF Server → Boltz → AF Server | Ensures glycan data survives Boltz roundtrip |
| `test_roundtrip_af3_server_chai_af3_server` | AF Server → Chai → AF Server | Ensures glycan data survives Chai roundtrip |

---

### Glycan Tests (`test_glycans_roundtrip.py`)

Tests for glycan string parsing and rendering in different dialects.

| Test | Description | Purpose |
|------|-------------|---------|
| `test_af3_server_roundtrip` | Parse and render AF Server glycan string | Ensures glycan structure survives parse/render cycle |
| `test_chai_roundtrip` | Parse and render Chai glycan string | Ensures branched glycan structure survives parse/render |

---

### MSA Path Tests (`test_msa_path.py`)

Tests for Multiple Sequence Alignment (MSA) path preservation across formats.

| Test | Description | Purpose |
|------|-------------|---------|
| `test_alphafold3_msa_path_roundtrip` | AF3 → Boltz → AF3 preserves MSA path | MSA path survives roundtrip conversion |
| `test_boltz_msa_path_roundtrip` | Boltz → AF3 → Boltz preserves MSA path | MSA path survives roundtrip conversion |
| `test_af3_server_msa_path_to_boltz` | AF Server with `msa_path` converts to Boltz | Non-standard extension support |
| `test_af3_server_msa_path_to_alphafold3` | AF Server with `msa_path` converts to AF3 | Non-standard extension support |
| `test_af3_server_output_never_includes_msa` | AF Server output has no MSA fields | Strict spec compliance for AF Server output |

---

### Chai1 Restraints Tests (`test_chai1_restraints.py`)

Tests for Chai1-specific restraints handling.

| Test | Description | Purpose |
|------|-------------|---------|
| `test_chai1_restraints_explicit_path` | Convert Chai1 with explicit restraints path | Validates `--restraints` option for Chai1 input |

---

## Test Data

The `test_data/` directory contains sample input files used by some tests:

| File | Description |
|------|-------------|
| `fold_dupilumab_with_il_4ra_truncated_with_man3_glycans_job_request.json` | AF Server format with glycans |
| `fold_dupilumab_with_il_4ra_truncated_with_man3_glycans_job_request_alphafold3.json` | AF3 format with glycans |
| `fold_dupilumab_with_il_4ra_truncated_with_man3_glycans_job_request_chai1.fasta` | Chai1 FASTA output |
| `fold_dupilumab_with_il_4ra_truncated_with_man3_glycans_job_request_chai1.constraints.csv` | Chai1 constraints |
| `fold_dupilumab_with_il_4ra_truncated_with_man3_glycans_job_request_to-boltz2.yaml` | Boltz2 YAML output |
| `fold_il_31_with_known_mab_job_request_with-msa.json` | AF Server format with MSA paths |
| `fold_il_31_with_known_mab_job_request_with-msa_alphafold3.json` | AF3 format with MSA paths |

---

## Configuration

Test configuration is managed in `conftest.py`, which sets up the Python path to include the `src/` directory for imports.


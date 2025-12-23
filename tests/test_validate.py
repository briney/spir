"""Tests for SPIR validation functionality."""

import json

import pytest

from spir.validate import ValidationResult, Severity, validate


class TestValidationResult:
    def test_empty_result_is_valid(self):
        result = ValidationResult()
        assert result.is_valid
        assert result.error_count == 0
        assert result.warning_count == 0

    def test_warnings_dont_invalidate(self):
        result = ValidationResult()
        result.add_warning("This is a warning")
        assert result.is_valid
        assert result.warning_count == 1
        assert result.error_count == 0

    def test_errors_invalidate(self):
        result = ValidationResult()
        result.add_error("This is an error")
        assert not result.is_valid
        assert result.error_count == 1
        assert result.warning_count == 0

    def test_merge(self):
        result1 = ValidationResult()
        result1.add_error("Error 1")
        result2 = ValidationResult()
        result2.add_warning("Warning 1")
        result1.merge(result2)
        assert result1.error_count == 1
        assert result1.warning_count == 1


class TestAlphaFold3Validation:
    def test_valid_alphafold3(self, tmp_path):
        payload = {
            "name": "test",
            "modelSeeds": [1],
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": "MLKK",
                    }
                }
            ],
            "dialect": "alphafold3",
            "version": 4,
        }
        path = tmp_path / "valid.json"
        path.write_text(json.dumps(payload))

        result = validate(str(path), "alphafold3")
        assert result.is_valid

    def test_missing_model_seeds(self, tmp_path):
        payload = {
            "name": "test",
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": "MLKK",
                    }
                }
            ],
            "dialect": "alphafold3",
            "version": 4,
        }
        path = tmp_path / "missing_seeds.json"
        path.write_text(json.dumps(payload))

        result = validate(str(path), "alphafold3")
        assert not result.is_valid
        assert any("modelSeeds" in i.message for i in result.issues)

    def test_empty_model_seeds(self, tmp_path):
        payload = {
            "name": "test",
            "modelSeeds": [],
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": "MLKK",
                    }
                }
            ],
            "dialect": "alphafold3",
            "version": 4,
        }
        path = tmp_path / "empty_seeds.json"
        path.write_text(json.dumps(payload))

        result = validate(str(path), "alphafold3")
        assert not result.is_valid
        assert any("at least one seed" in i.message for i in result.issues)

    def test_missing_sequence(self, tmp_path):
        payload = {
            "name": "test",
            "modelSeeds": [1],
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                    }
                }
            ],
            "dialect": "alphafold3",
            "version": 4,
        }
        path = tmp_path / "missing_seq.json"
        path.write_text(json.dumps(payload))

        result = validate(str(path), "alphafold3")
        assert not result.is_valid
        assert any("sequence" in i.message.lower() for i in result.issues)

    def test_duplicate_entity_id(self, tmp_path):
        payload = {
            "name": "test",
            "modelSeeds": [1],
            "sequences": [
                {"protein": {"id": "A", "sequence": "MLKK"}},
                {"protein": {"id": "A", "sequence": "GGGG"}},
            ],
            "dialect": "alphafold3",
            "version": 4,
        }
        path = tmp_path / "dup_id.json"
        path.write_text(json.dumps(payload))

        result = validate(str(path), "alphafold3")
        assert not result.is_valid
        assert any("Duplicate" in i.message for i in result.issues)

    def test_invalid_bond_entity_ref(self, tmp_path):
        payload = {
            "name": "test",
            "modelSeeds": [1],
            "sequences": [{"protein": {"id": "A", "sequence": "MLKK"}}],
            "bondedAtomPairs": [[["A", 1, "CA"], ["X", 1, "CA"]]],
            "dialect": "alphafold3",
            "version": 4,
        }
        path = tmp_path / "bad_bond.json"
        path.write_text(json.dumps(payload))

        result = validate(str(path), "alphafold3")
        assert not result.is_valid
        assert any("Entity 'X' not found" in i.message for i in result.issues)

    def test_wrong_dialect_warning(self, tmp_path):
        payload = {
            "name": "test",
            "modelSeeds": [1],
            "sequences": [{"protein": {"id": "A", "sequence": "MLKK"}}],
            "dialect": "alphafold2",
            "version": 4,
        }
        path = tmp_path / "wrong_dialect.json"
        path.write_text(json.dumps(payload))

        result = validate(str(path), "alphafold3")
        # Should pass but with warnings
        assert result.is_valid
        assert result.warning_count > 0
        assert any("dialect" in i.message.lower() for i in result.issues)

    def test_invalid_json(self, tmp_path):
        path = tmp_path / "invalid.json"
        path.write_text("{ not valid json }")

        result = validate(str(path), "alphafold3")
        assert not result.is_valid
        assert any("Failed to parse JSON" in i.message for i in result.issues)


class TestBoltz2Validation:
    def test_valid_boltz2(self, tmp_path):
        payload = """
version: 1
sequences:
  - protein:
      id: A
      sequence: MLKK
"""
        path = tmp_path / "valid.yaml"
        path.write_text(payload)

        result = validate(str(path), "boltz2")
        assert result.is_valid

    def test_missing_sequences(self, tmp_path):
        payload = """
version: 1
"""
        path = tmp_path / "no_seq.yaml"
        path.write_text(payload)

        result = validate(str(path), "boltz2")
        assert not result.is_valid
        assert any("sequences" in i.message for i in result.issues)

    def test_invalid_bond_entity(self, tmp_path):
        payload = """
version: 1
sequences:
  - protein:
      id: A
      sequence: MLKK
constraints:
  - bond:
      atom1: [A, 1, CA]
      atom2: [B, 1, CA]
"""
        path = tmp_path / "bad_bond.yaml"
        path.write_text(payload)

        result = validate(str(path), "boltz2")
        assert not result.is_valid
        assert any("Entity 'B' not found" in i.message for i in result.issues)


class TestAlphaFoldServerValidation:
    def test_valid_server_format(self, tmp_path):
        payload = [
            {
                "name": "test",
                "sequences": [{"proteinChain": {"sequence": "MLKK", "count": 1}}],
            }
        ]
        path = tmp_path / "valid.json"
        path.write_text(json.dumps(payload))

        result = validate(str(path), "alphafoldserver")
        assert result.is_valid

    def test_not_a_list(self, tmp_path):
        payload = {
            "name": "test",
            "sequences": [{"proteinChain": {"sequence": "MLKK", "count": 1}}],
        }
        path = tmp_path / "not_list.json"
        path.write_text(json.dumps(payload))

        result = validate(str(path), "alphafoldserver")
        assert not result.is_valid
        assert any("list of jobs" in i.message for i in result.issues)

    def test_empty_sequence(self, tmp_path):
        payload = [
            {
                "name": "test",
                "sequences": [{"proteinChain": {"sequence": "", "count": 1}}],
            }
        ]
        path = tmp_path / "empty_seq.json"
        path.write_text(json.dumps(payload))

        result = validate(str(path), "alphafoldserver")
        assert not result.is_valid
        assert any("empty" in i.message.lower() for i in result.issues)


class TestProtenixValidation:
    def test_valid_protenix(self, tmp_path):
        payload = [
            {
                "name": "test",
                "sequences": [{"proteinChain": {"sequence": "MLKK", "count": 1}}],
            }
        ]
        path = tmp_path / "valid.json"
        path.write_text(json.dumps(payload))

        result = validate(str(path), "protenix")
        assert result.is_valid

    def test_invalid_bond_entity_ref(self, tmp_path):
        payload = [
            {
                "name": "test",
                "sequences": [{"proteinChain": {"sequence": "MLKK", "count": 1}}],
                "covalent_bonds": [
                    {"entity1": 1, "position1": 1, "entity2": 5, "position2": 1}
                ],
            }
        ]
        path = tmp_path / "bad_bond.json"
        path.write_text(json.dumps(payload))

        result = validate(str(path), "protenix")
        assert not result.is_valid
        assert any("out of range" in i.message for i in result.issues)


class TestChai1Validation:
    def test_valid_chai1(self, tmp_path):
        fasta = """>protein|chain_A
MLKK
>ligand|ATP
ATP
"""
        path = tmp_path / "valid.fasta"
        path.write_text(fasta)

        result = validate(str(path), "chai1")
        assert result.is_valid

    def test_invalid_entity_type(self, tmp_path):
        fasta = """>unknown|chain_A
MLKK
"""
        path = tmp_path / "bad_type.fasta"
        path.write_text(fasta)

        result = validate(str(path), "chai1")
        assert not result.is_valid
        assert any("Invalid entity type" in i.message for i in result.issues)

    def test_empty_sequence(self, tmp_path):
        fasta = """>protein|chain_A

"""
        path = tmp_path / "empty.fasta"
        path.write_text(fasta)

        result = validate(str(path), "chai1")
        assert not result.is_valid
        assert any("empty" in i.message.lower() for i in result.issues)

    def test_valid_with_restraints(self, tmp_path):
        fasta = """>protein|chain_A
MLKK
>protein|chain_B
GGGG
"""
        restraints = """restraint_id,chainA,res_idxA,chainB,res_idxB,connection_type,confidence,min_distance_angstrom,max_distance_angstrom,comment
bond1,A,M1@CA,B,G1@CA,covalent,1.0,0.0,0.0,test
"""
        fasta_path = tmp_path / "valid.fasta"
        fasta_path.write_text(fasta)
        restraints_path = tmp_path / "restraints.csv"
        restraints_path.write_text(restraints)

        result = validate(str(fasta_path), "chai1", str(restraints_path))
        assert result.is_valid

    def test_invalid_connection_type(self, tmp_path):
        fasta = """>protein|chain_A
MLKK
"""
        restraints = """restraint_id,chainA,res_idxA,chainB,res_idxB,connection_type,confidence,min_distance_angstrom,max_distance_angstrom,comment
bond1,A,M1@CA,B,G1@CA,invalid_type,1.0,0.0,0.0,test
"""
        fasta_path = tmp_path / "test.fasta"
        fasta_path.write_text(fasta)
        restraints_path = tmp_path / "restraints.csv"
        restraints_path.write_text(restraints)

        result = validate(str(fasta_path), "chai1", str(restraints_path))
        assert not result.is_valid
        assert any("Invalid connection_type" in i.message for i in result.issues)


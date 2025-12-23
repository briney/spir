"""Tests for MSA path support across dialects."""

import json

import yaml

from spir.convert import ConvertOptions, convert


def test_alphafold3_msa_path_roundtrip(tmp_path):
    """Test AF3 -> Boltz -> AF3 preserves unpairedMsaPath."""
    payload = {
        "name": "msa_test",
        "modelSeeds": [42],
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": "MVLSPADKTN",
                    "unpairedMsaPath": "./msa/protein_a.a3m",
                }
            },
            {
                "rna": {
                    "id": "B",
                    "sequence": "ACGU",
                    "unpairedMsaPath": "./msa/rna_b.a3m",
                }
            },
        ],
        "dialect": "alphafold3",
        "version": 4,
    }
    in_path = tmp_path / "input.json"
    in_path.write_text(json.dumps(payload), encoding="utf-8")

    # AF3 -> Boltz
    boltz_prefix = tmp_path / "boltz"
    convert(str(in_path), "alphafold3", str(boltz_prefix), "boltz2", ConvertOptions())
    boltz_path = tmp_path / "boltz.yaml"
    boltz_data = yaml.safe_load(boltz_path.read_text(encoding="utf-8"))

    # Boltz should have msa for protein only
    protein_seq = boltz_data["sequences"][0]["protein"]
    assert protein_seq["msa"] == "./msa/protein_a.a3m"
    # RNA should not have msa in Boltz (Boltz only supports MSA for proteins)
    rna_seq = boltz_data["sequences"][1]["rna"]
    assert "msa" not in rna_seq

    # Boltz -> AF3
    out_prefix = tmp_path / "output"
    convert(str(boltz_path), "boltz2", str(out_prefix), "alphafold3", ConvertOptions())
    out_path = tmp_path / "output.json"
    out_data = json.loads(out_path.read_text(encoding="utf-8"))

    # AF3 output should have unpairedMsaPath for protein
    protein_out = out_data["sequences"][0]["protein"]
    assert protein_out["unpairedMsaPath"] == "./msa/protein_a.a3m"
    assert protein_out["pairedMsa"] == ""


def test_boltz_msa_path_roundtrip(tmp_path):
    """Test Boltz -> AF3 -> Boltz preserves msa path."""
    payload = {
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": "MVLSPADKTN",
                    "msa": "./examples/msa/seq1.a3m",
                }
            },
            {
                "ligand": {
                    "id": "L",
                    "ccd": "ATP",
                }
            },
        ],
    }
    in_path = tmp_path / "input.yaml"
    in_path.write_text(yaml.dump(payload), encoding="utf-8")

    # Boltz -> AF3
    af3_prefix = tmp_path / "af3"
    convert(str(in_path), "boltz2", str(af3_prefix), "alphafold3", ConvertOptions())
    af3_path = tmp_path / "af3.json"
    af3_data = json.loads(af3_path.read_text(encoding="utf-8"))

    protein_seq = af3_data["sequences"][0]["protein"]
    assert protein_seq["unpairedMsaPath"] == "./examples/msa/seq1.a3m"
    assert protein_seq["pairedMsa"] == ""

    # AF3 -> Boltz
    out_prefix = tmp_path / "output"
    convert(str(af3_path), "alphafold3", str(out_prefix), "boltz2", ConvertOptions())
    out_path = tmp_path / "output.yaml"
    out_data = yaml.safe_load(out_path.read_text(encoding="utf-8"))

    protein_out = out_data["sequences"][0]["protein"]
    assert protein_out["msa"] == "./examples/msa/seq1.a3m"


def test_af3_server_msa_path_to_boltz(tmp_path):
    """Test AF3 Server with non-standard msa_path converts to Boltz."""
    payload = [
        {
            "name": "msa_test",
            "modelSeeds": [42],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": "MVLSPADKTN",
                        "count": 1,
                        "msa_path": "./msa/protein.a3m",
                    }
                },
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        }
    ]
    in_path = tmp_path / "input.json"
    in_path.write_text(json.dumps(payload), encoding="utf-8")

    # AF3 Server -> Boltz
    boltz_prefix = tmp_path / "boltz"
    convert(str(in_path), "alphafoldserver", str(boltz_prefix), "boltz2", ConvertOptions())
    boltz_path = tmp_path / "boltz.yaml"
    boltz_data = yaml.safe_load(boltz_path.read_text(encoding="utf-8"))

    protein_seq = boltz_data["sequences"][0]["protein"]
    assert protein_seq["msa"] == "./msa/protein.a3m"


def test_af3_server_msa_path_to_alphafold3(tmp_path):
    """Test AF3 Server with non-standard msa_path converts to AF3."""
    payload = [
        {
            "name": "msa_test",
            "modelSeeds": [42],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": "MVLSPADKTN",
                        "count": 1,
                        "msa_path": "./msa/protein.a3m",
                    }
                },
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        }
    ]
    in_path = tmp_path / "input.json"
    in_path.write_text(json.dumps(payload), encoding="utf-8")

    # AF3 Server -> AF3
    af3_prefix = tmp_path / "af3"
    convert(str(in_path), "alphafoldserver", str(af3_prefix), "alphafold3", ConvertOptions())
    af3_path = tmp_path / "af3.json"
    af3_data = json.loads(af3_path.read_text(encoding="utf-8"))

    protein_seq = af3_data["sequences"][0]["protein"]
    assert protein_seq["unpairedMsaPath"] == "./msa/protein.a3m"
    assert protein_seq["pairedMsa"] == ""


def test_af3_server_output_never_includes_msa(tmp_path):
    """Test that AF3 Server output strictly adheres to spec (no MSA output)."""
    # Input with MSA paths
    payload = {
        "name": "msa_test",
        "modelSeeds": [42],
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": "MVLSPADKTN",
                    "unpairedMsaPath": "./msa/protein.a3m",
                }
            },
        ],
        "dialect": "alphafold3",
        "version": 4,
    }
    in_path = tmp_path / "input.json"
    in_path.write_text(json.dumps(payload), encoding="utf-8")

    # AF3 -> AF3 Server
    out_prefix = tmp_path / "output"
    convert(str(in_path), "alphafold3", str(out_prefix), "alphafoldserver", ConvertOptions())
    out_path = tmp_path / "output.json"
    out_data = json.loads(out_path.read_text(encoding="utf-8"))

    # AF3 Server output should NOT have any MSA fields
    protein_chain = out_data[0]["sequences"][0]["proteinChain"]
    assert "msa_path" not in protein_chain
    assert "unpairedMsaPath" not in protein_chain
    assert "unpairedMsa" not in protein_chain


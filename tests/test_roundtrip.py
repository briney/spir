import json
import os

from spir.convert import ConvertOptions, convert


def _extract_af3_server_glycans(payload):
    job = payload[0]
    glycans = []
    for entry in job.get("sequences", []):
        if "proteinChain" in entry:
            chain = entry["proteinChain"]
            for g in chain.get("glycans") or []:
                glycans.append((g["residues"], g["position"]))
    return glycans


def _af3_server_payload():
    return [
        {
            "name": "job1",
            "modelSeeds": [],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": "NVT",
                        "count": 1,
                        "glycans": [
                            {"residues": "NAG(NAG)(BMA)", "position": 1}
                        ],
                    }
                }
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        }
    ]


def test_roundtrip_af3_server_boltz_af3_server(tmp_path):
    in_path = tmp_path / "input.json"
    boltz_prefix = tmp_path / "out"
    out_prefix = tmp_path / "roundtrip"
    in_path.write_text(json.dumps(_af3_server_payload()), encoding="utf-8")

    convert(str(in_path), "alphafoldserver", str(boltz_prefix), "boltz2", ConvertOptions())
    boltz_path = tmp_path / "out.yaml"
    convert(str(boltz_path), "boltz2", str(out_prefix), "alphafoldserver", ConvertOptions())

    out_path = tmp_path / "roundtrip.json"
    out = json.loads(out_path.read_text(encoding="utf-8"))
    glycans = _extract_af3_server_glycans(out)
    assert glycans == [("NAG(NAG)(BMA)", 1)]


def test_roundtrip_af3_server_chai_af3_server(tmp_path):
    in_path = tmp_path / "input.json"
    chai_dir = tmp_path / "chai"
    chai_prefix = chai_dir / "case"
    out_prefix = tmp_path / "roundtrip"
    in_path.write_text(json.dumps(_af3_server_payload()), encoding="utf-8")

    convert(str(in_path), "alphafoldserver", str(chai_prefix), "chai1", ConvertOptions())
    assert os.path.exists(chai_dir / "case.fasta")
    assert os.path.exists(chai_dir / "case.constraints.csv")
    convert(str(chai_dir), "chai1", str(out_prefix), "alphafoldserver", ConvertOptions())

    out_path = tmp_path / "roundtrip.json"
    out = json.loads(out_path.read_text(encoding="utf-8"))
    glycans = _extract_af3_server_glycans(out)
    assert glycans == [("NAG(NAG)(BMA)", 1)]

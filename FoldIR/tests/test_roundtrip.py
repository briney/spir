import json
import os

from foldir.convert import ConvertOptions, convert


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
    boltz_path = tmp_path / "out.yaml"
    out_path = tmp_path / "roundtrip.json"
    in_path.write_text(json.dumps(_af3_server_payload()), encoding="utf-8")

    convert(str(in_path), "alphafoldserver", str(boltz_path), "boltz2", ConvertOptions())
    convert(str(boltz_path), "boltz2", str(out_path), "alphafoldserver", ConvertOptions())

    out = json.loads(out_path.read_text(encoding="utf-8"))
    glycans = _extract_af3_server_glycans(out)
    assert glycans == [("NAG(NAG)(BMA)", 1)]


def test_roundtrip_af3_server_chai_af3_server(tmp_path):
    in_path = tmp_path / "input.json"
    chai_dir = tmp_path / "chai"
    out_path = tmp_path / "roundtrip.json"
    in_path.write_text(json.dumps(_af3_server_payload()), encoding="utf-8")

    convert(str(in_path), "alphafoldserver", str(chai_dir), "chai1", ConvertOptions())
    assert os.path.exists(chai_dir / "chains.fasta")
    convert(str(chai_dir), "chai1", str(out_path), "alphafoldserver", ConvertOptions())

    out = json.loads(out_path.read_text(encoding="utf-8"))
    glycans = _extract_af3_server_glycans(out)
    assert glycans == [("NAG(NAG)(BMA)", 1)]

import json

from spir.convert import ConvertOptions, convert


def test_af3_server_to_af3(tmp_path):
    in_payload = [
        {
            "name": "job1",
            "modelSeeds": [],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": "N",
                        "count": 1,
                        "glycans": [
                            {"residues": "NAG(NAG)", "position": 1}
                        ],
                    }
                }
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        }
    ]
    in_path = tmp_path / "input.json"
    out_prefix = tmp_path / "output"
    in_path.write_text(json.dumps(in_payload), encoding="utf-8")

    convert(str(in_path), "alphafoldserver", str(out_prefix), "alphafold3", ConvertOptions())

    out_path = tmp_path / "output.json"
    out = json.loads(out_path.read_text(encoding="utf-8"))
    assert out["dialect"] == "alphafold3"
    assert out["modelSeeds"], "AF3 output should include modelSeeds"
    sequences = out["sequences"]
    ligands = [s["ligand"] for s in sequences if "ligand" in s]
    assert any(l.get("ccdCodes") == ["NAG", "NAG"] for l in ligands)
    bonded = out.get("bondedAtomPairs") or []
    assert len(bonded) == 2

from __future__ import annotations

import json

import yaml

from spir import convert_af3_server_glycan

from ..io import read
from ..ir.glycan import parse_af3_server_glycan, render_chai_glycan
from ..ir.ir import IR


def test_minimal_af3_snippet_contains_ccdcodes():
    s = convert_af3_server_glycan("NAG(NAG)(BMA)", to="af3", detail="minimal")
    obj = json.loads(s)
    assert obj["ligand"]["ccdCodes"] == ["NAG", "NAG", "BMA"]


def test_minimal_af3_server_is_echo():
    s = convert_af3_server_glycan(
        "  NAG ( NAG ) ( BMA )  ", to="af3-server", detail="minimal"
    )
    assert s == "NAG ( NAG ) ( BMA )".strip()


def test_minimal_protenix_has_joined_ccd():
    s = convert_af3_server_glycan("NAG(NAG)(BMA)", to="protenix", detail="minimal")
    obj = json.loads(s)
    assert obj["ligand"]["ligand"] == "CCD_NAG_NAG_BMA"


def test_minimal_boltz_uses_first_ccd():
    s = convert_af3_server_glycan("NAG(NAG)(BMA)", to="boltz", detail="minimal")
    assert "ccd: NAG" in s
    assert "components: [NAG, NAG, BMA]" in s


def test_minimal_chai_outputs_full_glycan():
    inp = "NAG(NAG)(BMA)"
    s = convert_af3_server_glycan(inp, to="chai", detail="minimal")
    expected = render_chai_glycan(parse_af3_server_glycan(inp))
    assert s.strip().splitlines() == [">glycan|G", expected]


def _split_chai_combo(text: str):
    parts = text.split("\n---- restraints.csv ----\n")
    if len(parts) == 1:
        return parts[0], ""
    return parts[0], parts[1]


def test_full_linear_conversion_all_formats():
    gly = "NAG(NAG(MAN(MAN(MAN))))"  # 5 components linear
    n = 5

    # AF3
    s = convert_af3_server_glycan(gly, to="af3", detail="full")
    obj = json.loads(s)
    assert obj["dialect"] == "alphafold3"
    # ligand with all CCDs
    ligs = [e["ligand"] for e in obj["sequences"] if "ligand" in e]
    assert any(
        l.get("id") == "G" and l.get("ccdCodes") and len(l["ccdCodes"]) == n
        for l in ligs
    )
    # bonds = n-1 intra + 1 anchor
    assert "bondedAtomPairs" in obj and len(obj["bondedAtomPairs"]) == (n - 1 + 1)

    # anchor must involve protein A@ND2 to glycan root @C1 (any glycan chain id)
    def _is_anchor(pair):
        (e1, r1, a1), (e2, r2, a2) = pair
        return (e1 == "A" and r1 == 1 and a1 == "ND2" and a2 == "C1" and r2 == 1) or (
            e2 == "A" and r2 == 1 and a2 == "ND2" and a1 == "C1" and r1 == 1
        )

    assert any(_is_anchor(pair) for pair in obj["bondedAtomPairs"])

    # Boltz
    y = convert_af3_server_glycan(gly, to="boltz", detail="full")
    yobj = yaml.safe_load(y)
    # n ligands + 1 protein
    lig_seq = [e for e in yobj.get("sequences", []) if "ligand" in e]
    assert len(lig_seq) == n
    bonds = [c for c in yobj.get("constraints", []) if "bond" in c]
    assert len(bonds) == (n - 1 + 1)
    assert any(
        b["bond"]["atom1"][0] == "A" and b["bond"]["atom1"][2] == "ND2" for b in bonds
    )

    # Chai
    c = convert_af3_server_glycan(gly, to="chai", detail="full")
    fasta, restraints = _split_chai_combo(c)
    # FASTA glycan must equal rendered spec
    expected = render_chai_glycan(parse_af3_server_glycan(gly))
    flines = [ln for ln in fasta.strip().splitlines() if not ln.startswith(">protein")]
    assert any(ln == expected for ln in flines)
    # Restraints: only one data row and it is the anchor to @C1
    rdata = [
        ln
        for ln in restraints.strip().splitlines()
        if ln and not ln.startswith("chainA,")
    ]
    assert (
        len(rdata) == 1 and ",".join(["A", "N1@N", "G", "@C1", "covalent"]) in rdata[0]
    )

    # Protenix
    p = convert_af3_server_glycan(gly, to="protenix", detail="full")
    parr = json.loads(p)
    pobj = parr[0]
    cb = pobj.get("covalent_bonds", [])
    assert len(cb) == (n - 1 + 1)
    assert any(
        b.get("atom1") == "ND2" and b.get("atom2") == "C1" and b.get("position2") == "1"
        for b in cb
    )


def test_full_branched_conversion_all_formats():
    gly = "NAG(NAG(MAN(MAN(MAN)(MAN(MAN)(MAN)))))"
    # Count CCD tokens:
    n = sum(
        1 for ch in gly if ch == "A" or ch == "M"
    )  # crude but stable for NAG/MAN examples
    # AF3
    s = convert_af3_server_glycan(gly, to="af3", detail="full")
    obj = json.loads(s)
    assert obj["dialect"] == "alphafold3"
    ligs = [e["ligand"] for e in obj["sequences"] if "ligand" in e]
    assert any(
        l.get("id") == "G" and l.get("ccdCodes") and len(l["ccdCodes"]) >= 5
        for l in ligs
    )
    # bonds at least n-1 intra + 1 anchor (parser yields a tree)
    assert "bondedAtomPairs" in obj and len(obj["bondedAtomPairs"]) >= 5

    def _is_anchor2(pair):
        (e1, r1, a1), (e2, r2, a2) = pair
        return (e1 == "A" and r1 == 1 and a1 == "ND2" and a2 == "C1" and r2 == 1) or (
            e2 == "A" and r2 == 1 and a2 == "ND2" and a1 == "C1" and r1 == 1
        )

    assert any(_is_anchor2(pair) for pair in obj["bondedAtomPairs"])
    # Boltz
    y = convert_af3_server_glycan(gly, to="boltz", detail="full")
    yobj = yaml.safe_load(y)
    bonds = [c for c in yobj.get("constraints", []) if "bond" in c]
    assert any(
        b["bond"]["atom1"][0] == "A" and b["bond"]["atom1"][2] == "ND2" for b in bonds
    )
    # Chai
    c = convert_af3_server_glycan(gly, to="chai", detail="full")
    fasta2, restraints2 = _split_chai_combo(c)
    expected2 = render_chai_glycan(parse_af3_server_glycan(gly))
    flines2 = [
        ln for ln in fasta2.strip().splitlines() if not ln.startswith(">protein")
    ]
    assert any(ln == expected2 for ln in flines2)
    rdata2 = [
        ln
        for ln in restraints2.strip().splitlines()
        if ln and not ln.startswith("chainA,")
    ]
    assert len(rdata2) == 1 and "A,N1@N,G,@C1,covalent" in rdata2[0]
    # Protenix
    p = convert_af3_server_glycan(gly, to="protenix", detail="full")
    parr = json.loads(p)
    pobj = parr[0]
    cb = pobj.get("covalent_bonds", [])
    assert any(b.get("atom1") == "ND2" and b.get("atom2") == "C1" for b in cb)


def test_af3_server_to_all_formats_preserves_bonds_linear():
    # Build AF3-Server JSON with one protein and one glycan
    job = [
        {
            "name": "bonded",
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": "ACDENST",
                        "count": 1,
                        "glycans": [
                            {"residues": "NAG(NAG(MAN(MAN(MAN))))", "position": 5}
                        ],
                    }
                }
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        }
    ]
    text = json.dumps(job)
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "in.af3server.json"
        p.write_text(text)
        ir = read(str(p))

        # AF3
        af3 = ir.write_alphafold3(td)
        af3_obj = json.loads(Path(af3).read_text())
        assert af3_obj.get("bondedAtomPairs")

        def _is_anchor_lin(pair):
            (e1, r1, a1), (e2, r2, a2) = pair
            return (
                e1 == "A" and r1 == 5 and a1 == "ND2" and a2 == "C1" and r2 == 1
            ) or (e2 == "A" and r2 == 5 and a2 == "ND2" and a1 == "C1" and r1 == 1)

        assert any(_is_anchor_lin(pair) for pair in af3_obj["bondedAtomPairs"])

        # Protenix
        px = ir.write_protenix(td)
        px_arr = json.loads(Path(px).read_text())
        assert px_arr[0].get("covalent_bonds")

        # Boltz
        bol = ir.write_boltz(td)
        bol_yaml = yaml.safe_load(Path(bol).read_text())
        bonds = [c for c in bol_yaml.get("constraints", []) if "bond" in c]
        assert bonds and any(b["bond"]["atom1"][2] == "ND2" for b in bonds)

        # Chai
        fasta_path, restraints_path = ir.write_chai(td)
        r = Path(restraints_path).read_text()
        assert "covalent" in r and "A,N5@N" in r

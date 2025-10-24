from __future__ import annotations

from ..io import read


def test_smoke_af3_roundtrip():
    inp = {
        "name": "hello",
        "modelSeeds": [1],
        "sequences": [
            {"protein": {"id": "A", "sequence": "ACDE"}},
            {"ligand": {"id": "L", "ccdCodes": ["ATP"]}},
        ],
        "dialect": "alphafold3",
        "version": 4,
    }
    text = __import__("json").dumps(inp)
    # Write temp file-less path: load via adapter directly, then use IR writer for smoke
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "in.af3.json"
        p.write_text(text)
        ir = read(str(p))
        out_path = Path(td) / "out"
        out_path.mkdir(parents=True, exist_ok=True)
        path = ir.write_alphafold3(out_path)
        out = path.read_text()
    assert "alphafold3" in out


def test_smoke_af3_server_write_with_glycan():
    # Build a minimal IR via AF3 reader and then attach a glycan via bonds
    import tempfile
    from pathlib import Path

    from ..adapters.af3 import AF3ReaderWriter
    from ..ir.model import AtomRef, Bond, ComplexInput, Ligand, PolymerChain

    ci = ComplexInput(
        name="hello-server",
        seeds=[123],
        proteins=[PolymerChain(type="protein", ids=["A"], sequence="ACDEFGHIK")],
        ligands=[
            Ligand(
                ids=["G"],
                ccd_codes=["NAG", "NAG", "MAN", "MAN", "MAN"],
                af3_residues="NAG(NAG(MAN(MAN(MAN))))",
            )
        ],
        bonds=[
            Bond(
                atom1=AtomRef(chain_id="A", residue_index=8),
                atom2=AtomRef(chain_id="G", component_index=1),
            )
        ],
    )

    from ..ir.ir import IR

    ir = IR(ci=ci, source_paths=(), source_format="af3")
    with tempfile.TemporaryDirectory() as td:
        out_dir = Path(td)
        path = ir.write_alphafoldserver(out_dir)
        text = path.read_text()
        assert "alphafoldserver" in text
        assert '"residues": "NAG(NAG(MAN(MAN(MAN))))"' in text

        # Also write Chai and verify FASTA glycan and single anchor row exist
        fasta_path, restraints_path = ir.write_chai(out_dir)
        fasta_text = fasta_path.read_text()
        assert ">glycan|G" in fasta_text and "NAG(NAG(MAN(MAN(MAN))))" in fasta_text
        if restraints_path:
            r_text = restraints_path.read_text()
            data_rows = [
                ln
                for ln in r_text.strip().splitlines()
                if ln and not ln.startswith("chainA,")
            ]
            assert len(data_rows) == 1
            # AF3-Server write will choose chain id 'G' for glycan in this smoke test
            row = data_rows[0]
            # Expect ResidueLetter+Position (I8 here); atom letter may be absent for non-glyco residues
            assert row.startswith("A,I8") and (
                ",G,@C1,covalent" in row or ",B,@C1,covalent" in row
            )

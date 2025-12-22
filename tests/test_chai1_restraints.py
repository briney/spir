from foldir.convert import ConvertOptions, convert


def test_chai1_restraints_explicit_path(tmp_path):
    fasta_path = tmp_path / "input.fasta"
    restraints_path = tmp_path / "restraints.csv"
    out_prefix = tmp_path / "out"

    fasta_path.write_text(
        ">protein|chain1\nACDE\n>protein|chain2\nFGHI\n",
        encoding="utf-8",
    )
    restraints_path.write_text(
        "restraint_id,chainA,res_idxA,chainB,res_idxB,connection_type,confidence,"
        "min_distance_angstrom,max_distance_angstrom,comment\n"
        "restraint1,A,A1@CA,B,B1@CA,contact,1.0,0.0,8.0,test\n",
        encoding="utf-8",
    )

    convert(
        str(fasta_path),
        "chai1",
        str(out_prefix),
        "chai1",
        ConvertOptions(),
        restraints_path=str(restraints_path),
    )

    out_constraints = tmp_path / "out.constraints.csv"
    assert out_constraints.exists()
    assert "contact" in out_constraints.read_text(encoding="utf-8")

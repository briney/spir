from foldir.ir.glycans.parse_af3_server import parse_af3_server_glycan_string
from foldir.ir.glycans.parse_chai import parse_chai_glycan_string
from foldir.ir.glycans.render_af3_server import render_af3_server_glycan_string
from foldir.ir.glycans.render_chai import render_chai_glycan_string


def test_af3_server_roundtrip():
    s = "NAG(NAG)(BMA)"
    g = parse_af3_server_glycan_string("g1", s)
    out = render_af3_server_glycan_string(g, "g1.n0")
    assert out == s


def test_chai_roundtrip():
    s = "NAG(4-1 NAG(4-1 BMA(3-1 MAN)(6-1 MAN)))"
    g = parse_chai_glycan_string("g2", s)
    out = render_chai_glycan_string(g, "g2.n0")
    assert out == s

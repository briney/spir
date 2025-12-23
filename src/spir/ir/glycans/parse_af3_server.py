from __future__ import annotations

import re
from typing import List

from spir.ir.models import Glycan, GlycanEdge, GlycanNode

CCD_RE = re.compile(r"[A-Za-z0-9]{3}")


class ParseError(ValueError):
    pass


def parse_af3_server_glycan_string(glycan_id: str, s: str) -> Glycan:
    """
    Parse AF3 Server compact glycan tree string into a Glycan graph.
    Linkage atoms are unknown in this format, so parent_atom is None.
    """

    i = 0
    nodes: List[GlycanNode] = []
    edges: List[GlycanEdge] = []
    node_counter = 0

    def parse_node() -> str:
        nonlocal i, node_counter
        m = CCD_RE.match(s, i)
        if not m:
            raise ParseError(f"Expected CCD code at offset {i}: ...{s[i:i+10]!r}")
        ccd = m.group(0)
        i = m.end()

        node_id = f"{glycan_id}.n{node_counter}"
        node_counter += 1
        nodes.append(GlycanNode(node_id=node_id, ccd=ccd))

        while i < len(s) and s[i] == "(":
            i += 1
            child_id = parse_node()
            if i >= len(s) or s[i] != ")":
                raise ParseError(f"Missing ')' at offset {i}")
            i += 1
            edges.append(GlycanEdge(parent=node_id, child=child_id, parent_atom=None, child_atom="C1"))

        return node_id

    parse_node()
    if i != len(s):
        raise ParseError(f"Trailing junk at offset {i}: {s[i:]!r}")

    return Glycan(glycan_id=glycan_id, nodes=nodes, edges=edges, attachments=[])

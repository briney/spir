from __future__ import annotations

import re
from typing import List

from foldir.ir.models import Glycan, GlycanEdge, GlycanNode

CCD_RE = re.compile(r"[A-Za-z0-9]{3}")
INT_RE = re.compile(r"\d+")


class ParseError(ValueError):
    pass


def parse_chai_glycan_string(glycan_id: str, s: str) -> Glycan:
    i = 0
    nodes: List[GlycanNode] = []
    edges: List[GlycanEdge] = []
    node_counter = 0

    def read_ccd() -> str:
        nonlocal i
        m = CCD_RE.match(s, i)
        if not m:
            raise ParseError(f"Expected CCD at {i}")
        i = m.end()
        return m.group(0)

    def read_int() -> int:
        nonlocal i
        m = INT_RE.match(s, i)
        if not m:
            raise ParseError(f"Expected int at {i}")
        i = m.end()
        return int(m.group(0))

    def skip_ws() -> None:
        nonlocal i
        while i < len(s) and s[i].isspace():
            i += 1

    def parse_node() -> str:
        nonlocal i, node_counter
        ccd = read_ccd()
        node_id = f"{glycan_id}.n{node_counter}"
        node_counter += 1
        nodes.append(GlycanNode(node_id=node_id, ccd=ccd))

        while i < len(s) and s[i] == "(":
            i += 1
            parent_pos = read_int()
            if i >= len(s) or s[i] != "-":
                raise ParseError(f"Expected '-' after parent_pos at {i}")
            i += 1
            child_pos = read_int()
            skip_ws()
            child_id = parse_node()
            if i >= len(s) or s[i] != ")":
                raise ParseError(f"Missing ')' at {i}")
            i += 1

            edges.append(
                GlycanEdge(
                    parent=node_id,
                    child=child_id,
                    parent_atom=f"O{parent_pos}",
                    child_atom=f"C{child_pos}",
                )
            )

        return node_id

    parse_node()
    if i != len(s):
        raise ParseError(f"Trailing junk: {s[i:]!r}")

    return Glycan(glycan_id=glycan_id, nodes=nodes, edges=edges, attachments=[])

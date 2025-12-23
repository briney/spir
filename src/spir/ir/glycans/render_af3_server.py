from __future__ import annotations

from collections import defaultdict

from spir.ir.models import Glycan


class RenderError(ValueError):
    pass


def render_af3_server_glycan_string(g: Glycan, root_node_id: str) -> str:
    nodes = {n.node_id: n for n in g.nodes}
    children = defaultdict(list)
    parents = {}
    for e in g.edges:
        children[e.parent].append(e.child)
        if e.child in parents:
            raise RenderError("Not a tree (child has multiple parents).")
        parents[e.child] = e.parent

    if len(g.nodes) > 8:
        raise RenderError("AF3 Server supports up to 8 glycan residues.")
    for _, kids in children.items():
        if len(kids) > 2:
            raise RenderError("AF3 Server glycan nodes may have at most 2 children.")

    def rec(nid: str) -> str:
        s = nodes[nid].ccd
        for kid in children.get(nid, []):
            s += f"({rec(kid)})"
        return s

    return rec(root_node_id)

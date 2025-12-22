from __future__ import annotations

from collections import defaultdict

from foldir.ir.models import Glycan


class RenderError(ValueError):
    pass


def _atom_to_pos(atom: str) -> int:
    if len(atom) < 2 or not atom[1:].isdigit():
        raise RenderError(f"Cannot convert atom to pos: {atom}")
    return int(atom[1:])


def render_chai_glycan_string(g: Glycan, root_node_id: str) -> str:
    nodes = {n.node_id: n for n in g.nodes}
    children = defaultdict(list)
    parents = {}
    for e in g.edges:
        if e.parent_atom is None or e.child_atom is None:
            raise RenderError("Chai requires explicit linkage positions/atoms.")
        children[e.parent].append((e, e.child))
        if e.child in parents:
            raise RenderError("Not a tree (child has multiple parents).")
        parents[e.child] = e.parent

    def rec(nid: str) -> str:
        s = nodes[nid].ccd
        for edge, kid in children.get(nid, []):
            ppos = _atom_to_pos(edge.parent_atom)
            cpos = _atom_to_pos(edge.child_atom)
            s += f"({ppos}-{cpos} {rec(kid)})"
        return s

    return rec(root_node_id)

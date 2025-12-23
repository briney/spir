from __future__ import annotations

from typing import Protocol

from spir.ir.models import Glycan


class LinkageResolver(Protocol):
    def resolve(self, parent_ccd: str, child_ccd: str) -> tuple[str, str]:
        """
        Return (parent_atom, child_atom), e.g. ("O4", "C1").
        """


class DefaultSugarLinkageResolver:
    """
    Conservative defaults:
    - use child_atom = C1 (common)
    - use parent_atom = O4 (common)
    """

    def __init__(self, default_parent_atom: str = "O4", default_child_atom: str = "C1"):
        self.default_parent_atom = default_parent_atom
        self.default_child_atom = default_child_atom

    def resolve(self, parent_ccd: str, child_ccd: str) -> tuple[str, str]:
        return (self.default_parent_atom, self.default_child_atom)


def fill_missing_linkages(g: Glycan, resolver: LinkageResolver) -> Glycan:
    nodes = {n.node_id: n.ccd for n in g.nodes}
    new_edges = []
    for e in g.edges:
        if e.parent_atom is None or e.child_atom is None:
            p_atom, c_atom = resolver.resolve(nodes[e.parent], nodes[e.child])
            new_edges.append(e.model_copy(update={"parent_atom": p_atom, "child_atom": c_atom}))
        else:
            new_edges.append(e)
    return g.model_copy(update={"edges": new_edges})

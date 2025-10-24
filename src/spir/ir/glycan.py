from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple


@dataclass
class GlycanComponent:
    ccd: str


@dataclass
class GlycanBond:
    # indices are 1-based component indices in the glycan component list
    parent_index: int
    parent_atom: str  # e.g., "O4"
    child_index: int
    child_atom: str  # e.g., "C1"


@dataclass
class GlycanGraph:
    components: List[GlycanComponent]
    bonds: List[GlycanBond]

    def components_as_ccd_list(self) -> List[str]:
        return [c.ccd for c in self.components]


class GlycanParseError(ValueError):
    pass


def parse_chai_glycan(spec: str) -> GlycanGraph:
    """Parse a Chai glycan specification like
    NAG(4-1 NAG(4-1 BMA(3-1 MAN)(6-1 MAN)))

    Returns a GlycanGraph with components in preorder and bonds as annotated edges
    from parent to child using the link numbers (O<left>-C<right> convention).
    """
    tokens = _tokenize(spec)
    index = 0
    components: List[GlycanComponent] = []
    bonds: List[GlycanBond] = []

    def parse_node(parent_idx: Optional[int] = None) -> int:
        nonlocal index
        # Expect CCD code
        if index >= len(tokens) or tokens[index][0] != "CCD":
            raise GlycanParseError("Expected CCD code at position %d" % index)
        ccd = tokens[index][1]
        index += 1
        this_idx = len(components) + 1
        components.append(GlycanComponent(ccd=ccd))
        # Parse zero or more branches in parentheses
        while index < len(tokens) and tokens[index][0] == "LPAREN":
            # '(', link, child, ')'
            index += 1
            # Link token like '4-1'
            if index >= len(tokens) or tokens[index][0] != "LINK":
                raise GlycanParseError(
                    "Expected link like '4-1' after '(' at %d" % index
                )
            left, right = tokens[index][1]
            index += 1
            # Parse child subtree
            child_idx = parse_node(parent_idx=this_idx)
            # Record bond: O<left> (parent) to C<right> (child)
            bonds.append(
                GlycanBond(
                    parent_index=this_idx,
                    parent_atom=f"O{left}",
                    child_index=child_idx,
                    child_atom=f"C{right}",
                )
            )
            # Expect ')'
            if index >= len(tokens) or tokens[index][0] != "RPAREN":
                raise GlycanParseError("Expected ')' after child at %d" % index)
            index += 1
        return this_idx

    root_idx = parse_node(None)
    if index != len(tokens):
        raise GlycanParseError("Unexpected trailing tokens starting at %d" % index)
    return GlycanGraph(components=components, bonds=bonds)


def render_chai_glycan(graph: GlycanGraph) -> str:
    """Render a GlycanGraph back to Chai glycan syntax.
    We assume the bonds form a tree rooted at component 1 and that bonds are oriented parent->child.
    """
    children: List[List[Tuple[int, str, str]]] = [
        list() for _ in range(len(graph.components) + 1)
    ]
    for b in graph.bonds:
        children[b.parent_index].append((b.child_index, b.parent_atom, b.child_atom))

    def render(idx: int) -> str:
        ccd = graph.components[idx - 1].ccd
        branches = []
        for child_idx, parent_atom, child_atom in children[idx]:
            # parent_atom is like O4; child_atom is like C1
            l = parent_atom[1:] if parent_atom.startswith("O") else parent_atom
            r = child_atom[1:] if child_atom.startswith("C") else child_atom
            branches.append(f"({l}-{r} {render(child_idx)})")
        return ccd + "".join(branches)

    return render(1) if graph.components else ""


def _tokenize(spec: str):
    s = spec.strip()
    i = 0
    tokens = []
    while i < len(s):
        ch = s[i]
        if ch.isspace():
            i += 1
            continue
        if ch == "(":
            tokens.append(("LPAREN", ch))
            i += 1
            continue
        if ch == ")":
            tokens.append(("RPAREN", ch))
            i += 1
            continue
        # LINK like 4-1
        if ch.isdigit():
            j = i
            while j < len(s) and s[j].isdigit():
                j += 1
            if j < len(s) and s[j] == "-":
                k = j + 1
                while k < len(s) and s[k].isdigit():
                    k += 1
                if k > j + 1:
                    tokens.append(("LINK", (s[i:j], s[j + 1 : k])))
                    i = k
                    continue
            raise GlycanParseError(f"Invalid link starting at position {i}")
        # CCD code: uppercase letters/numbers/"-"
        if ch.isalpha():
            j = i
            while j < len(s) and (s[j].isalnum() or s[j] in {"-", "_"}):
                j += 1
            ccd = s[i:j]
            tokens.append(("CCD", ccd))
            i = j
            continue
        raise GlycanParseError(f"Unexpected character '{ch}' at position {i}")
    return tokens


# ---- AF3-Server glycan parsing (parentheses-only; no explicit link tokens) ----


def parse_af3_server_glycan(spec: str) -> GlycanGraph:
    """Parse an AF3-Server style glycan string like
    NAG(NAG(MAN(MAN(MAN)))) or NAG(NAG)(BMA)

    Returns a GlycanGraph where components are listed in preorder and bonds are
    inferred deterministically per linkage policy, oriented parent->child with
    parent O[n] ↔ child C1.
    """
    tokens = _tokenize_af3server(spec)
    index = 0
    components: List[GlycanComponent] = []
    bonds: List[GlycanBond] = []

    def parse_node() -> int:
        nonlocal index
        # Expect CCD token
        if index >= len(tokens) or tokens[index][0] != "CCD":
            raise GlycanParseError("Expected CCD code at position %d" % index)
        ccd = tokens[index][1]
        index += 1
        this_idx = len(components) + 1
        components.append(GlycanComponent(ccd=ccd))

        child_count = 0
        # Zero or more child subtrees in parentheses
        while index < len(tokens) and tokens[index][0] == "LPAREN":
            index += 1  # consume '('
            child_idx = parse_node()
            # Infer linkage once child is parsed
            parent_ccd = components[this_idx - 1].ccd
            child_ccd = components[child_idx - 1].ccd
            is_trunk = child_count == 0
            branch_idx = 0 if is_trunk else child_count  # 1-based for branches
            parent_atom = _infer_parent_atom_for_link(
                parent_ccd, child_ccd, is_trunk, branch_idx
            )
            bonds.append(
                GlycanBond(
                    parent_index=this_idx,
                    parent_atom=parent_atom,
                    child_index=child_idx,
                    child_atom="C1",
                )
            )
            # expect ')'
            if index >= len(tokens) or tokens[index][0] != "RPAREN":
                raise GlycanParseError(
                    "Expected ')' after child at token index %d" % index
                )
            index += 1
            child_count += 1
        return this_idx

    root_idx = parse_node()
    if index != len(tokens):
        raise GlycanParseError("Unexpected trailing tokens starting at %d" % index)
    # Sanity: require tree
    if root_idx != 1:
        raise GlycanParseError("Invalid glycan root index")
    return GlycanGraph(components=components, bonds=bonds)


def _tokenize_af3server(spec: str):
    s = spec.strip()
    i = 0
    tokens = []
    while i < len(s):
        ch = s[i]
        if ch.isspace():
            i += 1
            continue
        if ch == "(":
            tokens.append(("LPAREN", ch))
            i += 1
            continue
        if ch == ")":
            tokens.append(("RPAREN", ch))
            i += 1
            continue
        # CCD code: letters/numbers/"-"/"_"
        if ch.isalpha():
            j = i
            while j < len(s) and (s[j].isalnum() or s[j] in {"-", "_"}):
                j += 1
            ccd = s[i:j]
            tokens.append(("CCD", ccd))
            i = j
            continue
        raise GlycanParseError(f"Unexpected character '{ch}' at position {i}")
    return tokens


def _infer_parent_atom_for_link(
    parent_ccd: str, child_ccd: str, is_trunk: bool, branch_idx: int
) -> str:
    """Deterministically assign parent atom O[n] for the parent→child bond.

    - For trunk edges: use specific defaults; MAN→MAN = O2; NAG→NAG/NAG→MAN = O4; else O4.
    - For branch edges: if parent is mannose (MAN/BMA), use [3, 6, 4, 2, 4, 2, ...]; else [3, 6, 4, 2, 4, 2, ...].
    Returns the parent atom name like "O4".
    """
    p = parent_ccd.upper()
    c = child_ccd.upper()

    def is_mannose(code: str) -> bool:
        return code in {"MAN", "BMA"}

    if is_trunk:
        if is_mannose(p) and is_mannose(c):
            n = 2
        elif p == "NAG" and (c == "NAG" or is_mannose(c)):
            n = 4
        else:
            n = 4
        return f"O{n}"

    # branch_idx is 1-based for branches
    idx = max(1, branch_idx)
    if is_mannose(p):
        if idx == 1:
            n = 3
        elif idx == 2:
            n = 6
        else:
            n = 4 if (idx - 3) % 2 == 0 else 2
    else:
        if idx == 1:
            n = 3
        elif idx == 2:
            n = 6
        else:
            n = 4 if (idx - 3) % 2 == 0 else 2
    return f"O{n}"

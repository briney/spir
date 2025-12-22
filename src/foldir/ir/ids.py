from __future__ import annotations

from typing import Iterable, List, Set, TypeVar

T = TypeVar("T")


def _dedupe_ids(items: Iterable[T], seen: Set[str], prefix: str, attr: str = "id") -> List[T]:
    out: List[T] = []
    counter = 1
    for item in items:
        value = getattr(item, attr)
        new_value = value
        if not new_value or new_value in seen:
            base = new_value if new_value else prefix
            while True:
                candidate = f"{base}{counter}"
                counter += 1
                if candidate not in seen:
                    new_value = candidate
                    break
        seen.add(new_value)
        if new_value != value:
            out.append(item.model_copy(update={attr: new_value}))
        else:
            out.append(item)
    return out


def ensure_unique_entity_ids(
    polymers: Iterable[T],
    ligands: Iterable[T],
    ions: Iterable[T],
) -> tuple[List[T], List[T], List[T]]:
    seen: Set[str] = set()
    new_polymers = _dedupe_ids(polymers, seen, "P")
    new_ligands = _dedupe_ids(ligands, seen, "L")
    new_ions = _dedupe_ids(ions, seen, "I")
    return new_polymers, new_ligands, new_ions


def ensure_unique_glycan_ids(glycans: Iterable[T]) -> List[T]:
    seen: Set[str] = set()
    return _dedupe_ids(glycans, seen, "G", attr="glycan_id")

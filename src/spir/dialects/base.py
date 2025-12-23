from __future__ import annotations

from typing import Protocol

from spir.ir.models import DocumentIR


class Dialect(Protocol):
    name: str

    def parse(self, path: str) -> DocumentIR:
        ...

    def render(self, doc: DocumentIR, out_path: str) -> None:
        ...

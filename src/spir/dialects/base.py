from __future__ import annotations

from typing import TYPE_CHECKING, Protocol

from spir.ir.models import DocumentIR

if TYPE_CHECKING:
    from spir.validate import ValidationResult


class Dialect(Protocol):
    name: str

    def parse(self, path: str) -> DocumentIR:
        ...

    def render(self, doc: DocumentIR, out_path: str) -> None:
        ...

    def validate(self, path: str) -> ValidationResult:
        ...

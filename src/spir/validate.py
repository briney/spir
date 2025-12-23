from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING, List, Optional

import typer

if TYPE_CHECKING:
    pass


class Severity(str, Enum):
    ERROR = "error"
    WARNING = "warning"


@dataclass
class ValidationIssue:
    """A single validation issue found in an input file."""

    severity: Severity
    message: str
    location: Optional[str] = None

    def __str__(self) -> str:
        loc = f" at {self.location}" if self.location else ""
        return f"[{self.severity.value.upper()}]{loc} {self.message}"


@dataclass
class ValidationResult:
    """Result of validating an input file."""

    issues: List[ValidationIssue] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        """Returns True if there are no errors (warnings are acceptable)."""
        return not any(i.severity == Severity.ERROR for i in self.issues)

    @property
    def error_count(self) -> int:
        return sum(1 for i in self.issues if i.severity == Severity.ERROR)

    @property
    def warning_count(self) -> int:
        return sum(1 for i in self.issues if i.severity == Severity.WARNING)

    def add_error(self, message: str, location: Optional[str] = None) -> None:
        self.issues.append(ValidationIssue(Severity.ERROR, message, location))

    def add_warning(self, message: str, location: Optional[str] = None) -> None:
        self.issues.append(ValidationIssue(Severity.WARNING, message, location))

    def merge(self, other: ValidationResult) -> None:
        """Merge another ValidationResult into this one."""
        self.issues.extend(other.issues)


def validate(
    path: str, dialect_name: str, restraints_path: Optional[str] = None
) -> ValidationResult:
    """
    Validate an input file against a dialect's schema and rules.

    Args:
        path: Path to the input file
        dialect_name: Name of the dialect to validate against
        restraints_path: Optional path to restraints file (for Chai-1)

    Returns:
        ValidationResult containing any issues found
    """
    # Import here to avoid circular imports
    from spir.dialects import get_dialect

    dialect = get_dialect(dialect_name)
    if restraints_path:
        return dialect.validate(path, restraints_path)
    return dialect.validate(path)


def print_validation_result(result: ValidationResult, path: str, dialect: str) -> None:
    """Print validation results to the console."""
    typer.echo(f"Validating: {path} (dialect: {dialect})")
    typer.echo()

    if not result.issues:
        typer.secho("✓ Validation passed", fg=typer.colors.GREEN)
        return

    for issue in result.issues:
        if issue.severity == Severity.ERROR:
            typer.secho(str(issue), fg=typer.colors.RED)
        else:
            typer.secho(str(issue), fg=typer.colors.YELLOW)

    typer.echo()
    if result.is_valid:
        typer.secho(
            f"Validation passed with {result.warning_count} warning(s)",
            fg=typer.colors.YELLOW,
        )
    else:
        typer.secho(
            f"Validation failed: {result.error_count} error(s), {result.warning_count} warning(s)",
            fg=typer.colors.RED,
        )

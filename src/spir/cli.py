from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import click

from .io import read
from .io import validate as validate_api


@click.group()
def main() -> None:
    pass


@main.command()
@click.option("--src", "src_paths", multiple=True, required=True, help="Input path(s)")
@click.option(
    "--to",
    "to_format",
    type=click.Choice(["af3", "af3-server", "boltz", "chai", "protenix"]),
)
@click.option(
    "--out", "out_path", required=False, help="Output file path (or prefix for Chai)"
)
def convert(
    src_paths: Tuple[str, ...], to_format: str, out_path: Optional[str]
) -> None:
    ir = read(*src_paths)
    _dump(ir, to_format, out_path)


def _load(fmt: str, src_paths: Tuple[str, ...]):
    # Deprecated: handled by read(); kept for backward compatibility if referenced elsewhere
    ir = read(*src_paths, format=fmt)
    return ir.ci


def _dump(ir, to_format: str, out_path: Optional[str]) -> None:
    # Compute base directory and optional filename override for CLI semantics
    if to_format == "chai":
        # out_path acts as base prefix; default to "out.chai" (no extension)
        if out_path:
            base = Path(out_path)
            directory = base.parent if base.parent != Path("") else Path(".")
            filename = base.name
        else:
            directory = Path(".")
            filename = "out.chai"
        fasta_path, restraints_path = ir.write_chai(directory, filename)
        click.echo(str(fasta_path))
        if restraints_path is not None:
            click.echo(str(restraints_path))
        return

    # Single-file formats: out_path is a full filename when provided
    if out_path:
        out_p = Path(out_path)
        directory = out_p.parent if out_p.parent != Path("") else Path(".")
        filename = out_p.stem
    else:
        directory = Path(".")
        filename = None

    if to_format == "af3":
        path = ir.write_alphafold3(directory, filename)
        click.echo(str(path))
        return
    if to_format == "af3-server":
        path = ir.write_alphafoldserver(directory, filename)
        click.echo(str(path))
        return
    if to_format == "boltz":
        path = ir.write_boltz(directory, filename)
        click.echo(str(path))
        return
    if to_format == "protenix":
        path = ir.write_protenix(directory, filename)
        click.echo(str(path))
        return
    raise click.UsageError(f"Unsupported output format: {to_format}")


def _write_output(text: str, path: str) -> None:
    Path(path).write_text(text)
    click.echo(path)


@main.command()
@click.option("--src", "src_paths", multiple=True, required=True, help="Input path(s)")
@click.option(
    "--format",
    "force_format",
    type=click.Choice(["af3", "af3-server", "boltz", "chai", "protenix"]),
    required=False,
    help="Override auto-detected format",
)
@click.option("--explain/--no-explain", default=False, help="Print detailed issues")
def validate(
    src_paths: Tuple[str, ...], force_format: Optional[str], explain: bool
) -> None:
    """Validate that the given file(s) are correctly formatted."""
    msg_or_bool = validate_api(*src_paths, format=force_format, explain=True)
    valid = msg_or_bool == ""
    click.echo(f"Valid: {'yes' if valid else 'no'}")
    if explain or not valid:
        click.echo(f"Issues: {msg_or_bool or '(none)'}")

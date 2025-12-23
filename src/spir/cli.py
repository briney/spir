from typing import Optional

import typer

from spir.convert import ConvertOptions, convert
from spir.dialects import dialect_help, get_dialect
from spir.validate import validate as validate_file, print_validation_result

app = typer.Typer(no_args_is_help=True, help="SPIR: Protein folding input format converter and validator")


@app.command(name="convert")
def convert_cmd(
    in_path: str = typer.Argument(..., help="Input file path"),
    in_dialect: str = typer.Option(
        ...,
        "--from",
        help=f"Input dialect. Supported: {dialect_help()}.",
    ),
    out_prefix: str = typer.Argument(
        ..., help="Output prefix (no extension); the correct extension is added automatically."
    ),
    out_dialect: str = typer.Option(
        ...,
        "--to",
        help=(
            f"Output dialect. Supported: {dialect_help()}. "
            "Chai outputs <prefix>.fasta and <prefix>.constraints.csv."
        ),
    ),
    restraints: Optional[str] = typer.Option(
        None,
        "--restraints",
        help=(
            "Optional Chai-1 restraints CSV path (use with --from chai1; "
            "INPUT_FILE should be FASTA)."
        ),
    ),
) -> None:
    """Convert between supported dialects."""
    opts = ConvertOptions()
    convert(in_path, in_dialect, out_prefix, out_dialect, opts, restraints_path=restraints)


@app.command(name="validate")
def validate_cmd(
    input_file: str = typer.Argument(..., help="Input file to validate"),
    dialect: str = typer.Option(
        ...,
        "--dialect",
        "-d",
        help=f"Dialect to validate against. Supported: {dialect_help()}.",
    ),
    restraints: Optional[str] = typer.Option(
        None,
        "--restraints",
        help=(
            "Optional Chai-1 restraints CSV path (use with --dialect chai1; "
            "INPUT_FILE should be FASTA)."
        ),
    ),
) -> None:
    """Validate an input file against a dialect's schema and rules."""
    result = validate_file(input_file, dialect, restraints_path=restraints)
    print_validation_result(result, input_file, dialect)
    if not result.is_valid:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()

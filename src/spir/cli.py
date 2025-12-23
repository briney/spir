from typing import Optional

import typer

from spir.convert import ConvertOptions, convert
from spir.dialects import dialect_help

app = typer.Typer(no_args_is_help=True)


@app.command()
def convert_cmd(
    in_path: str = typer.Argument(...),
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


if __name__ == "__main__":
    app()

import typer

from foldir.convert import ConvertOptions, convert

app = typer.Typer(no_args_is_help=True)


@app.command()
def convert_cmd(
    in_path: str = typer.Argument(...),
    in_dialect: str = typer.Option(..., "--from"),
    out_path: str = typer.Argument(...),
    out_dialect: str = typer.Option(..., "--to"),
) -> None:
    opts = ConvertOptions()
    convert(in_path, in_dialect, out_path, out_dialect, opts)


if __name__ == "__main__":
    app()

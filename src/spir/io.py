from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple, Union

from .adapters.af3 import AF3ReaderWriter
from .adapters.af3_server import AF3ServerReaderWriter
from .adapters.boltz import BoltzReaderWriter
from .adapters.chai import ChaiReaderWriter
from .adapters.protenix import ProtenixReaderWriter
from .detect import Format, detect_from_pathnames
from .ir.ir import IR

# Kept imports above minimal; ComplexInput not required directly in this module


def read(
    *paths: str, format: Optional[Format] = None, name: Optional[str] = None
) -> IR:
    if not paths:
        raise ValueError("read() requires at least one path")
    path_objs: Tuple[Path, ...] = tuple(Path(p) for p in paths)
    for p in path_objs:
        if not p.exists():
            raise FileNotFoundError(str(p))

    fmt: Format = format or detect_from_pathnames(*[str(p) for p in path_objs])

    if fmt == "af3":
        text = path_objs[0].read_text()
        ci = AF3ReaderWriter.load_str(text)
    elif fmt == "af3-server":
        text = path_objs[0].read_text()
        ci = AF3ServerReaderWriter.load_str(text)
    elif fmt == "boltz":
        text = path_objs[0].read_text()
        ci = BoltzReaderWriter.load_str(text)
    elif fmt == "protenix":
        text = path_objs[0].read_text()
        ci = ProtenixReaderWriter.load_str(text)
    elif fmt == "chai":
        fasta = path_objs[0].read_text()
        restraints = path_objs[1].read_text() if len(path_objs) > 1 else None
        ci = ChaiReaderWriter.load_fasta_and_restraints(fasta, restraints)
    else:
        raise ValueError(f"Unsupported format: {fmt}")

    # Set name if provided or missing
    if name is not None:
        ci.name = name
    elif not ci.name:
        ci.name = path_objs[0].stem

    return IR(ci=ci, source_paths=path_objs, source_format=fmt)


def validate(
    *paths: str, format: Optional[Format] = None, explain: bool = False
) -> Union[bool, str]:
    """Validate that the given path(s) conform to a supported input format.

    When explain=False (default): returns True if valid, False otherwise.
    When explain=True: returns an empty string if valid, otherwise a human-readable
    description of the problem.
    """
    try:
        if not paths:
            raise ValueError("validate() requires at least one path")
        path_objs: Tuple[Path, ...] = tuple(Path(p) for p in paths)
        for p in path_objs:
            if not p.exists():
                raise FileNotFoundError(str(p))

        fmt: Format = format or detect_from_pathnames(*[str(p) for p in path_objs])

        # Attempt to parse using the appropriate adapter; rely on pydantic for
        # structural validation where applicable.
        if fmt == "af3":
            text = path_objs[0].read_text()
            _ = AF3ReaderWriter.load_str(text)
        elif fmt == "af3-server":
            text = path_objs[0].read_text()
            _ = AF3ServerReaderWriter.load_str(text)
        elif fmt == "boltz":
            text = path_objs[0].read_text()
            _ = BoltzReaderWriter.load_str(text)
        elif fmt == "protenix":
            text = path_objs[0].read_text()
            _ = ProtenixReaderWriter.load_str(text)
        elif fmt == "chai":
            fasta = path_objs[0].read_text()
            restraints = path_objs[1].read_text() if len(path_objs) > 1 else None
            _ = ChaiReaderWriter.load_fasta_and_restraints(fasta, restraints)
        else:
            raise ValueError(f"Unsupported format: {fmt}")

        return "" if explain else True
    except Exception as e:
        msg = f"{type(e).__name__}: {e}"
        return msg if explain else False

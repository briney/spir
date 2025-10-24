from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional, Tuple, Union

from ..adapters.af3 import AF3ReaderWriter
from ..adapters.af3_server import AF3ServerReaderWriter
from ..adapters.boltz import BoltzReaderWriter
from ..adapters.chai import ChaiReaderWriter
from ..adapters.protenix import ProtenixReaderWriter
from ..constraints.normalize import normalize_multi_ccd_for_boltz
from ..detect import Format
from .model import ComplexInput

PathLike = Union[str, Path]


class IR:
    def __init__(
        self,
        ci: ComplexInput,
        source_paths: Tuple[Path, ...],
        source_format: Format,
    ) -> None:
        self.ci = ci
        self.source_paths = source_paths
        self.source_format = source_format

    @property
    def stem(self) -> str:
        if self.source_paths and len(self.source_paths) > 0:
            first = self.source_paths[0]
            try:
                return first.stem
            except Exception:
                pass
        if self.ci.name:
            return self.ci.name
        return "spir-job"

    # ---- Writers (explicit) ----
    def write_alphafold3(
        self, directory: PathLike, filename: Optional[str] = None
    ) -> Path:
        out_dir = _ensure_dir(directory)
        name = filename or self.stem
        path = out_dir / f"{name}.af3.json"
        text = AF3ReaderWriter.dump(self.ci)
        path.write_text(text)
        return path

    def write_boltz(self, directory: PathLike, filename: Optional[str] = None) -> Path:
        out_dir = _ensure_dir(directory)
        name = filename or self.stem
        path = out_dir / f"{name}.boltz.yaml"
        # Normalize multi-CCD glycans into per-component ligands for Boltz constraints
        ci_norm = normalize_multi_ccd_for_boltz(self.ci)
        text = BoltzReaderWriter.dump(ci_norm)
        path.write_text(text)
        return path

    def write_protenix(
        self, directory: PathLike, filename: Optional[str] = None
    ) -> Path:
        out_dir = _ensure_dir(directory)
        name = filename or self.stem
        path = out_dir / f"{name}.protenix.json"
        text = ProtenixReaderWriter.dump(self.ci)
        path.write_text(text)
        return path

    def write_chai(
        self, directory: PathLike, filename: Optional[str] = None
    ) -> Tuple[Path, Optional[Path]]:
        out_dir = _ensure_dir(directory)
        name = filename or self.stem
        fasta_text, restraints_text = ChaiReaderWriter.dump(self.ci)
        fasta_path = out_dir / f"{name}.fasta"
        fasta_path.write_text(fasta_text)
        restraints_path: Optional[Path] = None
        if restraints_text:
            restraints_path = out_dir / f"{name}.restraints.csv"
            restraints_path.write_text(restraints_text)
        return fasta_path, restraints_path

    def write_alphafoldserver(
        self, directory: PathLike, filename: Optional[str] = None
    ) -> Path:
        out_dir = _ensure_dir(directory)
        name = filename or self.stem
        path = out_dir / f"{name}.af3server.json"
        text = AF3ServerReaderWriter.dump(self.ci)
        path.write_text(text)
        return path

    # ---- Convenience dispatcher ----
    def write(
        self,
        format: Literal["af3", "af3-server", "boltz", "chai", "protenix"],
        directory: PathLike,
        filename: Optional[str] = None,
    ):
        if format == "af3":
            return self.write_alphafold3(directory, filename)
        if format == "af3-server":
            return self.write_alphafoldserver(directory, filename)
        if format == "boltz":
            return self.write_boltz(directory, filename)
        if format == "chai":
            return self.write_chai(directory, filename)
        if format == "protenix":
            return self.write_protenix(directory, filename)
        raise ValueError(f"Unsupported format: {format}")


def _ensure_dir(path: PathLike) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p

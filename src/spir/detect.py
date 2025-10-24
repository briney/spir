from __future__ import annotations

import json
from typing import Literal, Tuple

import yaml

Format = Literal["af3", "af3-server", "boltz", "chai", "protenix"]


def detect_from_pathnames(*paths: str) -> Format:
    # Simple heuristics by extension and content
    if len(paths) == 1:
        p = paths[0]
        low = p.lower()
        if low.endswith((".yaml", ".yml")):
            # Boltz or others; sniff keys
            with open(p, "r") as f:
                try:
                    y = yaml.safe_load(f)
                except Exception:
                    raise
            if isinstance(y, dict) and "sequences" in y:
                return "boltz"
        if low.endswith(".json"):
            with open(p, "r") as f:
                j = json.load(f)
            if isinstance(j, dict) and j.get("dialect") == "alphafold3":
                return "af3"
            if isinstance(j, list):
                # AF3 server vs Protenix; AF3 server entries contain dialect in job dict
                if (
                    j
                    and isinstance(j[0], dict)
                    and j[0].get("dialect") == "alphafoldserver"
                ):
                    return "af3-server"
                # Protenix uses [ { name, sequences, covalent_bonds? } ]
                if j and isinstance(j[0], dict) and "sequences" in j[0]:
                    return "protenix"
        if low.endswith((".fasta", ".fa")):
            return "chai"
    # Multiple paths -> Chai if FASTA + CSV present
    exts = {p.lower().rsplit(".", 1)[-1] for p in paths if "." in p}
    if "fasta" in exts or "fa" in exts:
        return "chai"
    raise ValueError("Unable to auto-detect format from given paths")

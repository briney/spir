from __future__ import annotations

from typing import Dict

from foldir.dialects.alphafold3 import AlphaFold3Dialect
from foldir.dialects.alphafold3_server import AlphaFold3ServerDialect
from foldir.dialects.boltz2 import Boltz2Dialect
from foldir.dialects.chai1 import Chai1Dialect
from foldir.dialects.protenix import ProtenixDialect


_DIALECTS: Dict[str, object] = {
    "alphafold3": AlphaFold3Dialect(),
    "alphafold3server": AlphaFold3ServerDialect(),
    "alphafoldserver": AlphaFold3ServerDialect(),
    "boltz2": Boltz2Dialect(),
    "chai1": Chai1Dialect(),
    "protenix": ProtenixDialect(),
}


def get_dialect(name: str):
    key = name.lower()
    if key not in _DIALECTS:
        raise ValueError(f"Unknown dialect: {name}")
    return _DIALECTS[key]

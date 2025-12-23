from __future__ import annotations

from typing import Any

import yaml


def read_yaml(path: str) -> Any:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def write_yaml(path: str, payload: Any) -> None:
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False)

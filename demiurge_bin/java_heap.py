from __future__ import annotations

import re
from typing import Any


JAVA_HEAP_PATTERN = re.compile(r"[1-9][0-9]*[MG]", re.IGNORECASE)
DEFAULT_JAVA_HEAP = "4G"
_UNIT_BYTES = {"M": 1024**2, "G": 1024**3}


def normalize_java_heap(value: Any) -> str:
    """Return a JVM-safe, canonical maximum heap such as ``512M`` or ``4G``."""
    if not isinstance(value, str) or JAVA_HEAP_PATTERN.fullmatch(value) is None:
        raise ValueError(
            "java heap must be a positive integer followed by M or G "
            "(for example 512M, 2G, 4G, or 8G)"
        )
    return value.upper()


def java_heap_bytes(value: Any) -> int:
    normalized = normalize_java_heap(value)
    return int(normalized[:-1]) * _UNIT_BYTES[normalized[-1]]

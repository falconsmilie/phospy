from __future__ import annotations

import pandas as pd

_PANDAS_MAJOR_VERSION_TOKEN = str(pd.__version__).split(".", maxsplit=1)[0]
try:
    _PANDAS_MAJOR_VERSION = int(_PANDAS_MAJOR_VERSION_TOKEN)
except ValueError:
    _PANDAS_MAJOR_VERSION = 0

_COW_DETACH_AVAILABLE = _PANDAS_MAJOR_VERSION >= 3


def detached_frame_copy(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a detached frame copy at public ownership boundaries.

    On pandas 3+, Copy-on-Write guarantees mutation isolation for shallow
    frame copies, so we avoid eager deep copies. Older pandas versions still
    require deep copies for the same semantics.
    """

    return frame.copy(deep=not _COW_DETACH_AVAILABLE)


def detached_series_copy(series: pd.Series) -> pd.Series:
    """Return a detached Series copy at public ownership boundaries."""

    return series.copy(deep=not _COW_DETACH_AVAILABLE)


__all__ = [
    "detached_frame_copy",
    "detached_series_copy",
]

"""Dataset builder contracts."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

DatasetInput = pd.DataFrame | str | Path

__all__ = ["DatasetInput"]

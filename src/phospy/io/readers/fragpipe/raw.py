"""Raw FragPipe table loading and source-shape checks."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

import pandas as pd

from phospy.io.readers._table_parsing import require_non_empty_unique_columns
from phospy.io.readers.importers import _read_upstream_table


def read_fragpipe_source(source: object) -> pd.DataFrame:
    """Read a FragPipe source object into a table without domain conversion."""

    return _read_upstream_table(source)


def require_fragpipe_source_columns(source: pd.DataFrame) -> None:
    require_non_empty_unique_columns(source, importer_label="FragPipe")


__all__ = ["read_fragpipe_source", "require_fragpipe_source_columns"]

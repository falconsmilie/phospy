"""Reader ownership for tabular I/O primitives."""

from phospy.io.readers.tables import (
    read_table,
    table_format_from_path,
    table_suffix_for_format,
    write_table,
)

__all__ = [
    "read_table",
    "table_format_from_path",
    "table_suffix_for_format",
    "write_table",
]

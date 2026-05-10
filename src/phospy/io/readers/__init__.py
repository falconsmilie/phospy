"""Reader ownership for tabular I/O primitives."""

from phospy.io.readers.tables import (
    read_contrast_matrix,
    read_design_matrix,
    read_phospho_matrix,
    read_sample_metadata,
    read_site_metadata,
    read_table,
    read_total_matrix,
    table_format_from_path,
    table_suffix_for_format,
    write_table,
)

__all__ = [
    "read_contrast_matrix",
    "read_design_matrix",
    "read_phospho_matrix",
    "read_sample_metadata",
    "read_site_metadata",
    "read_table",
    "read_total_matrix",
    "table_format_from_path",
    "table_suffix_for_format",
    "write_table",
]

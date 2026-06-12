"""Reader ownership for tabular I/O primitives."""

from __future__ import annotations

from typing import TYPE_CHECKING

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

if TYPE_CHECKING:
    from phospy.io.readers.importers import (
        ColumnMappedPhosphositeImporter,
        MappedPhosphositeTableImporter,
    )

__all__ = [
    "ColumnMappedPhosphositeImporter",
    "MappedPhosphositeTableImporter",
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


def __getattr__(name: str) -> object:
    if name == "ColumnMappedPhosphositeImporter":
        from phospy.io.readers.importers import ColumnMappedPhosphositeImporter

        return ColumnMappedPhosphositeImporter
    if name == "MappedPhosphositeTableImporter":
        from phospy.io.readers.importers import MappedPhosphositeTableImporter

        return MappedPhosphositeTableImporter
    raise AttributeError(name)

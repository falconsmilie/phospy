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
    from phospy.io.readers.maxquant import (
        MaxQuantColumnMapping,
        MaxQuantPhosphositeImporter,
        MaxQuantPhosphositeImportRequest,
    )

__all__ = [
    "ColumnMappedPhosphositeImporter",
    "MappedPhosphositeTableImporter",
    "MaxQuantColumnMapping",
    "MaxQuantPhosphositeImporter",
    "MaxQuantPhosphositeImportRequest",
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
    if name == "MaxQuantColumnMapping":
        from phospy.io.readers.maxquant import MaxQuantColumnMapping

        return MaxQuantColumnMapping
    if name == "MaxQuantPhosphositeImporter":
        from phospy.io.readers.maxquant import MaxQuantPhosphositeImporter

        return MaxQuantPhosphositeImporter
    if name == "MaxQuantPhosphositeImportRequest":
        from phospy.io.readers.maxquant import MaxQuantPhosphositeImportRequest

        return MaxQuantPhosphositeImportRequest
    raise AttributeError(name)

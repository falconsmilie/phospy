"""Structured data access and output persistence domain.

This package owns shared table reading, mapping-file loading, and path-oriented
I/O helpers that are not specific scientific behaviours. Domain logic should
not accumulate here.
"""

from .mappings import load_grouped_mapping, load_string_mapping
from .publishing import OutputPublisher, RunManifestWriter, package_version
from .readers import (
    DEFAULT_TEXT_ENCODING,
    clean_columns,
    clean_table_columns,
    default_text_encoding,
    infer_text_encoding,
    load_phospho_table,
    load_pred_mat,
    load_total_table,
    read_table,
    read_table_raw,
)
from .writers import (
    CoreOutputArtifact,
    CoreOutputFormat,
    CoreOutputWriter,
    CoreProcessingResultWriter,
    DelimitedTabularWriter,
    KinaseActivityResultWriter,
    KinaseActivityWriter,
    ParquetTabularWriter,
    TabularOutputWriter,
)

__all__ = [
    "DEFAULT_TEXT_ENCODING",
    "CoreOutputArtifact",
    "CoreOutputFormat",
    "CoreOutputWriter",
    "CoreProcessingResultWriter",
    "DelimitedTabularWriter",
    "KinaseActivityResultWriter",
    "KinaseActivityWriter",
    "OutputPublisher",
    "ParquetTabularWriter",
    "RunManifestWriter",
    "TabularOutputWriter",
    "clean_columns",
    "clean_table_columns",
    "default_text_encoding",
    "infer_text_encoding",
    "load_grouped_mapping",
    "load_phospho_table",
    "load_pred_mat",
    "load_string_mapping",
    "load_total_table",
    "package_version",
    "read_table",
    "read_table_raw",
]

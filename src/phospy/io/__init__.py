"""Structured data access and output persistence domain.

This package owns shared table reading, mapping-file loading, and path-oriented
I/O helpers that are not specific scientific behaviours. Domain logic should
not accumulate here.
"""

from .mappings import load_grouped_mapping, load_string_mapping
from .publishing import (
    OutputPublisher,
    RunManifestWriter,
    load_simple_kinase_workflow_output_bundle,
    load_simple_kinase_workflow_output_bundle_metadata,
    package_version,
    save_simple_kinase_workflow_output_bundle,
)
from .readers import (
    DEFAULT_TEXT_ENCODING,
    SimpleKinaseWorkflowBundleReader,
    SimpleKinaseWorkflowOutputBundle,
    clean_columns,
    clean_table_columns,
    default_text_encoding,
    load_phospho_table,
    load_pred_mat,
    load_total_table,
    read_table,
    read_table_raw,
)
from .writers import (
    BundleTableArtifact,
    CoreOutputArtifact,
    CoreOutputFormat,
    CoreOutputWriter,
    CoreProcessingResultWriter,
    DelimitedTabularWriter,
    KinaseActivityResultWriter,
    KinaseActivityWriter,
    ParquetTabularWriter,
    SimpleKinaseWorkflowBundleWriter,
    TabularOutputWriter,
)

__all__ = [
    "DEFAULT_TEXT_ENCODING",
    "BundleTableArtifact",
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
    "SimpleKinaseWorkflowBundleReader",
    "SimpleKinaseWorkflowBundleWriter",
    "SimpleKinaseWorkflowOutputBundle",
    "TabularOutputWriter",
    "clean_columns",
    "clean_table_columns",
    "default_text_encoding",
    "load_grouped_mapping",
    "load_phospho_table",
    "load_pred_mat",
    "load_simple_kinase_workflow_output_bundle",
    "load_simple_kinase_workflow_output_bundle_metadata",
    "load_string_mapping",
    "load_total_table",
    "package_version",
    "read_table",
    "read_table_raw",
    "save_simple_kinase_workflow_output_bundle",
]

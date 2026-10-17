"""Logical-owner policy support governed by ADR-0053 and ADR-0054."""

from __future__ import annotations

import csv
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

from import_graph_support import ImportGraph, ImportRecord

OWNERS = frozenset(
    {
        "activities",
        "datasets",
        "preprocessing",
        "preparation_records",
        "differential",
        "enrichment",
        "kinase",
        "signalome",
        "phosphosites",
        "taxonomy",
        "quantitative",
        "evidence",
        "references",
        "statistics",
        "provenance",
        "run_lineage",
        "errors",
        "frames",
        "policies",
        "validation",
        "io",
        "public_facades",
    }
)

CAPABILITY_OWNERS = frozenset(
    {
        "activities",
        "datasets",
        "preprocessing",
        "differential",
        "enrichment",
        "kinase",
        "signalome",
    }
)
SHARED_OWNERS = frozenset(
    {
        "phosphosites",
        "taxonomy",
        "quantitative",
        "evidence",
        "references",
        "statistics",
        "provenance",
    }
)
PASSIVE_RECORD_OWNERS = frozenset({"preparation_records"})
FOUNDATION_OWNERS = frozenset({"errors", "frames", "policies", "validation"})


@dataclass(frozen=True, slots=True)
class MixedModuleReview:
    """Exhaustive ownership review for one provisionally single-owner module."""

    primary_owner: str
    declarations: frozenset[str]
    foreign_owners: Mapping[str, str]


SHARED_ALLOWED_TARGETS: Mapping[str, frozenset[str]] = {
    "errors": frozenset(),
    "taxonomy": frozenset({"errors"}),
    "policies": frozenset({"errors"}),
    "frames": frozenset({"errors"}),
    "validation": frozenset({"errors", "frames", "policies"}),
    "provenance": frozenset({"errors", "frames", "taxonomy"}),
    "phosphosites": frozenset(
        {"errors", "frames", "policies", "taxonomy", "validation"}
    ),
    "quantitative": frozenset(
        {"errors", "frames", "policies", "provenance", "validation"}
    ),
    "statistics": frozenset({"errors", "frames", "policies", "validation"}),
    "references": frozenset(
        {
            "errors",
            "frames",
            "phosphosites",
            "policies",
            "provenance",
            "taxonomy",
            "validation",
        }
    ),
    "evidence": frozenset(
        {
            "errors",
            "frames",
            "phosphosites",
            "policies",
            "provenance",
            "quantitative",
            "taxonomy",
            "validation",
        }
    ),
    "preparation_records": frozenset(
        {"errors", "frames", "policies", "provenance", "quantitative", "validation"}
    ),
}

# Closed current symbol-owner map for passive preparation records. The live
# module map classifies each source file by its direct physical owner.
PREPARATION_RECORD_ALLOWED_SYMBOLS: Mapping[str, frozenset[str]] = {
    "phospy.preparation_records.policies": frozenset(
        {
            "MissingDataPolicy",
            "NormalisationPolicy",
            "SiteSequenceConflictPolicy",
            "SiteSequenceResolutionMode",
            "TotalProteinCorrectionIdentityMatchingPolicy",
            "TotalProteinCorrectionPolicy",
        }
    ),
    "phospy.preparation_records.imputation": frozenset(
        {
            "DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICIES",
            "DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_COLUMN_MEAN_WITH_CAVEAT",
            "DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_ERROR",
            "DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_VERSION",
            "IMPUTATION_INPUT_SCALE_SOURCE_CALLER_SELECTED",
            "IMPUTATION_INPUT_SCALE_SOURCE_METHOD_REQUIRED",
            "IMPUTATION_OPERATION_ORDERS",
            "IMPUTATION_OPERATION_ORDER_AFTER_INTENSITY_TRANSFORM",
            "IMPUTATION_OPERATION_ORDER_BEFORE_INTENSITY_TRANSFORM",
            "IMPUTATION_OPERATION_ORDER_NO_INTENSITY_TRANSFORM",
        }
    ),
    "phospy.preparation_records.report_schema": frozenset(
        {
            "COMPARISON_GROUP_STATS_COLUMNS",
            "COMPARISON_PAIR_STATS_COLUMNS",
            "DUPLICATE_SITE_RESOLUTION_COLUMNS",
            "METADATA_CONFLICT_COLUMNS",
            "OPERATIONS_COLUMNS",
            "ROW_AUDIT_COLUMNS",
            "ROW_COUNTS_COLUMNS",
            "ComparisonGroupStatsRow",
            "ComparisonPairStatsRow",
            "DuplicateSiteResolutionRow",
            "MetadataConflictRow",
            "PreprocessingOperationRow",
            "PreprocessingRowAuditRow",
            "PreprocessingRowCountRow",
            "_DataclassInstance",
            "_ReportRowT",
            "_columns_for_dataclass",
            "comparison_group_stats_rows_from_dataframe",
            "comparison_pair_stats_rows_from_dataframe",
            "dataframe_from_comparison_group_stats_rows",
            "dataframe_from_comparison_pair_stats_rows",
            "dataframe_from_duplicate_site_resolution_rows",
            "dataframe_from_metadata_conflict_rows",
            "dataframe_from_operation_rows",
            "dataframe_from_row_audit_rows",
            "dataframe_from_row_count_rows",
            "dataframe_from_rows",
            "duplicate_site_resolution_rows_from_dataframe",
            "metadata_conflict_rows_from_dataframe",
            "missing_columns",
            "operation_rows_from_dataframe",
            "reorder_columns",
            "row_audit_rows_from_dataframe",
            "row_count_rows_from_dataframe",
            "rows_from_dataframe",
        }
    ),
    "phospy.preparation_records.batch_correction": frozenset(
        {
            "BATCH_CORRECTION_CONFOUNDING_CONFOUNDED",
            "BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE",
            "BATCH_CORRECTION_CONFOUNDING_NOT_CHECKED",
            "BATCH_CORRECTION_CONFOUNDING_PASSED",
            "BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS",
            "BATCH_CORRECTION_STATUS_APPLIED",
            "BATCH_CORRECTION_STATUS_DISABLED",
            "BATCH_CORRECTION_STATUS_REJECTED",
            "BatchCorrectionConfoundingCheckStatus",
            "BatchCorrectionDiagnostics",
            "BatchCorrectionPolicy",
            "BatchCorrectionReport",
            "BatchCorrectionStatus",
            "MatrixShape",
            "_normalize_condition_columns",
            "_optional_text",
            "_shape_to_payload",
        }
    ),
    "phospy.preparation_records.protein_aware": frozenset(
        {
            "DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES",
            "DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_ALLOW_MISSING_WITH_REPORT",
            "DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS",
            "DATASET_PROTEIN_AWARE_PREPARATION_POLICIES",
            "DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED",
            "DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS",
            "DatasetProteinAwarePreparationMappingPolicy",
            "DatasetProteinAwarePreparationPolicy",
            "PROTEIN_AWARE_AMBIGUOUS_MAPPING_DIAGNOSTIC_COLUMNS",
            "PROTEIN_AWARE_MATCHED_PAIR_COLUMNS",
            "PROTEIN_AWARE_MISSING_PROTEIN_ABUNDANCE_DIAGNOSTIC_COLUMNS",
            "PROTEIN_AWARE_SITE_ELIGIBILITY_COLUMNS",
            "ProteinAwareMappingDiagnostics",
            "ProteinAwarePreparationEligibility",
            "ProteinAwareSampleAlignmentDiagnostics",
            "ProteinAwareSiteEligibility",
            "ProteinAwareTransformationStateDiagnostics",
            "ProteinMappingStatus",
            "_dataframe_records_payload",
            "_normalize_string_tuple",
            "_optional_non_empty_string",
            "_owned_diagnostic_table",
            "_payload_value",
            "_require_non_empty_string",
        }
    ),
    "phospy.preparation_records.group_evidence": frozenset(
        {
            "_canonical_mask",
            "hash_group_aware_imputation_mask",
            "hash_group_aware_mechanism_mask",
        }
    ),
    "phospy.preparation_records.validation.sample_groups": frozenset(
        {
            "ResolvedSampleGroups",
            "SampleGroupMetadataResolver",
            "_format_label_preview",
            "_format_positions",
            "_group_samples",
            "_is_missing_value",
            "_normalize_label_index",
            "_require_column",
            "_require_no_missing_or_extra_metadata_rows",
            "_require_unique_labels",
            "_resolve_groups_by_sample",
        }
    ),
}

# Public and compatibility facades may import only these reviewed,
# owner-defined contract/export modules and sanctioned composition entrypoints.
# Approval is exact: a sibling or descendant implementation module does not
# inherit permission merely because its name avoids an implementation keyword.
PUBLIC_FACADE_ALLOWED_TARGETS: Mapping[str, frozenset[str]] = {
    "datasets": frozenset(
        {
            "phospy.datasets.builder",
            "phospy.datasets.contracts",
            "phospy.datasets.models",
            "phospy.datasets.schemas",
        }
    ),
    "differential": frozenset(
        {
            "phospy.differential.config",
            "phospy.differential.contracts",
            "phospy.differential.design.contrast_helpers",
            "phospy.differential.design.models",
            "phospy.differential.diagnostics",
            "phospy.differential.moderation.config",
            "phospy.differential.policies",
            "phospy.differential.public",
            "phospy.differential.results",
            "phospy.differential.tables",
        }
    ),
    "enrichment": frozenset(
        {
            "phospy.enrichment.contracts.collections",
            "phospy.enrichment.contracts.config",
            "phospy.enrichment.contracts.identifier_sets",
            "phospy.enrichment.contracts.requests",
            "phospy.enrichment.contracts.results",
            "phospy.enrichment.contracts.vocabulary",
            "phospy.enrichment.workflow",
        }
    ),
    "errors": frozenset(
        {
            "phospy.errors",
            "phospy.errors.base",
            "phospy.errors.build",
            "phospy.errors.input",
            "phospy.errors.provenance",
            "phospy.errors.references",
            "phospy.errors.transformations",
            "phospy.errors.validation",
            "phospy.errors.workflows",
        }
    ),
    "evidence": frozenset(
        {
            "phospy.evidence.config",
        }
    ),
    "frames": frozenset({"phospy.frames.table_schema"}),
    "io": frozenset(
        {
            "phospy.io.importers.contracts",
            "phospy.io.bundles.kinase_library",
            "phospy.references.acquisition",
            "phospy.io.publishers.workflows",
            "phospy.io.readers.dataset_inputs",
        }
    ),
    "kinase": frozenset(
        {
            "phospy.kinase.activity.config",
            "phospy.kinase.activity.diagnostics",
            "phospy.kinase.activity.results",
            "phospy.kinase.results",
            "phospy.kinase.requests",
            "phospy.kinase.result_boundary",
            "phospy.kinase.policies",
            "phospy.kinase.prediction.config",
            "phospy.kinase.prediction.results",
            "phospy.kinase.references.models",
            "phospy.kinase.scoring.config",
            "phospy.kinase.scoring.policies",
            "phospy.kinase.scoring.results",
            "phospy.kinase.activity.tables",
            "phospy.kinase.scoring.tables",
            "phospy.kinase.scoring.sequence_validation",
            "phospy.kinase.public",
        }
    ),
    "preprocessing": frozenset(
        {
            "phospy.preprocessing.contracts.configuration",
            "phospy.preprocessing.contracts.configuration.dataset",
            "phospy.preprocessing.control_sites.contracts",
            "phospy.preprocessing.control_sites.discovery_contracts",
            "phospy.preprocessing.control_sites.discovery_provenance",
            "phospy.preprocessing.control_sites.discovery_results",
            "phospy.preprocessing.control_sites.workflow",
            "phospy.preprocessing.batch_correction.contracts.configuration",
            "phospy.preprocessing.batch_correction.contracts.missingness",
            "phospy.preprocessing.batch_correction.contracts.stage",
        }
    ),
    "provenance": frozenset({"phospy.provenance.result_caveats"}),
    "quantitative": frozenset(
        {
            "phospy.quantitative.evidence",
            "phospy.quantitative.models",
        }
    ),
    "references": frozenset(
        {
            "phospy.references.builder",
            "phospy.references.config",
            "phospy.references.models",
            "phospy.references.tables",
        }
    ),
    "signalome": frozenset(
        {
            "phospy.signalome.config",
            "phospy.signalome.errors",
            "phospy.signalome.public",
            "phospy.signalome.requests",
            "phospy.signalome.results",
            "phospy.signalome.tables",
        }
    ),
    "statistics": frozenset({"phospy.statistics.multiple_testing"}),
    "taxonomy": frozenset({"phospy.taxonomy.organisms"}),
}

PUBLIC_FACADE_STEADY_STATE_TARGET_OWNERS = frozenset(
    {
        "datasets",
        "differential",
        "enrichment",
        "errors",
        "evidence",
        "io",
        "kinase",
        "preprocessing",
        "provenance",
        "quantitative",
        "references",
        "signalome",
        "statistics",
        "taxonomy",
    }
)

# IO adapters may depend inward only on these exact, reviewed passive contract,
# model, provenance, and serialization modules.  Exact membership is
# intentional: siblings and descendants are implementation paths unless they
# receive their own review and entry here.
IO_ADAPTER_ALLOWED_TARGETS: Mapping[str, frozenset[str]] = {
    "activities": frozenset({"phospy.activities.semantics"}),
    "datasets": frozenset(
        {
            "phospy.datasets.contracts",
            "phospy.datasets.trusted_assertions",
            "phospy.datasets.models",
            "phospy.datasets.processing_state",
        }
    ),
    "enrichment": frozenset(
        {
            "phospy.enrichment.contracts.collections",
            "phospy.enrichment.contracts.vocabulary",
        }
    ),
    "evidence": frozenset(
        {
            "phospy.evidence.config",
            "phospy.evidence.contracts",
            "phospy.evidence.localisation",
            "phospy.evidence.modified_peptides",
            "phospy.evidence.multi_site",
            "phospy.evidence.quantification_depth",
        }
    ),
    "kinase": frozenset(
        {
            "phospy.kinase.references.loading",
            "phospy.kinase.references.models",
            "phospy.kinase.references.projection_contracts",
            "phospy.kinase.requests",
            "phospy.kinase.results",
            "phospy.kinase.scoring.config",
            "phospy.kinase.scoring.policies",
            "phospy.kinase.scoring.results",
            "phospy.kinase.prediction.config",
            "phospy.kinase.prediction.results",
            "phospy.kinase.prediction.tables",
            "phospy.kinase.activity.config",
            "phospy.kinase.activity.diagnostics",
            "phospy.kinase.activity.membership",
            "phospy.kinase.activity.models",
            "phospy.kinase.activity.results",
        }
    ),
    "phosphosites": frozenset({"phospy.phosphosites.identifiers"}),
    "preparation_records": frozenset({"phospy.preparation_records.policies"}),
    "preprocessing": frozenset(
        {
            "phospy.preprocessing.contracts",
            "phospy.preprocessing.contracts.configuration.dataset",
            "phospy.preprocessing.contracts.policies",
            "phospy.preprocessing.results",
            "phospy.preprocessing.serialization",
        }
    ),
    "run_lineage": frozenset(
        {
            "phospy.run_lineage.models",
            "phospy.run_lineage.serialization",
        }
    ),
    "provenance": frozenset(
        {
            "phospy.provenance.hashing",
            "phospy.provenance.immutability",
            "phospy.provenance.models",
            "phospy.provenance.result_caveats",
            "phospy.provenance.serialization",
        }
    ),
    "quantitative": frozenset(
        {
            "phospy.quantitative.assertions",
            "phospy.quantitative.models",
            "phospy.quantitative.serialization",
        }
    ),
    "references": frozenset(
        {
            "phospy.references.acquisition",
            "phospy.references.config",
            "phospy.references.models",
        }
    ),
    "signalome": frozenset(
        {
            "phospy.signalome.config",
            "phospy.signalome.constants",
            "phospy.signalome.context",
            "phospy.signalome.models",
            "phospy.signalome.requests",
            "phospy.signalome.result_validation",
            "phospy.signalome.results",
            "phospy.signalome.tables",
            "phospy.signalome.tables.assignments",
            "phospy.signalome.tables.context",
            "phospy.signalome.tables.modules",
            "phospy.signalome.tables.network",
        }
    ),
    "taxonomy": frozenset({"phospy.taxonomy.organisms"}),
}

IO_ADAPTER_STEADY_STATE_TARGET_OWNERS = frozenset(
    {
        "activities",
        "datasets",
        "enrichment",
        "evidence",
        "errors",
        "frames",
        "kinase",
        "phosphosites",
        "policies",
        "preparation_records",
        "preprocessing",
        "provenance",
        "run_lineage",
        "quantitative",
        "references",
        "signalome",
        "taxonomy",
        "validation",
    }
)

# Every top-level declaration in every reviewed mixed module is registered.
# Declarations absent from ``foreign_owners`` belong to the provisional primary
# owner.  This makes a newly added declaration fail review instead of silently
# inheriting the file owner.
REVIEWED_MIXED_MODULES: Mapping[str, MixedModuleReview] = {
    "src/phospy/errors/workflows.py": MixedModuleReview(
        primary_owner="errors",
        declarations=frozenset(
            {
                "class:PhosPyWorkflowError",
                "class:WorkflowBoundaryError",
                "class:WorkflowStageError",
            }
        ),
        foreign_owners={},
    ),
}

# This is the permanent, reviewed logical-owner graph. Path-level helpers below
# further narrow the capability, IO, and facade edges. Temporary architecture
# exceptions are deliberately absent: they ratchet the current tree without
# becoming permitted owner directions.
PERMITTED_LOGICAL_OWNER_EDGES = frozenset(
    {
        (source, target)
        for source, targets in SHARED_ALLOWED_TARGETS.items()
        for target in targets
    }
    | {
        (source, target)
        for source in CAPABILITY_OWNERS
        for target in SHARED_OWNERS | FOUNDATION_OWNERS
    }
    | {
        ("datasets", "preprocessing"),
        ("datasets", "preparation_records"),
        ("preprocessing", "preparation_records"),
        ("differential", "datasets"),
        ("differential", "preprocessing"),
        ("differential", "preparation_records"),
        ("kinase", "datasets"),
        ("kinase", "activities"),
        ("signalome", "kinase"),
        ("signalome", "datasets"),
        ("datasets", "run_lineage"),
        ("differential", "run_lineage"),
        ("enrichment", "run_lineage"),
        ("kinase", "run_lineage"),
        ("signalome", "run_lineage"),
        ("run_lineage", "errors"),
        ("run_lineage", "evidence"),
        ("run_lineage", "frames"),
        ("run_lineage", "preprocessing"),
        ("run_lineage", "provenance"),
        ("run_lineage", "references"),
        ("run_lineage", "taxonomy"),
    }
    | {("io", target) for target in IO_ADAPTER_STEADY_STATE_TARGET_OWNERS}
    | {
        ("public_facades", target)
        for target in PUBLIC_FACADE_STEADY_STATE_TARGET_OWNERS
    }
)

DATASET_PATH_READER_EDGE = (
    "phospy.datasets.builder",
    "phospy.io.readers.dataset_inputs",
)

DATASET_PREPROCESSING_ALLOWED_IMPORTS: Mapping[tuple[str, str], frozenset[str]] = {
    (
        "phospy.datasets.builder",
        "phospy.preprocessing.dataset_construction",
    ): frozenset(
        {
            "SpsRuvStyleBatchCorrectionRunner",
            "create_default_dataset_preprocessing_runner",
        }
    ),
    (
        "phospy.datasets.contracts",
        "phospy.preprocessing.contracts.configuration.dataset",
    ): frozenset({"DatasetPreprocessingConfig"}),
    (
        "phospy.datasets.validation.request",
        "phospy.preprocessing.contracts.configuration.comparisons",
    ): frozenset({"DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS"}),
    (
        "phospy.datasets.validation.request",
        "phospy.preprocessing.contracts.configuration.total_protein",
    ): frozenset({"DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE"}),
    (
        "phospy.datasets.validation.request",
        "phospy.preprocessing.dataset_construction",
    ): frozenset({"validate_corrected_preprocessing_output"}),
    (
        "phospy.datasets.validation.request",
        "phospy.preprocessing.workflow",
    ): frozenset({"DatasetPreprocessingConfigValidator"}),
    (
        "phospy.datasets.construction.contracts",
        "phospy.preprocessing.contracts",
    ): frozenset({"PreprocessingPlan", "PreprocessingStageExecution"}),
    (
        "phospy.datasets.construction.contracts",
        "phospy.preprocessing.contracts.configuration.dataset",
    ): frozenset({"DatasetPreprocessingConfig"}),
    (
        "phospy.datasets.construction.contracts",
        "phospy.preprocessing.dataset_construction",
    ): frozenset({"CorrectedPreprocessingOutput", "ResolvedBatchCorrectionMetadata"}),
    (
        "phospy.datasets.construction.executor",
        "phospy.preprocessing.contracts",
    ): frozenset({"reject_external_corrected_output_after_downstream_preprocessing"}),
    (
        "phospy.datasets.construction.executor",
        "phospy.preprocessing.workflow",
    ): frozenset({"AdditivePreprocessingScaleGuard"}),
    (
        "phospy.datasets.construction.input_resolution",
        "phospy.preprocessing.contracts.policies",
    ): frozenset({"SiteMatrixDuplicateSitePolicy", "SiteMatrixPolicy"}),
    (
        "phospy.datasets.construction.interpreter",
        "phospy.preprocessing.dataset_construction",
    ): frozenset({"CorrectedPreprocessingOutput"}),
    (
        "phospy.datasets.construction.planning",
        "phospy.preprocessing.contracts",
    ): frozenset({"PreprocessingPlan"}),
    (
        "phospy.datasets.construction.planning",
        "phospy.preprocessing.contracts.policies",
    ): frozenset({"SiteMatrixPolicy"}),
    (
        "phospy.datasets.construction.planning",
        "phospy.preprocessing.planning.boundary",
    ): frozenset({"PreprocessingPlanInterpreter"}),
    (
        "phospy.datasets.construction.planning",
        "phospy.preprocessing.workflow",
    ): frozenset({"GroupCoverageFilterMetadataValidator"}),
    (
        "phospy.datasets.construction.preprocessing",
        "phospy.preprocessing.contracts",
    ): frozenset({"PreprocessingPlan", "PreprocessingStageExecution"}),
    (
        "phospy.datasets.construction.preprocessing",
        "phospy.preprocessing.dataset_construction",
    ): frozenset(
        {
            "BatchCorrectionAdequacyValidatorProtocol",
            "BatchDesignMetadataValidatorProtocol",
            "CorrectedPreprocessingOutput",
            "CorrectedPreprocessingOutputIntegratorProtocol",
            "DatasetConstructionPreprocessingRunner",
            "DatasetConstructionProteinAwareRunner",
            "PreprocessingPipelineProtocol",
            "PreprocessingProvenanceAdapterProtocol",
            "ProteinAwarePreparationStageProtocol",
            "SpsRuvStyleBatchCorrectionRunner",
        }
    ),
    (
        "phospy.datasets.construction.protein_aware",
        "phospy.preprocessing.dataset_construction",
    ): frozenset({"ProteinAwarePreparationStageOutput"}),
    (
        "phospy.datasets.construction.provenance",
        "phospy.preprocessing.contracts",
    ): frozenset(
        {
            "DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION",
            "DATASET_PREPROCESSING_STAGE_COMPARISONS",
            "DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER",
            "DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM",
            "DATASET_PREPROCESSING_STAGE_LOCALISATION",
            "DATASET_PREPROCESSING_STAGE_MISSING_DATA",
            "DATASET_PREPROCESSING_STAGE_NORMALISATION",
            "DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT",
            "DATASET_PREPROCESSING_STAGE_SITE_MATRIX",
            "DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION",
            "DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION",
            "PreprocessingPlan",
            "PreprocessingStageExecution",
            "PreprocessingStageOrderResolution",
            "TotalProteinCorrectionIdentityPolicy",
        }
    ),
    (
        "phospy.datasets.construction.provenance",
        "phospy.preprocessing.workflow",
    ): frozenset(
        {
            "DeterminismKind",
            "PreprocessingStageOrderPolicy",
            "PreprocessingStageProvenance",
            "build_duplicate_site_resolution_policy",
        }
    ),
    (
        "phospy.datasets.construction.report_assembly",
        "phospy.preprocessing.contracts",
    ): frozenset(
        {
            "DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM",
            "DATASET_PREPROCESSING_STAGE_PROTEIN_AWARE_PREPARATION",
            "DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION",
            "PreprocessingPlan",
            "PreprocessingStageExecution",
        }
    ),
    (
        "phospy.datasets.construction.report_assembly",
        "phospy.preprocessing.dataset_construction",
    ): frozenset({"ResolvedBatchCorrectionMetadata", "levels_in_sample_order"}),
    (
        "phospy.datasets.construction.site_sequence_validation",
        "phospy.preprocessing.contracts",
    ): frozenset(
        {
            "DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION",
            "PreprocessingStageExecution",
        }
    ),
    (
        "phospy.datasets.construction.state_assembly",
        "phospy.preprocessing.contracts",
    ): frozenset(
        {
            "DATASET_PREPROCESSING_STAGE_MISSING_DATA",
            "DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION",
            "DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION",
            "PreprocessingPlan",
            "PreprocessingStageExecution",
        }
    ),
    (
        "phospy.datasets.construction.state_assembly",
        "phospy.preprocessing.contracts.policies",
    ): frozenset({"ComparisonBuildingPolicy", "SiteMatrixPolicy"}),
    (
        "phospy.datasets.construction.state_assembly",
        "phospy.preprocessing.dataset_construction",
    ): frozenset(
        {
            "ProcessingTraceDiagnostics",
            "QuantitativeOperationEvidenceContext",
            "dataset_construction_stage_metadata",
            "parse_processing_trace_diagnostics",
            "validate_dataset_quantitative_operation_evidence",
        }
    ),
    (
        "phospy.datasets.construction.transformation_state",
        "phospy.preprocessing.contracts",
    ): frozenset(
        {
            "DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM",
            "PreprocessingPlan",
            "PreprocessingStageExecution",
        }
    ),
    (
        "phospy.datasets.construction.transformation_state",
        "phospy.preprocessing.contracts.policies",
    ): frozenset({"IntensityTransformPolicy"}),
}

DATASET_PREPROCESSING_ALLOWED_TARGETS: Mapping[str, frozenset[str]] = {
    source: frozenset(
        target
        for candidate_source, target in DATASET_PREPROCESSING_ALLOWED_IMPORTS
        if candidate_source == source
    )
    for source in {source for source, _target in DATASET_PREPROCESSING_ALLOWED_IMPORTS}
}


@dataclass(frozen=True, slots=True)
class ModuleOwner:
    source_path: str
    module: str
    owner: str
    rationale: str


@dataclass(frozen=True, slots=True)
class ClassifiedImport:
    record: ImportRecord
    source_owner: str
    target_owner: str
    allowed: bool


def load_module_owners(path: Path) -> tuple[ModuleOwner, ...]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle)
        return tuple(
            ModuleOwner(
                source_path=row["source_path"],
                module=row["module"],
                owner=row["owner"],
                rationale=row["rationale"],
            )
            for row in rows
        )


def mixed_module_declaration_problems(
    *,
    actual_declarations: Iterable[str],
    review: MixedModuleReview,
) -> list[str]:
    """Return completeness and ownership problems for one reviewed mixed file."""

    actual = frozenset(actual_declarations)
    problems = [
        f"unregistered mixed-module declaration: {declaration}"
        for declaration in sorted(actual - review.declarations)
    ]
    problems.extend(
        f"stale mixed-module declaration: {declaration}"
        for declaration in sorted(review.declarations - actual)
    )
    problems.extend(
        f"foreign owner recorded for undeclared definition: {declaration}"
        for declaration in sorted(set(review.foreign_owners) - review.declarations)
    )
    if review.primary_owner not in OWNERS:
        problems.append(f"unknown mixed-module primary owner: {review.primary_owner}")
    for declaration, owner in sorted(review.foreign_owners.items()):
        if owner not in OWNERS:
            problems.append(
                f"unknown mixed-module foreign owner {owner!r}: {declaration}"
            )
        elif owner == review.primary_owner:
            problems.append(
                f"foreign declaration repeats primary owner {owner!r}: {declaration}"
            )
    return problems


def module_inventory_problems(
    *,
    live_paths: Mapping[str, str],
    entries: Iterable[ModuleOwner],
) -> list[str]:
    """Return exact-coverage, uniqueness, and schema problems for an owner map."""

    rows = tuple(entries)
    problems: list[str] = []
    path_counts = Counter(row.source_path for row in rows)
    module_counts = Counter(row.module for row in rows)
    problems.extend(
        f"multiply mapped source path: {path} ({count} entries)"
        for path, count in sorted(path_counts.items())
        if count != 1
    )
    problems.extend(
        f"multiply mapped module: {module} ({count} entries)"
        for module, count in sorted(module_counts.items())
        if count != 1
    )
    mapped_paths = set(path_counts)
    live_path_set = set(live_paths)
    problems.extend(
        f"unmapped source module: {path}"
        for path in sorted(live_path_set - mapped_paths)
    )
    problems.extend(
        f"stale mapped source module: {path}"
        for path in sorted(mapped_paths - live_path_set)
    )
    for row in rows:
        expected_module = live_paths.get(row.source_path)
        if expected_module is not None and row.module != expected_module:
            problems.append(
                f"module/path mismatch: {row.source_path} maps {row.module}; "
                f"live module is {expected_module}"
            )
        if row.owner not in OWNERS:
            problems.append(f"unknown owner {row.owner!r}: {row.source_path}")
        if not row.rationale.strip():
            problems.append(f"missing ownership rationale: {row.source_path}")
    return problems


def classify_imports(
    graph: ImportGraph,
    *,
    project_root: Path,
    owners_by_module: Mapping[str, str],
) -> tuple[ClassifiedImport, ...]:
    classified: list[ClassifiedImport] = []
    for record in graph.records:
        source_owner = owners_by_module[record.source_module]
        target_owner = owners_by_module[record.target]
        allowed = is_allowed_import(
            source_owner=source_owner,
            target_owner=target_owner,
            source_module=record.source_module,
            target_module=record.target,
        )
        classified.append(
            ClassifiedImport(
                record=record,
                source_owner=source_owner,
                target_owner=target_owner,
                allowed=allowed,
            )
        )
    return tuple(classified)


def logical_owner_edges(
    module_edges: Mapping[str, Iterable[str]],
    *,
    owners_by_module: Mapping[str, str],
) -> frozenset[tuple[str, str]]:
    """Collapse module imports to cross-owner dependency edges."""

    return frozenset(
        (owners_by_module[source], owners_by_module[target])
        for source, targets in module_edges.items()
        for target in targets
        if owners_by_module[source] != owners_by_module[target]
    )


def logical_owner_graph_problems(
    owner_edges: Iterable[tuple[str, str]],
) -> list[str]:
    """Return schema and cycle problems in a proposed logical-owner graph."""

    edges = frozenset(owner_edges)
    problems: list[str] = []
    for source, target in sorted(edges):
        if source not in OWNERS:
            problems.append(f"unknown source owner in permitted edge: {source!r}")
        if target not in OWNERS:
            problems.append(f"unknown target owner in permitted edge: {target!r}")

    adjacency: dict[str, set[str]] = {}
    for source, target in edges:
        if source not in OWNERS or target not in OWNERS or source == target:
            continue
        adjacency.setdefault(source, set()).add(target)
        adjacency.setdefault(target, set())

    index = 0
    stack: list[str] = []
    on_stack: set[str] = set()
    indices: dict[str, int] = {}
    lowlinks: dict[str, int] = {}
    components: list[tuple[str, ...]] = []

    def strong_connect(owner: str) -> None:
        nonlocal index
        indices[owner] = index
        lowlinks[owner] = index
        index += 1
        stack.append(owner)
        on_stack.add(owner)

        for target in sorted(adjacency.get(owner, set())):
            if target not in indices:
                strong_connect(target)
                lowlinks[owner] = min(lowlinks[owner], lowlinks[target])
            elif target in on_stack:
                lowlinks[owner] = min(lowlinks[owner], indices[target])

        if lowlinks[owner] != indices[owner]:
            return
        component: list[str] = []
        while True:
            target = stack.pop()
            on_stack.remove(target)
            component.append(target)
            if target == owner:
                break
        components.append(tuple(sorted(component)))

    for owner in sorted(adjacency):
        if owner not in indices:
            strong_connect(owner)

    problems.extend(
        f"logical-owner dependency cycle: {', '.join(component)}"
        for component in sorted(components)
        if len(component) > 1
    )
    return problems


def is_allowed_import(
    *,
    source_owner: str,
    target_owner: str,
    source_module: str,
    target_module: str,
) -> bool:
    """Apply owner direction and contract-path rules to one resolved import."""

    if source_owner == target_owner:
        return True
    if (source_module, target_module) == DATASET_PATH_READER_EDGE:
        return True
    if (source_owner, target_owner) not in PERMITTED_LOGICAL_OWNER_EDGES:
        return False
    if source_owner == "public_facades":
        return _public_facade_target(target_owner, target_module)
    if source_owner in SHARED_ALLOWED_TARGETS:
        return True
    if target_owner in PASSIVE_RECORD_OWNERS:
        return _preparation_record_target(source_owner, source_module, target_module)
    if source_owner in CAPABILITY_OWNERS and target_owner in (
        SHARED_OWNERS | FOUNDATION_OWNERS
    ):
        return True
    if source_owner == "kinase" and target_owner == "activities":
        return target_module in {
            "phospy.activities.contracts",
            "phospy.activities.semantics",
            "phospy.activities.validation",
        }
    if source_owner == "datasets" and target_owner == "preprocessing":
        return _dataset_to_preprocessing_contract(source_module, target_module)
    if source_owner == "run_lineage":
        return _run_lineage_target(target_owner, target_module)
    if source_owner == "datasets" and target_owner == "run_lineage":
        return source_module != "phospy.datasets" and target_module in {
            "phospy.run_lineage.models",
            "phospy.run_lineage.serialization",
        }
    if target_owner == "run_lineage":
        return target_module in {
            "phospy.run_lineage",
            "phospy.run_lineage.models",
            "phospy.run_lineage.serialization",
        }
    if source_owner == "differential" and target_owner == "datasets":
        return _dataset_contract_target(target_module)
    if source_owner == "differential" and target_owner == "preprocessing":
        return target_module in {
            "phospy.preprocessing.contracts.configuration.total_protein",
            "phospy.preprocessing.contracts.policies",
            "phospy.preprocessing.protein_aware.alignment",
        }
    if source_owner == "kinase" and target_owner == "datasets":
        return _dataset_contract_target(target_module)
    if source_owner == "signalome" and target_owner == "kinase":
        return _kinase_result_contract_target(target_module)
    if source_owner == "signalome" and target_owner == "datasets":
        return _signalome_dataset_contract_target(source_module, target_module)
    if source_owner == "io" and target_owner in FOUNDATION_OWNERS:
        return True
    if source_owner == "io":
        return _io_adapter_target(target_owner, target_module)
    return False


def _preparation_record_target(
    source_owner: str,
    source: str,
    target: str,
) -> bool:
    if source_owner in {"datasets", "preprocessing"}:
        return target in PREPARATION_RECORD_ALLOWED_SYMBOLS
    if source_owner == "differential":
        return target in {
            "phospy.differential.planning.protein_aware_inputs": {
                "phospy.preparation_records.policies",
                "phospy.preparation_records.protein_aware",
            },
            "phospy.differential.result_assembly": {
                "phospy.preparation_records.batch_correction",
            },
        }.get(source, set())
    if source_owner == "io":
        return source == "phospy.io.bundles._shared.processing_state" and target == (
            "phospy.preparation_records.policies"
        )
    return False


def _dataset_to_preprocessing_contract(source: str, target: str) -> bool:
    return target in DATASET_PREPROCESSING_ALLOWED_TARGETS.get(source, frozenset())


def _signalome_dataset_contract_target(source: str, target: str) -> bool:
    return target in {
        "phospy.signalome.results": {
            "phospy.datasets.models",
            "phospy.datasets.processing_state",
        },
        "phospy.signalome.result_validation": {
            "phospy.datasets.internal_view",
            "phospy.datasets.models",
        },
        "phospy.signalome.contracts": {
            "phospy.datasets.internal_view",
            "phospy.datasets.models",
        },
        "phospy.signalome.validator": {"phospy.datasets.internal_view"},
    }.get(source, set())


def _dataset_contract_target(target: str) -> bool:
    published_modules = {
        "phospy.datasets.contracts",
        "phospy.datasets.internal_view",
        "phospy.datasets.models",
        "phospy.datasets.processing_state",
        "phospy.datasets.derived_quantitative",
        "phospy.datasets.imputation_metadata",
        "phospy.datasets.trusted_construction_caveats",
    }
    return target in published_modules


def _preprocessing_evidence_contract_target(target: str) -> bool:
    published_modules = {
        "phospy.preprocessing.contracts",
        "phospy.preprocessing.contracts.configuration.dataset",
        "phospy.preprocessing.contracts.configuration.total_protein",
        "phospy.preprocessing.contracts.policies",
        "phospy.preprocessing.evidence.missing_data",
        "phospy.preprocessing.planning.boundary",
        "phospy.preprocessing.protein_aware.alignment",
        "phospy.preprocessing.provenance.composition",
        "phospy.preprocessing.provenance.models",
        "phospy.preprocessing.provenance.policies",
        "phospy.preprocessing.results.diagnostics",
        "phospy.preprocessing.results.protein_aware",
        "phospy.preprocessing.results.quantitative_evidence",
        "phospy.preparation_records.report_schema",
        "phospy.preprocessing.validation.configuration",
        "phospy.preprocessing.validation.group_coverage",
        "phospy.preparation_records.validation.sample_groups",
        "phospy.preprocessing.workflow",
    }
    return target in published_modules


def _run_lineage_target(target_owner: str, target: str) -> bool:
    allowed = {
        "evidence": {"phospy.evidence.serialization.historical_peptide"},
        "preprocessing": {
            "phospy.preprocessing.provenance.models",
            "phospy.preprocessing.provenance.serialization",
        },
        "provenance": {
            "phospy.datasets.derived_quantitative",
            "phospy.provenance.environment",
            "phospy.provenance.immutability",
            "phospy.provenance.models",
            "phospy.provenance.models._shared",
            "phospy.provenance.models.environment",
            "phospy.provenance.models.tables",
            "phospy.provenance.scientific_policy_models",
            "phospy.provenance.serialization",
            "phospy.provenance.serialization._payload",
            "phospy.provenance.serialization.environment",
            "phospy.provenance.serialization.tables",
        },
        "references": {
            "phospy.references.provenance.context",
            "phospy.references.provenance.models",
            "phospy.references.provenance.serialization",
        },
    }
    if target_owner in {"errors", "frames", "taxonomy"}:
        return True
    return target in allowed.get(target_owner, set())


def _kinase_result_contract_target(target: str) -> bool:
    published_modules = {
        "phospy.kinase.contracts",
        "phospy.kinase.requests",
        "phospy.kinase.result_boundary",
        "phospy.kinase.results",
        "phospy.kinase.scoring.results",
        "phospy.kinase.prediction.tables",
    }
    return target in published_modules


def _io_adapter_target(target_owner: str, target: str) -> bool:
    return target in IO_ADAPTER_ALLOWED_TARGETS.get(target_owner, frozenset())


def _public_facade_target(target_owner: str, target: str) -> bool:
    return target in PUBLIC_FACADE_ALLOWED_TARGETS.get(target_owner, frozenset())

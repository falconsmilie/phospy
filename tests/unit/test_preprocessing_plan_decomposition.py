from __future__ import annotations

import ast
from dataclasses import fields
from pathlib import Path

import pandas as pd
import pytest

from phospy.errors.input import PhosPyInputError
from phospy.science.configs.dataset import DatasetPreprocessingConfig
from phospy.science.configs.preprocessing import (
    CorrectionMissingnessPolicy,
    DatasetTotalProteinCorrectionIdentityConfig,
    InternalBatchCorrectionControlSiteMode,
    InternalBatchCorrectionControlSiteSource,
    InternalBatchCorrectionImputationPolicy,
    InternalBatchCorrectionMethod,
    InternalBatchCorrectionMissingValuePolicy,
    InternalBatchCorrectionRequest,
    InternalBatchCorrectionStageOrder,
)
from phospy.science.datasets.preprocessing.imputation_scale_policy import (
    resolve_imputation_scale_policy,
)
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.plan_assembly import (
    PreprocessingPlanAssembler,
)
from phospy.science.datasets.preprocessing.plan_constants import (
    PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
    PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA,
)
from phospy.science.datasets.preprocessing.plan_interpreter import (
    PreprocessingPlanInterpreter,
)
from phospy.science.datasets.preprocessing.plan_resolved import (
    ResolvedBatchCorrectionPlanFields,
    ResolvedCoreTransformPlanFields,
    ResolvedGroupCoveragePlanFields,
    ResolvedImputationScalePlanFields,
    ResolvedLocalisationPlanFields,
    ResolvedPreprocessingPlanFields,
    ResolvedRuvReadinessPlanFields,
    ResolvedSiteMatrixComparisonPlanFields,
    ResolvedSiteSequencePlanFields,
    ResolvedStageOrderPlanFields,
    ResolvedTotalProteinCorrectionPlanFields,
)
from phospy.science.datasets.preprocessing.plan_rules import (
    PreprocessingBatchCorrectionPlanRuleFamily,
    PreprocessingCorePlanPolicyRuleFamily,
    PreprocessingRuvReadinessPlanRuleFamily,
    PreprocessingSiteMatrixComparisonPlanRuleFamily,
    PreprocessingSiteSequencePlanRuleFamily,
)
from phospy.science.datasets.preprocessing.plan_stage_order import (
    PreprocessingStageOrderResolution,
)
from phospy.science.datasets.preprocessing.policy_models import (
    ComparisonBuildingPolicy,
    ImputationInputScale,
    IntensityTransformPolicy,
    LocalisationEligibilityMode,
    MissingDataPolicy,
    NormalisationPolicy,
    SiteMatrixDuplicateSitePolicy,
    SiteMatrixMissingDataPolicy,
    SiteMatrixPolicy,
    SiteSequenceConflictPolicy,
    SiteSequenceResolutionMode,
    TotalProteinCorrectionPolicy,
)
from phospy.science.datasets.preprocessing.stages.comparisons import ComparisonsStage
from phospy.science.datasets.preprocessing.total_protein_identity import (
    DEFAULT_TOTAL_PROTEIN_CORRECTION_IDENTITY_POLICY,
    TotalProteinCorrectionIdentityResolver,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def test_core_rule_family_resolves_transform_normalisation_and_missing_policy() -> None:
    resolved = PreprocessingCorePlanPolicyRuleFamily().run(
        intensity_transform_policy="log2",  # type: ignore[arg-type]
        normalisation_policy="median_center",  # type: ignore[arg-type]
        missing_data_policy="impute_row_median",  # type: ignore[arg-type]
    )

    assert resolved.intensity_transform_policy is IntensityTransformPolicy.LOG2
    assert resolved.normalisation_policy is NormalisationPolicy.MEDIAN_CENTER
    assert resolved.missing_data_policy is MissingDataPolicy.IMPUTE_ROW_MEDIAN


def test_imputation_scale_resolution_is_policy_specific() -> None:
    minprob = resolve_imputation_scale_policy(
        missing_data_policy=MissingDataPolicy.IMPUTE_MINPROB,
        requested_input_scale=None,
        intensity_transform_policy=IntensityTransformPolicy.LOG2,
    )
    row_median = resolve_imputation_scale_policy(
        missing_data_policy=MissingDataPolicy.IMPUTE_ROW_MEDIAN,
        requested_input_scale="linear",
        intensity_transform_policy=IntensityTransformPolicy.LOG2,
    )
    knn = resolve_imputation_scale_policy(
        missing_data_policy=MissingDataPolicy.IMPUTE_KNN,
        requested_input_scale="log2",
        intensity_transform_policy=IntensityTransformPolicy.LOG2,
    )

    assert minprob.input_scale is ImputationInputScale.LOG2
    assert minprob.input_scale_source == "method_required"
    assert minprob.operation_order == "after_intensity_transform"
    assert row_median.input_scale is ImputationInputScale.LINEAR
    assert row_median.operation_order == "before_intensity_transform"
    assert knn.input_scale is ImputationInputScale.LOG2
    assert knn.operation_order == "after_intensity_transform"


def test_site_sequence_rule_resolves_mode_conflict_and_columns() -> None:
    resolved = PreprocessingSiteSequencePlanRuleFamily().run(
        site_sequence_resolution_enabled=True,
        site_sequence_resolution_fasta_path="reference.fasta",
        site_sequence_resolution_mode="replace_existing",  # type: ignore[arg-type]
        site_sequence_resolution_conflict_policy="replace_existing",  # type: ignore[arg-type]
        site_sequence_resolution_flank_size=7,
        site_sequence_resolution_accession_column=" protein_accession ",
        site_sequence_resolution_site_column=" site ",
    )

    assert resolved.site_sequence_resolution_enabled is True
    assert resolved.site_sequence_resolution_mode is (
        SiteSequenceResolutionMode.REPLACE_EXISTING
    )
    assert resolved.site_sequence_resolution_conflict_policy is (
        SiteSequenceConflictPolicy.REPLACE_EXISTING
    )
    assert resolved.site_sequence_resolution_accession_column == "protein_accession"
    assert resolved.site_sequence_resolution_site_column == "site"


def test_site_matrix_and_comparison_rule_resolves_known_fields() -> None:
    resolved = PreprocessingSiteMatrixComparisonPlanRuleFamily().run(
        site_matrix_policy="build_from_metadata",  # type: ignore[arg-type]
        site_matrix_duplicate_site_policy="aggregate_mean",  # type: ignore[arg-type]
        site_matrix_missing_data_policy="drop_any_missing",  # type: ignore[arg-type]
        site_matrix_minimum_observed_values=None,
        comparison_building_policy="sample_metadata_pairs",  # type: ignore[arg-type]
        comparison_sample_group_column=" comparison_group ",
        comparison_pairs=(("control", "treated"),),
    )

    assert resolved.site_matrix_policy is SiteMatrixPolicy.BUILD_FROM_METADATA
    assert resolved.site_matrix_duplicate_site_policy is (
        SiteMatrixDuplicateSitePolicy.AGGREGATE_MEAN
    )
    assert resolved.site_matrix_missing_data_policy is (
        SiteMatrixMissingDataPolicy.DROP_ANY_MISSING
    )
    assert resolved.comparison_building_policy is (
        ComparisonBuildingPolicy.SAMPLE_METADATA_PAIRS
    )
    assert resolved.comparison_sample_group_column == "comparison_group"
    assert resolved.comparison_pairs == (("control", "treated"),)


def test_total_protein_direct_identity_resolver_normalises_keys() -> None:
    resolved = TotalProteinCorrectionIdentityResolver().run(
        DatasetTotalProteinCorrectionIdentityConfig(
            phosphosite_key=" gene_symbol ",
            total_protein_key=" __index__ ",
        )
    )

    assert resolved.mode == "direct"
    assert resolved.phosphosite_key == "gene_symbol"
    assert resolved.total_protein_key == "__index__"
    assert resolved.mapping_table is None
    assert resolved.mapping_table_fingerprint is None


def test_total_protein_mapping_table_identity_resolver_normalises_and_fingerprints() -> (
    None
):
    mapping_table = pd.DataFrame(
        {
            "phosphosite_id": [" MAPK14 ", pd.NA],
            "protein_id": [" MAPK14_TOTAL ", " AKT1_TOTAL "],
        }
    )

    resolved = TotalProteinCorrectionIdentityResolver().run(
        DatasetTotalProteinCorrectionIdentityConfig(
            mode="mapping_table",
            mapping_table=mapping_table,
            mapping_phosphosite_key="phosphosite_id",
            mapping_total_protein_key="protein_id",
        )
    )

    assert resolved.mapping_table == (
        ("MAPK14", "MAPK14_TOTAL"),
        ("", "AKT1_TOTAL"),
    )
    assert resolved.mapping_phosphosite_key == "phosphosite_id"
    assert resolved.mapping_total_protein_key == "protein_id"
    assert isinstance(resolved.mapping_table_fingerprint, str)
    assert len(resolved.mapping_table_fingerprint) == 64


def test_total_protein_identity_resolver_rejects_tampered_unknown_mode() -> None:
    config = DatasetTotalProteinCorrectionIdentityConfig()
    object.__setattr__(config, "mode", "future_mode")

    with pytest.raises(PhosPyInputError, match="unsupported mode"):
        TotalProteinCorrectionIdentityResolver().run(config)


def test_total_protein_identity_resolver_rejects_missing_mapping_column() -> None:
    config = DatasetTotalProteinCorrectionIdentityConfig(
        mode="mapping_table",
        mapping_table=pd.DataFrame({"phosphosite_id": ["MAPK14"]}),
        mapping_phosphosite_key="phosphosite_id",
        mapping_total_protein_key="protein_id",
    )

    with pytest.raises(PhosPyInputError, match="missing column 'protein_id'"):
        TotalProteinCorrectionIdentityResolver().run(config)


def test_ruv_readiness_rule_resolves_report_only_fields() -> None:
    resolved = PreprocessingRuvReadinessPlanRuleFamily().run(
        ruv_readiness_enabled=True,
        ruv_readiness_control_feature_column=" is_control ",
        ruv_readiness_replicate_group_column=" replicate_group ",
        ruv_readiness_batch_column=" batch ",
    )

    assert resolved.ruv_readiness_enabled is True
    assert resolved.ruv_readiness_control_feature_column == "is_control"
    assert resolved.ruv_readiness_replicate_group_column == "replicate_group"
    assert resolved.ruv_readiness_batch_column == "batch"


def test_fixed_effect_batch_rule_resolves_metadata_fields() -> None:
    resolved = PreprocessingBatchCorrectionPlanRuleFamily().run(
        batch_correction_method="linear_residualize_batch",
        batch_correction_batch_column=" batch ",
        batch_correction_condition_column=" condition ",
        batch_correction_condition_columns=(" condition ",),
        batch_correction_replicate_column=None,
        batch_correction_control_site_set=None,
        batch_correction_missingness_policy=None,
        batch_correction_internal_request=None,
        batch_correction_preserve_condition_effects=True,
    )

    assert resolved.batch_correction_method == "linear_residualize_batch"
    assert resolved.batch_correction_batch_column == "batch"
    assert resolved.batch_correction_condition_column == "condition"
    assert resolved.batch_correction_condition_columns == ("condition",)
    assert resolved.batch_correction_preserve_condition_effects is True


def test_sps_ruv_batch_rule_resolves_control_missingness_and_internal_request() -> None:
    request = _internal_sps_request()
    missingness_policy = CorrectionMissingnessPolicy()
    control_site_set = object()

    resolved = PreprocessingBatchCorrectionPlanRuleFamily().run(
        batch_correction_method="sps_ruv_style",
        batch_correction_batch_column="batch",
        batch_correction_condition_column="condition",
        batch_correction_condition_columns=("condition",),
        batch_correction_replicate_column="replicate",
        batch_correction_control_site_set=control_site_set,
        batch_correction_missingness_policy=missingness_policy,
        batch_correction_internal_request=request,
        batch_correction_preserve_condition_effects=True,
    )

    assert resolved.batch_correction_method == "sps_ruv_style"
    assert resolved.batch_correction_replicate_column == "replicate"
    assert resolved.batch_correction_control_site_set is control_site_set
    assert resolved.batch_correction_missingness_policy is missingness_policy
    assert resolved.batch_correction_internal_request is request


def test_sps_ruv_batch_rule_reports_deterministic_first_missing_requirement() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="requires batch_correction_control_site_set",
    ):
        PreprocessingBatchCorrectionPlanRuleFamily().run(
            batch_correction_method="sps_ruv_style",
            batch_correction_batch_column="batch",
            batch_correction_condition_column="condition",
            batch_correction_condition_columns=("condition",),
            batch_correction_replicate_column=None,
            batch_correction_control_site_set=None,
            batch_correction_missingness_policy=None,
            batch_correction_internal_request=None,
            batch_correction_preserve_condition_effects=True,
        )


def test_direct_plan_rejects_inconsistent_resolved_imputation_state() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing_data_input_scale_source .* inconsistent",
    ):
        PreprocessingPlan(
            missing_data_policy=MissingDataPolicy.IMPUTE_ROW_MEDIAN,
            missing_data_input_scale=ImputationInputScale.LINEAR,
            missing_data_input_scale_source="method_required",
            missing_data_min_observed_values=1,
            stage_order=("missing_data",),
        )


def test_direct_plan_rejects_incomplete_sps_ruv_batch_state() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="requires batch_correction_control_site_set",
    ):
        PreprocessingPlan(
            batch_correction_method="sps_ruv_style",
            stage_order=("batch_correction",),
        )


def test_direct_plan_rejects_invalid_total_protein_identity_state() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="total_protein_correction_identity_policy",
    ):
        PreprocessingPlan(
            total_protein_correction_policy=(
                TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
            ),
            total_protein_correction_identity_policy=object(),  # type: ignore[arg-type]
            stage_order=("total_protein_correction",),
        )


def test_comparison_stage_rejects_policy_without_sample_grouping_metadata() -> None:
    state = PreprocessingState(
        phospho=pd.DataFrame({"sample_a": [1.0], "sample_b": [2.0]}),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["SEQ"],
            }
        ),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            comparison_building_policy=(ComparisonBuildingPolicy.SAMPLE_METADATA_PAIRS),
            stage_order=("comparisons",),
        ),
    )

    with pytest.raises(PhosPyInputError, match="requires sample_metadata"):
        ComparisonsStage().run(state)


def test_explicit_assembler_matches_interpreter_for_default_plan() -> None:
    assembled = PreprocessingPlanAssembler().run(_default_resolved_fields())
    interpreted = PreprocessingPlanInterpreter().run(DatasetPreprocessingConfig())

    assert assembled == interpreted


def test_resolved_preprocessing_plan_fields_are_typed_sections() -> None:
    section_names = tuple(
        field.name for field in fields(ResolvedPreprocessingPlanFields)
    )

    assert section_names == (
        "core",
        "imputation",
        "localisation",
        "site_sequence",
        "group_coverage",
        "total_protein",
        "site_matrix_comparisons",
        "ruv_readiness",
        "batch_correction",
        "stage_order",
    )


def test_plan_model_no_longer_contains_dynamic_resolution_path() -> None:
    source = _source("src/phospy/science/datasets/preprocessing/plan.py")

    assert "_set_resolved_plan_fields" not in source
    assert "_resolve_total_correction_identity_policy" not in source
    assert "class PreprocessingPlanInterpreter" not in source
    assert "PreprocessingConfigPolicyResolver" not in source
    assert "PreprocessingBatchCorrectionPlanRuleFamily" not in source
    assert "hash_table_tolerance" not in source
    assert "pandas" not in source


def test_plan_assembler_does_not_use_dynamic_field_copying() -> None:
    source = _source("src/phospy/science/datasets/preprocessing/plan_assembly.py")

    assert "fields(" not in source
    assert "setattr" not in source
    assert "dict[" not in source
    assert "to_dict" not in source


def test_interpreter_run_coordinates_typed_resolution_and_assembly() -> None:
    source = _source("src/phospy/science/datasets/preprocessing/plan_interpreter.py")
    tree = ast.parse(source)
    run_node = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name == "run"
    )
    run_source = ast.get_source_segment(source, run_node)
    assert run_source is not None

    assert "ResolvedPreprocessingPlanFields" in run_source
    assert "PreprocessingPlanAssembler" in run_source
    assert "TotalProteinCorrectionIdentityResolver" in run_source
    assert "_resolve_total_correction_identity_policy" not in run_source
    assert "_set_resolved_plan_fields" not in run_source


def test_plan_config_resolution_no_longer_owns_total_identity_mode_errors() -> None:
    source = _source(
        "src/phospy/science/datasets/preprocessing/plan_config_resolution.py"
    )

    assert "reject_unknown_total_correction_identity_mode" not in source
    assert "identity contains an unsupported mode" not in source


def _default_resolved_fields() -> ResolvedPreprocessingPlanFields:
    return ResolvedPreprocessingPlanFields(
        core=ResolvedCoreTransformPlanFields(
            intensity_transform_policy=IntensityTransformPolicy.IDENTITY,
            intensity_transform_pseudocount=1.0,
            normalisation_policy=NormalisationPolicy.NONE,
            missing_data_policy=MissingDataPolicy.FORBID,
            missing_data_min_observed_values=None,
            missing_data_q=None,
            missing_data_width=None,
            missing_data_seed=None,
            missing_data_k=None,
            missing_data_distance=None,
            missing_data_max_missing_fraction_per_row=None,
        ),
        imputation=ResolvedImputationScalePlanFields(
            missing_data_input_scale=None,
            missing_data_input_scale_source=None,
            missing_data_imputation_operation_order=None,
        ),
        localisation=ResolvedLocalisationPlanFields(
            localisation_mode=LocalisationEligibilityMode.REQUIRE_THRESHOLD,
            localisation_min_confidence=0.75,
            localisation_confidence_column="localisation_confidence",
            localisation_waiver_reason=None,
        ),
        site_sequence=ResolvedSiteSequencePlanFields(
            site_sequence_resolution_enabled=False,
            site_sequence_resolution_fasta_path=None,
            site_sequence_resolution_mode=(
                SiteSequenceResolutionMode.VALIDATE_EXISTING_AND_FILL_MISSING
            ),
            site_sequence_resolution_conflict_policy=(
                SiteSequenceConflictPolicy.PRESERVE_EXISTING
            ),
            site_sequence_resolution_flank_size=7,
            site_sequence_resolution_accession_column="protein_accession",
            site_sequence_resolution_site_column="site",
        ),
        group_coverage=ResolvedGroupCoveragePlanFields(
            group_coverage_filter_enabled=False,
            group_coverage_filter_group_column=None,
            group_coverage_filter_min_finite_observations_per_group=None,
            group_coverage_filter_min_finite_fraction_per_group=None,
            group_coverage_filter_min_groups_passing_threshold=1,
        ),
        total_protein=ResolvedTotalProteinCorrectionPlanFields(
            total_protein_correction_policy=TotalProteinCorrectionPolicy.NONE,
            total_protein_correction_identity_policy=(
                DEFAULT_TOTAL_PROTEIN_CORRECTION_IDENTITY_POLICY
            ),
            protein_aware_preparation_policy="disabled",
            protein_aware_preparation_mapping_policy="require_unambiguous",
        ),
        site_matrix_comparisons=ResolvedSiteMatrixComparisonPlanFields(
            site_matrix_policy=SiteMatrixPolicy.AS_INPUT,
            site_matrix_duplicate_site_policy=SiteMatrixDuplicateSitePolicy.ERROR,
            site_matrix_missing_data_policy=(
                SiteMatrixMissingDataPolicy.DROP_ANY_MISSING
            ),
            site_matrix_minimum_observed_values=None,
            comparison_building_policy=ComparisonBuildingPolicy.NONE,
            comparison_sample_group_column="comparison_group",
            comparison_pairs=None,
        ),
        ruv_readiness=ResolvedRuvReadinessPlanFields(
            ruv_readiness_enabled=False,
            ruv_readiness_control_feature_column="is_control_feature",
            ruv_readiness_replicate_group_column="replicate_group",
            ruv_readiness_batch_column="batch",
        ),
        batch_correction=ResolvedBatchCorrectionPlanFields(
            batch_correction_method="none",
            batch_correction_batch_column="batch",
            batch_correction_condition_column="condition",
            batch_correction_condition_columns=("condition",),
            batch_correction_replicate_column=None,
            batch_correction_control_site_set=None,
            batch_correction_missingness_policy=None,
            batch_correction_internal_request=None,
            batch_correction_preserve_condition_effects=True,
        ),
        stage_order=ResolvedStageOrderPlanFields(
            stage_order=("localisation_confidence", "missing_data"),
            stage_order_resolution=(
                PreprocessingStageOrderResolution(
                    stage="localisation_confidence",
                    order_index=0,
                    rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
                ),
                PreprocessingStageOrderResolution(
                    stage="missing_data",
                    order_index=1,
                    rationale=(
                        PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA
                    ),
                ),
            ),
        ),
    )


def _internal_sps_request() -> InternalBatchCorrectionRequest:
    return InternalBatchCorrectionRequest(
        method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
        batch_column="batch",
        condition_columns=("condition",),
        replicate_column="replicate",
        control_site_source=InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED,
        control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
        missing_value_policy=InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING,
        imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
        n_unwanted_factors=1,
        stage_order=InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM,
        diagnostics_enabled=True,
    )


def _source(relative_path: str) -> str:
    return (PROJECT_ROOT / relative_path).read_text(encoding="utf-8")

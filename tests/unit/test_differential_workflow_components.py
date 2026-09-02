from __future__ import annotations

import inspect
from collections.abc import Mapping
from dataclasses import fields, replace
from typing import Any, cast

import numpy as np
import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.advanced import (
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    DifferentialAnalysisConfig,
    DifferentialProteinAwareModelConfig,
    MultipleTestingConfig,
    TechnicalReplicatePolicy,
)
from phospy.api import (
    ContinuousCovariate,
    Contrast,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    FixedEffectCovariate,
    Organism,
    SampleDesignRecord,
)
from phospy.api.results import DifferentialAnalysisResult
from phospy.contracts.configs import (
    DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE,
)
from phospy.errors import (
    PhosPyInputError,
    WorkflowBoundaryError,
    WorkflowValidationError,
)
from phospy.provenance.hashing import fingerprint_optional_table_strict
from phospy.provenance.serialization.tables import table_fingerprint_to_payload
from phospy.science.configs.differential import (
    DifferentialImputedValuePolicy,
    PairedDesignPolicy,
)
from phospy.science.datasets.internal_view import (
    DatasetInternalView,
    ProteinAwarePreparationInternalView,
)
from phospy.science.datasets.models import DatasetPreprocessingReport
from phospy.science.datasets.preprocessing.protein_aware_alignment import (
    PROTEIN_AWARE_REASON_AMBIGUOUS_PROTEIN_MAPPING,
    PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW,
    ProteinAwarePreparationEligibility,
    ProteinAwareSampleAlignmentDiagnostics,
    ProteinAwareTransformationStateDiagnostics,
)
from phospy.science.datasets.preprocessing.protein_aware_preparation import (
    ProteinAwarePreparationReport,
    ProteinAwarePreparationResult,
    ProteinAwareSiteEligibility,
)
from phospy.science.datasets.preprocessing.protein_mapping import ProteinMappingStatus
from phospy.science.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialComputationExecutor,
)
from phospy.science.differential.internal_view import (
    DifferentialComputationResultInternalView,
)
from phospy.science.differential.linear_model import decompose_differential_design
from phospy.science.differential.models import EmpiricalBayesConfig
from phospy.science.differential.models.tables import (
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK,
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
)
from phospy.science.transformations.models import (
    IntensityScaleState,
    MatrixIntensityScaleState,
)
from phospy.validation.workflows.differential import (
    ValidatedExperimentalDesignContract,
)
from phospy.workflows.differential.executor import DifferentialAnalysisExecutor
from phospy.workflows.differential.interpreter import DifferentialAnalysisInterpreter
from phospy.workflows.differential.models import (
    InterpretedDifferentialAnalysisRequest,
    ProteinAwareDifferentialResolvedInputs,
    ResolvedDifferentialExecutionConfig,
    ValidatedDifferentialAnalysisRequest,
)
from phospy.workflows.differential.protein_aware_inputs import (
    ProteinAwareDifferentialInputResolver,
)
from phospy.workflows.differential.provenance import (
    build_differential_policy_provenance,
)
from phospy.workflows.differential.replicates import (
    TechnicalReplicateAggregationGroup,
    TechnicalReplicateAggregationPlan,
)
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.processing_state import (
    imputed_processing_state as valid_imputed_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns
from tests.support.unsafe_dataset_states import (
    unsafe_mark_dataset_total_protein_correction_applied,
    unsafe_remove_dataset_total_matrix,
    unsafe_replace_dataset_intensity_scale_state,
    unsafe_replace_dataset_phospho_scale_state,
    unsafe_replace_dataset_total_scale_state,
)


def _dataset() -> AnalysisReadyPhosphoDataset:
    genes = ["MAPK14", "GSK3B", "AKT1"]
    sites = ["Y182", "S9", "T308"]
    site_index = protein_site_key_index(protein_identifiers=genes, sites=sites)
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 1.0],
            "A_2": [1.1, 2.1, 1.1],
            "B_1": [2.1, 2.0, 1.0],
            "B_2": [2.0, 2.2, 0.9],
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
            **site_key_context_columns(site_index),
            "gene_symbol": genes,
            "site": sites,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15) for site in sites
            ],
            "protein_id": genes,
        },
        index=phospho.index.copy(),
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _request() -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=_dataset(),
        design=ExperimentalDesign(
            samples=(
                SampleDesignRecord(
                    sample_id="A_1",
                    condition="A",
                    biological_replicate_id="A_r1",
                ),
                SampleDesignRecord(
                    sample_id="A_2",
                    condition="A",
                    biological_replicate_id="A_r2",
                ),
                SampleDesignRecord(
                    sample_id="B_1",
                    condition="B",
                    biological_replicate_id="B_r1",
                ),
                SampleDesignRecord(
                    sample_id="B_2",
                    condition="B",
                    biological_replicate_id="B_r2",
                ),
            )
        ),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
    )


def _protein_aware_preparation_for_dataset(
    dataset: AnalysisReadyPhosphoDataset,
) -> ProteinAwarePreparationResult:
    phospho = dataset.phospho
    site_metadata = dataset.site_metadata
    total = dataset.total
    if total is None:
        raise AssertionError("protein-aware preparation test dataset requires total")
    site_keys = tuple(phospho.index.astype(str).tolist())
    total_row_keys = tuple(
        str(site_metadata.loc[site_key, "protein_identifier"]) for site_key in site_keys
    )
    eligibility_status = (
        ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION
    )
    site_eligibility = tuple(
        ProteinAwareSiteEligibility(
            site_key=site_key,
            eligibility=eligibility_status,
            mapping_status=ProteinMappingStatus.MATCHED,
            protein_identifier=str(site_metadata.loc[site_key, "protein_identifier"]),
            total_protein_row_key=total_row_keys[position],
            reasons=("matched_protein_available",),
        )
        for position, site_key in enumerate(site_keys)
    )
    report = ProteinAwarePreparationReport(
        site_eligibility=site_eligibility,
        sample_alignment=ProteinAwareSampleAlignmentDiagnostics(
            phospho_sample_columns=tuple(phospho.columns.astype(str).tolist()),
            total_protein_sample_columns=tuple(phospho.columns.astype(str).tolist()),
            exact_sample_order_match=True,
            sample_order_compatible=True,
            reordered_sample_columns=False,
            allow_reordered_samples=False,
            missing_total_protein_samples=(),
            extra_total_protein_samples=(),
        ),
        transformation_state=ProteinAwareTransformationStateDiagnostics(
            compatible=True,
            phospho_transformation_state={
                "kind": "log2",
                "transformed": True,
                "established_by": "test.phospho",
            },
            total_protein_transformation_state={
                "kind": "log2",
                "transformed": True,
                "established_by": "test.total",
            },
        ),
        preparation_policy="prepare_model_inputs",
        protein_mapping_policy="require_unambiguous",
        policy_parameters={
            "preparation_mode": "aligned_model_input_preparation_only",
            "modifies_phospho_matrix": False,
            "performs_total_protein_subtraction": False,
            "performs_differential_model_adjustment": False,
            "dataset_binding_table_fingerprints": [
                table_fingerprint_to_payload(fingerprint)
                for fingerprint in (
                    fingerprint_optional_table_strict(
                        phospho,
                        name="dataset.phospho",
                    ),
                    fingerprint_optional_table_strict(
                        site_metadata,
                        name="dataset.site_metadata",
                    ),
                    fingerprint_optional_table_strict(total, name="dataset.total"),
                )
                if fingerprint is not None
            ],
        },
    )
    protein_covariates = total.loc[
        list(dict.fromkeys(total_row_keys)),
        phospho.columns.tolist(),
    ].copy(deep=True)
    return ProteinAwarePreparationResult(
        matched_pairs=pd.DataFrame(
            {
                "site_key": list(site_keys),
                "protein_identifier": [
                    str(site_metadata.loc[site_key, "protein_identifier"])
                    for site_key in site_keys
                ],
                "total_protein_row_key": list(total_row_keys),
            }
        ),
        protein_covariate_matrix=protein_covariates,
        report=report,
    )


def _dataset_with_protein_aware_preparation() -> tuple[
    AnalysisReadyPhosphoDataset, ProteinAwarePreparationResult
]:
    return _dataset_with_custom_protein_aware_preparation()


def _dataset_with_custom_protein_aware_preparation(
    *,
    phospho: pd.DataFrame | None = None,
    total: pd.DataFrame | None = None,
) -> tuple[AnalysisReadyPhosphoDataset, ProteinAwarePreparationResult]:
    base_dataset = _dataset()
    resolved_phospho = base_dataset.phospho if phospho is None else phospho
    site_metadata = base_dataset.site_metadata
    site_keys = tuple(resolved_phospho.index.astype(str).tolist())
    total_row_keys = tuple(
        str(site_metadata.loc[site_key, "protein_identifier"]) for site_key in site_keys
    )
    resolved_total = (
        pd.DataFrame(
            {
                sample_id: [
                    1000.0 + float(position * 100 + sample_position)
                    for position in range(len(total_row_keys))
                ]
                for sample_position, sample_id in enumerate(
                    resolved_phospho.columns.astype(str)
                )
            },
            index=pd.Index(total_row_keys, name="protein_id"),
        )
        if total is None
        else total
    )
    base_dataset_with_total = trusted_analysis_ready_dataset_from_tables(
        phospho=resolved_phospho,
        site_metadata=base_dataset.site_metadata,
        sample_metadata=base_dataset.sample_metadata,
        total=resolved_total,
        comparisons=base_dataset.comparisons,
        organism=base_dataset.organism,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=True
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=True),
    )
    preparation = _protein_aware_preparation_for_dataset(base_dataset_with_total)
    preprocessing_report = DatasetPreprocessingReport.from_rows(
        protein_aware_preparation=preparation.report
    )
    return (
        trusted_analysis_ready_dataset_from_tables(
            phospho=base_dataset_with_total.phospho,
            site_metadata=base_dataset_with_total.site_metadata,
            sample_metadata=base_dataset_with_total.sample_metadata,
            total=base_dataset_with_total.total,
            comparisons=base_dataset_with_total.comparisons,
            organism=base_dataset_with_total.organism,
            intensity_scale_state=base_dataset_with_total.intensity_scale_state,
            processing_state=base_dataset_with_total.processing_state,
            preprocessing_report=preprocessing_report,
            protein_aware_preparation=preparation,
        ),
        preparation,
    )


def _six_sample_dataset_with_protein_aware_preparation() -> tuple[
    AnalysisReadyPhosphoDataset, ProteinAwarePreparationResult
]:
    base_dataset = _dataset()
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 1.0],
            "A_2": [1.1, 2.1, 1.1],
            "A_3": [0.9, 1.9, 0.8],
            "B_1": [2.1, 2.0, 1.0],
            "B_2": [2.0, 2.2, 0.9],
            "B_3": [2.2, 2.3, 1.2],
        },
        index=base_dataset.phospho.index.copy(),
    )
    total = pd.DataFrame(
        {
            "A_1": [10.0, 20.0, 30.0],
            "A_2": [11.0, 21.0, 30.5],
            "A_3": [9.0, 19.0, 32.0],
            "B_1": [13.0, 24.0, 29.0],
            "B_2": [12.0, 23.0, 31.0],
            "B_3": [14.0, 25.0, 33.0],
        },
        index=pd.Index(["MAPK14", "GSK3B", "AKT1"], name="protein_id"),
    )
    return _dataset_with_custom_protein_aware_preparation(
        phospho=phospho,
        total=total,
    )


def _protein_aware_request(
    dataset: AnalysisReadyPhosphoDataset,
    *,
    sample_ids: tuple[str, ...] = ("A_1", "A_2", "B_1", "B_2"),
    paired_design_policy: str = "reject",
    allow_design_subset: bool = False,
    fixed_effects: tuple[FixedEffectCovariate, ...] = (),
    sample_covariates: Mapping[str, Mapping[str, str | int | float]] | None = None,
) -> DifferentialAnalysisRequest:
    sample_blocks = {
        "A_1": "block_1",
        "B_1": "block_1",
        "A_2": "block_2",
        "B_2": "block_2",
        "A_3": "block_3",
        "B_3": "block_3",
    }
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=ExperimentalDesign(
            samples=tuple(
                SampleDesignRecord(
                    sample_id=sample_id,
                    condition=sample_id.split("_", maxsplit=1)[0],
                    biological_replicate_id=f"{sample_id}_bio",
                    block_id=(
                        sample_blocks[sample_id]
                        if paired_design_policy
                        in {PAIRED_DESIGN_POLICY_FIXED_BLOCK, "duplicate_correlation"}
                        else None
                    ),
                    covariates=(
                        sample_covariates.get(sample_id, {})
                        if sample_covariates is not None
                        else {}
                    ),
                )
                for sample_id in sample_ids
            ),
            fixed_effects=fixed_effects,
        ),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
        config=DifferentialAnalysisConfig(
            paired_design_policy=paired_design_policy,  # type: ignore[arg-type]
            allow_design_subset=allow_design_subset,
            protein_aware_model=DifferentialProteinAwareModelConfig(),
        ),
    )


def _validated_protein_aware_request(
    dataset: AnalysisReadyPhosphoDataset,
    *,
    sample_ids: tuple[str, ...] = ("A_1", "A_2", "B_1", "B_2"),
    paired_design_policy: str = "reject",
    allow_design_subset: bool = False,
    fixed_effects: tuple[FixedEffectCovariate, ...] = (),
    sample_covariates: Mapping[str, Mapping[str, str | int | float]] | None = None,
) -> ValidatedDifferentialAnalysisRequest:
    return DifferentialAnalysisValidator().run(
        _protein_aware_request(
            dataset,
            sample_ids=sample_ids,
            paired_design_policy=paired_design_policy,
            allow_design_subset=allow_design_subset,
            fixed_effects=fixed_effects,
            sample_covariates=sample_covariates,
        )
    )


def _interpreted_protein_aware_request(
    dataset: AnalysisReadyPhosphoDataset,
) -> InterpretedDifferentialAnalysisRequest:
    return DifferentialAnalysisInterpreter().run(
        _validated_protein_aware_request(dataset)
    )


def _preparation_report_with(
    report: ProteinAwarePreparationReport,
    *,
    site_eligibility: tuple[ProteinAwareSiteEligibility, ...] | None = None,
    transformation_state: ProteinAwareTransformationStateDiagnostics | None = None,
    replace_transformation_state: bool = False,
    preparation_policy: str | None = None,
    schema_version: int | None = None,
) -> ProteinAwarePreparationReport:
    return ProteinAwarePreparationReport(
        site_eligibility=(
            report.site_eligibility if site_eligibility is None else site_eligibility
        ),
        mapping_diagnostics=report.mapping_diagnostics,
        sample_alignment=report.sample_alignment,
        transformation_state=(
            transformation_state
            if replace_transformation_state
            else report.transformation_state
        ),
        preparation_policy=(
            report.preparation_policy
            if preparation_policy is None
            else preparation_policy
        ),
        protein_mapping_policy=report.protein_mapping_policy,
        policy_parameters=dict(report.policy_parameters),
        provenance=report.provenance,
        schema_version=report.schema_version
        if schema_version is None
        else schema_version,
    )


def _preparation_with(
    preparation: ProteinAwarePreparationResult,
    *,
    matched_pairs: pd.DataFrame | None = None,
    protein_covariate_matrix: pd.DataFrame | None = None,
    report: ProteinAwarePreparationReport | None = None,
) -> ProteinAwarePreparationResult:
    return ProteinAwarePreparationResult(
        matched_pairs=(
            preparation.matched_pairs_dataframe()
            if matched_pairs is None
            else matched_pairs
        ),
        protein_covariate_matrix=(
            preparation.protein_covariate_matrix_dataframe()
            if protein_covariate_matrix is None
            else protein_covariate_matrix
        ),
        report=preparation.report if report is None else report,
    )


def _replace_eligibility_row(
    report: ProteinAwarePreparationReport,
    *,
    position: int = 0,
    **changes: object,
) -> tuple[ProteinAwareSiteEligibility, ...]:
    rows = list(report.site_eligibility)
    rows[position] = replace(rows[position], **changes)
    return tuple(rows)


def _dataset_replacing_protein_aware_preparation(
    dataset: AnalysisReadyPhosphoDataset,
    preparation: ProteinAwarePreparationResult,
) -> AnalysisReadyPhosphoDataset:
    return trusted_analysis_ready_dataset_from_tables(
        phospho=dataset.phospho,
        site_metadata=dataset.site_metadata,
        sample_metadata=dataset.sample_metadata,
        total=dataset.total,
        comparisons=dataset.comparisons,
        organism=dataset.organism,
        intensity_scale_state=dataset.intensity_scale_state,
        processing_state=dataset.processing_state,
        preprocessing_report=DatasetPreprocessingReport.from_rows(
            protein_aware_preparation=preparation.report
        ),
        protein_aware_preparation=preparation,
    )


def _public_methods(cls: type[Any]) -> set[str]:
    return {
        name
        for name, value in inspect.getmembers(cls)
        if callable(value) and not name.startswith("_")
    }


def test_differential_workflow_calls_validator_interpreter_executor_in_order() -> None:
    events: list[str] = []
    validated = object()
    interpreted = object()
    expected_result = object()

    class _Validator:
        def run(self, request: object) -> object:
            events.append("validator")
            return validated

    class _Interpreter:
        def run(self, request: object) -> object:
            events.append("interpreter")
            assert request is validated
            return interpreted

    class _Executor:
        def run(self, request: object) -> object:
            events.append("executor")
            assert request is interpreted
            return expected_result

    result = DifferentialAnalysisWorkflow._with_components(
        validator=_Validator(),  # type: ignore[arg-type]
        interpreter=_Interpreter(),  # type: ignore[arg-type]
        executor=_Executor(),  # type: ignore[arg-type]
    ).run(_request())

    assert events == ["validator", "interpreter", "executor"]
    assert result is expected_result


def test_differential_workflow_dependency_injection_supports_real_stage_contracts() -> (
    None
):
    workflow = DifferentialAnalysisWorkflow._with_components(
        validator=DifferentialAnalysisValidator(),
        interpreter=DifferentialAnalysisInterpreter(),
        executor=DifferentialAnalysisExecutor(),
    )
    result = workflow.run(_request())
    assert isinstance(result, DifferentialAnalysisResult)


def test_differential_workflow_uses_explicit_design_not_sample_metadata_conditions() -> (
    None
):
    base_request = _request()
    base_dataset = base_request.dataset
    dataset_with_passive_metadata = trusted_analysis_ready_dataset_from_tables(
        phospho=base_dataset.phospho,
        site_metadata=base_dataset.site_metadata,
        sample_metadata=pd.DataFrame(
            {
                "condition": ["metadata_only"] * 4,
                "batch": ["batch_1", "batch_1", "batch_2", "batch_2"],
            },
            index=base_dataset.phospho.columns.copy(),
        ),
        organism=base_dataset.organism,
        intensity_scale_state=base_dataset.intensity_scale_state,
        processing_state=base_dataset.processing_state,
    )

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset_with_passive_metadata,
            design=base_request.design,
            contrasts=base_request.contrasts,
        )
    )

    assert isinstance(result, DifferentialAnalysisResult)


def test_differential_request_has_no_protein_aware_preparation_channel() -> None:
    assert [field.name for field in fields(DifferentialAnalysisRequest)] == [
        "dataset",
        "design",
        "contrasts",
        "config",
    ]
    request = _request()
    forbidden_inputs = {
        "protein_aware_preparation": object(),
        "matched_pairs": pd.DataFrame(
            {"site_key": ["site_a"], "total_protein_row_key": ["protein_a"]}
        ),
        "protein_covariate_matrix": pd.DataFrame(
            {"sample_a": [1.0]},
            index=pd.Index(["protein_a"], name="total_protein_row_key"),
        ),
        "protein_matrix": pd.DataFrame(
            {"sample_a": [1.0]},
            index=pd.Index(["protein_a"], name="total_protein_row_key"),
        ),
    }
    request_type = cast(Any, DifferentialAnalysisRequest)

    for keyword, value in forbidden_inputs.items():
        with pytest.raises(TypeError, match="unexpected keyword argument"):
            request_type(
                dataset=request.dataset,
                design=request.design,
                contrasts=request.contrasts,
                **{keyword: value},
            )


def test_ordinary_differential_lane_ignores_dataset_owned_protein_aware_sidecar() -> (
    None
):
    dataset_with_preparation, preparation = _dataset_with_protein_aware_preparation()
    dataset_without_preparation = trusted_analysis_ready_dataset_from_tables(
        phospho=dataset_with_preparation.phospho,
        site_metadata=dataset_with_preparation.site_metadata,
        sample_metadata=dataset_with_preparation.sample_metadata,
        total=dataset_with_preparation.total,
        comparisons=dataset_with_preparation.comparisons,
        organism=dataset_with_preparation.organism,
        intensity_scale_state=dataset_with_preparation.intensity_scale_state,
        processing_state=dataset_with_preparation.processing_state,
    )
    base_request = _request()

    absent_result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset_without_preparation,
            design=base_request.design,
            contrasts=base_request.contrasts,
        )
    )
    present_result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset_with_preparation,
            design=base_request.design,
            contrasts=base_request.contrasts,
        )
    )

    assert present_result.to_payload() == absent_result.to_payload()
    assert present_result.protein_aware_diagnostics is None
    assert absent_result.protein_aware_diagnostics is None
    assert not hasattr(present_result, "protein_aware_preparation")
    assert not hasattr(present_result, "protein_covariate_matrix")
    assert present_result.workflow_provenance is not None
    assert present_result.workflow_provenance["input_intensity_scale"] == "log2"
    assert (
        present_result.workflow_provenance["input_intensity_scale_evidence_level"]
        == "declared_by_user"
    )
    assert (
        present_result.workflow_provenance["input_intensity_scale_source"]
        == "declared_by_user"
    )
    assert (
        present_result.input_dataset_preprocessing_report is not None
        and present_result.input_dataset_preprocessing_report.protein_aware_preparation
        is preparation.report
    )
    assert absent_result.input_dataset_preprocessing_report is None
    assert "protein_aware" not in repr(present_result.policy_provenance).lower()


def test_protein_aware_public_workflow_branch_consumes_dataset_owned_sidecar() -> None:
    dataset, preparation = _dataset_with_protein_aware_preparation()
    withheld_site = str(dataset.phospho.index[1])
    fallback_rows = _replace_eligibility_row(
        preparation.report,
        position=1,
        eligibility=ProteinAwarePreparationEligibility.FALLBACK_TO_PHOSPHO_ONLY,
        mapping_status=ProteinMappingStatus.MISSING_TOTAL_PROTEIN_ROW,
        total_protein_row_key=None,
        reasons=(PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW,),
    )
    matched_pairs = preparation.matched_pairs_dataframe()
    matched_pairs = (
        matched_pairs.loc[matched_pairs.loc[:, "site_key"] != withheld_site, :]
        .reset_index(drop=True)
        .copy(deep=True)
    )
    protein_covariates = preparation.protein_covariate_matrix_dataframe().drop(
        index="GSK3B"
    )
    modified_dataset = _dataset_replacing_protein_aware_preparation(
        dataset,
        _preparation_with(
            preparation,
            matched_pairs=matched_pairs,
            protein_covariate_matrix=protein_covariates,
            report=_preparation_report_with(
                preparation.report,
                site_eligibility=fallback_rows,
            ),
        ),
    )

    result = DifferentialAnalysisWorkflow().run(
        _protein_aware_request(modified_dataset)
    )

    table = result.table_for("B_vs_A")
    diagnostics = result.protein_aware_diagnostics
    assert diagnostics is not None
    assert diagnostics.tested_site_count == 2
    assert diagnostics.withheld_site_count == 1
    assert result.policy_provenance is not None
    assert result.policy_provenance.protein_aware is not None
    assert table.loc[withheld_site, DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE
    )
    assert table.loc[withheld_site, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] == (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK
    )
    assert pd.isna(table.loc[withheld_site, "logFC"])


def test_protein_aware_public_workflow_branch_runs_end_to_end() -> None:
    dataset, _ = _dataset_with_protein_aware_preparation()

    result = DifferentialAnalysisWorkflow().run(_protein_aware_request(dataset))

    diagnostics = result.protein_aware_diagnostics
    assert diagnostics is not None
    assert diagnostics.method_id == DifferentialProteinAwareModelConfig().method
    assert diagnostics.tested_site_count == 3
    assert diagnostics.withheld_site_count == 0
    assert diagnostics.execution_sample_order == ("A_1", "A_2", "B_1", "B_2")
    assert diagnostics.protein_covariate_centered is True
    assert diagnostics.protein_covariate_standardized is False
    assert result.policy_provenance is not None
    assert result.policy_provenance.protein_aware is not None

    table = result.table_for("B_vs_A")
    assert (
        table[DIFFERENTIAL_RESULT_STATUS_COLUMN].tolist()
        == [DIFFERENTIAL_RESULT_STATUS_TESTED] * 3
    )
    assert table.loc[:, ["logFC", "t", "P.Value", "adj.P.Val"]].notna().all().all()


def test_interpreted_request_rejects_inconsistent_protein_aware_boundary_state() -> (
    None
):
    dataset, _ = _dataset_with_protein_aware_preparation()
    interpreted = _interpreted_protein_aware_request(dataset)
    resolved_inputs = interpreted.protein_aware_inputs
    assert resolved_inputs is not None

    ordinary_execution_config = replace(
        interpreted.execution_config,
        protein_aware_method=None,
    )
    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.interpreter.unselected_protein_aware_inputs",
    ):
        replace(interpreted, execution_config=ordinary_execution_config)

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.interpreter.protein_aware_inputs",
    ):
        replace(interpreted, protein_aware_inputs=None)

    object.__setattr__(resolved_inputs, "method_id", "mismatched_method_for_test")
    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.interpreter.protein_aware_method",
    ) as exc_info:
        replace(interpreted, protein_aware_inputs=resolved_inputs)

    assert exc_info.value.details == {
        "execution_method": DifferentialProteinAwareModelConfig().method,
        "resolved_method": "mismatched_method_for_test",
    }


def test_protein_aware_executor_translates_kernel_input_errors_with_diagnostics() -> (
    None
):
    dataset, _ = _dataset_with_protein_aware_preparation()
    interpreted = _interpreted_protein_aware_request(dataset)

    class _FailingProteinAwareKernel:
        def run(self, request: object) -> object:
            raise PhosPyInputError(
                "protein-aware kernel rejected trusted test input",
                diagnostics={
                    "site_failure_diagnostics": pd.DataFrame(
                        {
                            DIFFERENTIAL_RESULT_STATUS_COLUMN: [
                                DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID
                            ],
                            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: [
                                DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE
                            ],
                        },
                        index=pd.Index(["site_a"], name="site_key"),
                    ),
                    "augmented_design_failure_diagnostics": pd.DataFrame(
                        {"failure_message": ["zero-variance protein covariate"]},
                        index=pd.Index(["protein_a"], name="total_protein_row_key"),
                    ),
                },
            )

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.executor.protein_aware_fit",
    ) as exc_info:
        DifferentialAnalysisExecutor(
            protein_aware_kernel=_FailingProteinAwareKernel(),  # type: ignore[arg-type]
        ).run(interpreted)

    assert exc_info.value.next_action is not None
    assert "protein-aware eligible phosphosite" in exc_info.value.next_action
    assert exc_info.value.details == {
        "error": "protein-aware kernel rejected trusted test input",
        "diagnostics_type": "dict",
        "diagnostic_keys": (
            "site_failure_diagnostics",
            "augmented_design_failure_diagnostics",
        ),
        "diagnostic_row_counts": {
            "site_failure_diagnostics": 1,
            "augmented_design_failure_diagnostics": 1,
        },
    }


def test_protein_aware_input_resolver_builds_execution_request_from_sidecar() -> None:
    dataset, _ = _dataset_with_protein_aware_preparation()
    sample_ids = ("B_2", "B_1", "A_2", "A_1")

    resolved = ProteinAwareDifferentialInputResolver().run(
        _validated_protein_aware_request(
            dataset,
            sample_ids=sample_ids,
            allow_design_subset=True,
        )
    )

    assert isinstance(resolved, ProteinAwareDifferentialResolvedInputs)
    assert resolved.sample_order == sample_ids
    assert resolved.full_site_ids == tuple(dataset.phospho.index.astype(str).tolist())
    assert resolved.tested_site_ids == resolved.full_site_ids
    assert resolved.status_counts == ((DIFFERENTIAL_RESULT_STATUS_TESTED, 3),)
    assert dict(resolved.eligibility_counts) == {
        "total_site_count": 3,
        "ordinary_testable_site_count": 3,
        "protein_preparation_candidate_site_count": 3,
        "protein_aware_tested_site_count": 3,
    }
    assert list(resolved.resolved_protein_covariates.columns.astype(str)) == list(
        sample_ids
    )
    assert list(
        resolved.computation_request.phosphosite_matrix.columns.astype(str)
    ) == list(sample_ids)
    assert (
        tuple(
            resolved.computation_request.matched_pairs.loc[:, "site_key"]
            .astype(str)
            .tolist()
        )
        == resolved.tested_site_ids
    )
    assert (
        tuple(resolved.candidate_matched_pairs.loc[:, "site_key"].astype(str).tolist())
        == resolved.full_site_ids
    )
    assert resolved.feature_eligibility_inputs.attach_to_result_tables is True
    assert (
        resolved.site_eligibility_metadata.loc[
            :, "protein_aware_centering_policy"
        ].unique()
        == ["mean_centered_no_standardization"]
    ).all()


def test_protein_aware_input_resolver_uses_analysis_sample_subset_only() -> None:
    dataset, _ = _six_sample_dataset_with_protein_aware_preparation()
    sample_ids = ("B_2", "A_2", "B_1", "A_1")

    resolved = ProteinAwareDifferentialInputResolver().run(
        _validated_protein_aware_request(
            dataset,
            sample_ids=sample_ids,
            allow_design_subset=True,
        )
    )

    assert tuple(resolved.resolved_protein_covariates.columns.astype(str)) == sample_ids
    assert (
        tuple(resolved.computation_request.resolved_protein_covariates.columns)
        == sample_ids
    )
    assert tuple(resolved.computation_request.phosphosite_matrix.columns) == sample_ids
    assert "A_3" not in resolved.resolved_protein_covariates.columns
    assert "B_3" not in resolved.resolved_protein_covariates.columns


def test_protein_aware_input_resolver_rejects_missing_sidecar() -> None:
    validated = _validated_protein_aware_request(_dataset())

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.sidecar_missing",
    ):
        ProteinAwareDifferentialInputResolver().run(validated)


def test_protein_aware_input_resolver_rejects_unsupported_sidecar_schema() -> None:
    dataset, preparation = _dataset_with_protein_aware_preparation()
    object.__setattr__(preparation.report, "schema_version", 99)

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.sidecar_schema",
    ):
        ProteinAwareDifferentialInputResolver().run(
            _validated_protein_aware_request(dataset)
        )


def test_protein_aware_input_resolver_rejects_unsupported_sidecar_policy() -> None:
    dataset, preparation = _dataset_with_protein_aware_preparation()
    object.__setattr__(preparation.report, "preparation_policy", "legacy_policy")

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.sidecar_policy",
    ):
        ProteinAwareDifferentialInputResolver().run(
            _validated_protein_aware_request(dataset)
        )


def test_protein_aware_input_resolver_rejects_missing_total_matrix() -> None:
    dataset, _ = _dataset_with_protein_aware_preparation()
    validated = _validated_protein_aware_request(dataset)
    unsafe_remove_dataset_total_matrix(validated.dataset)

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.total_matrix_missing",
    ):
        ProteinAwareDifferentialInputResolver().run(
            replace(validated, dataset_view=DatasetInternalView(validated.dataset))
        )


def test_protein_aware_input_resolver_rejects_non_log2_total_scale() -> None:
    dataset, _ = _dataset_with_protein_aware_preparation()
    validated = _validated_protein_aware_request(dataset)
    unsafe_replace_dataset_total_scale_state(
        validated.dataset,
        MatrixIntensityScaleState.linear(established_by="test.linear_total"),
    )

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.total_scale",
    ):
        ProteinAwareDifferentialInputResolver().run(validated)


def test_protein_aware_input_resolver_rejects_non_log2_phosphosite_scale() -> None:
    dataset, _ = _dataset_with_protein_aware_preparation()
    validated = _validated_protein_aware_request(dataset)
    unsafe_replace_dataset_phospho_scale_state(
        validated.dataset,
        MatrixIntensityScaleState.linear(established_by="test.linear_phospho"),
    )

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.phospho_scale",
    ):
        ProteinAwareDifferentialInputResolver().run(validated)


def test_protein_aware_input_resolver_rejects_non_log2_sidecar_covariate_evidence() -> (
    None
):
    dataset, preparation = _dataset_with_protein_aware_preparation()
    report = _preparation_report_with(
        preparation.report,
        transformation_state=ProteinAwareTransformationStateDiagnostics(
            compatible=True,
            phospho_transformation_state={
                "kind": "log2",
                "transformed": True,
                "established_by": "test.phospho",
            },
            total_protein_transformation_state={
                "kind": "linear",
                "transformed": False,
                "established_by": "test.linear_total",
            },
        ),
        replace_transformation_state=True,
    )
    modified_dataset = _dataset_replacing_protein_aware_preparation(
        dataset,
        _preparation_with(preparation, report=report),
    )

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.prepared_covariate_scale",
    ) as exc_info:
        ProteinAwareDifferentialInputResolver().run(
            _validated_protein_aware_request(modified_dataset)
        )

    assert (
        exc_info.value.details["total_protein_transformation_state"]["kind"] == "linear"
    )


def test_protein_aware_input_resolver_rejects_unestablished_sidecar_covariate_evidence() -> (
    None
):
    dataset, preparation = _dataset_with_protein_aware_preparation()
    report = _preparation_report_with(
        preparation.report,
        transformation_state=ProteinAwareTransformationStateDiagnostics(
            compatible=True,
            phospho_transformation_state={
                "kind": "log2",
                "transformed": True,
                "established_by": "test.phospho",
            },
            total_protein_transformation_state={
                "kind": "log2",
                "transformed": True,
                "established_by": "",
            },
        ),
        replace_transformation_state=True,
    )
    modified_dataset = _dataset_replacing_protein_aware_preparation(
        dataset,
        _preparation_with(preparation, report=report),
    )

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.prepared_covariate_scale",
    ):
        ProteinAwareDifferentialInputResolver().run(
            _validated_protein_aware_request(modified_dataset)
        )


def test_protein_aware_input_resolver_rejects_non_log2_sidecar_phospho_evidence() -> (
    None
):
    dataset, preparation = _dataset_with_protein_aware_preparation()
    report = _preparation_report_with(
        preparation.report,
        transformation_state=ProteinAwareTransformationStateDiagnostics(
            compatible=True,
            phospho_transformation_state={
                "kind": "linear",
                "transformed": False,
                "established_by": "test.linear_phospho",
            },
            total_protein_transformation_state={
                "kind": "log2",
                "transformed": True,
                "established_by": "test.total",
            },
        ),
        replace_transformation_state=True,
    )
    modified_dataset = _dataset_replacing_protein_aware_preparation(
        dataset,
        _preparation_with(preparation, report=report),
    )

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.prepared_phospho_scale",
    ):
        ProteinAwareDifferentialInputResolver().run(
            _validated_protein_aware_request(modified_dataset)
        )


def test_protein_aware_input_resolver_rejects_prior_total_subtraction() -> None:
    dataset, _ = _dataset_with_protein_aware_preparation()
    validated = _validated_protein_aware_request(dataset)
    unsafe_mark_dataset_total_protein_correction_applied(validated.dataset)

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.prior_total_protein_subtraction",
    ):
        ProteinAwareDifferentialInputResolver().run(validated)


def test_protein_aware_input_resolver_rejects_stale_sidecar_binding() -> None:
    dataset, preparation = _dataset_with_protein_aware_preparation()
    mutated_covariates = preparation.protein_covariate_matrix_dataframe()
    mutated_covariates.iloc[0, 0] = mutated_covariates.iloc[0, 0] + 100.0
    object.__setattr__(preparation, "_protein_covariate_matrix", mutated_covariates)

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.sidecar_binding",
    ):
        ProteinAwareDifferentialInputResolver().run(
            _validated_protein_aware_request(dataset)
        )


def test_protein_aware_input_resolver_all_fallback_fails_with_status_counts() -> None:
    dataset, preparation = _dataset_with_protein_aware_preparation()
    fallback_rows = tuple(
        replace(
            row,
            eligibility=ProteinAwarePreparationEligibility.FALLBACK_TO_PHOSPHO_ONLY,
            mapping_status=ProteinMappingStatus.MISSING_TOTAL_PROTEIN_ROW,
            total_protein_row_key=None,
            reasons=(PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW,),
        )
        for row in preparation.report.site_eligibility
    )
    modified_preparation = _preparation_with(
        preparation,
        matched_pairs=preparation.matched_pairs_dataframe().iloc[0:0, :],
        protein_covariate_matrix=preparation.protein_covariate_matrix_dataframe().iloc[
            0:0, :
        ],
        report=_preparation_report_with(
            preparation.report,
            site_eligibility=fallback_rows,
        ),
    )
    modified_dataset = _dataset_replacing_protein_aware_preparation(
        dataset,
        modified_preparation,
    )

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.all_sites_withheld",
    ) as exc_info:
        ProteinAwareDifferentialInputResolver().run(
            _validated_protein_aware_request(modified_dataset)
        )

    assert exc_info.value.details["status_counts"] == {
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE: 3
    }
    assert exc_info.value.details["reason_counts"] == {
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK: 3
    }
    assert exc_info.value.details["eligibility_counts"] == {
        "total_site_count": 3,
        "ordinary_testable_site_count": 3,
        "protein_preparation_candidate_site_count": 0,
        "protein_aware_tested_site_count": 0,
    }


def test_protein_aware_input_resolver_rejects_duplicate_correlation() -> None:
    sample_ids = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
    dataset, _ = _six_sample_dataset_with_protein_aware_preparation()

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.duplicate_correlation",
    ):
        ProteinAwareDifferentialInputResolver().run(
            _validated_protein_aware_request(
                dataset,
                sample_ids=sample_ids,
                paired_design_policy="duplicate_correlation",
            )
        )


def test_protein_aware_input_resolver_rejects_actual_technical_aggregation() -> None:
    dataset, _ = _dataset_with_protein_aware_preparation()
    validated = _validated_protein_aware_request(dataset)
    plan = TechnicalReplicateAggregationPlan(
        technical_replicate_policy=TechnicalReplicatePolicy.MEAN,
        groups=(
            TechnicalReplicateAggregationGroup(
                condition="A",
                biological_replicate_id="A_r1",
                output_sample_id="A_r1",
                input_sample_ids=("A_1", "A_2"),
                technical_replicate_ids=("tech_1", "tech_2"),
                batch=None,
                block_id=None,
                covariates={},
            ),
        ),
        aggregate_phospho=True,
        aggregate_total_protein=True,
    )

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.technical_replicate_aggregation",
    ):
        ProteinAwareDifferentialInputResolver().run(
            replace(validated, technical_replicate_aggregation_plan=plan)
        )


def test_protein_aware_input_resolver_accepts_noop_technical_aggregation_plan() -> None:
    dataset, _ = _dataset_with_protein_aware_preparation()
    validated = _validated_protein_aware_request(dataset)
    plan = TechnicalReplicateAggregationPlan(
        technical_replicate_policy=TechnicalReplicatePolicy.MEAN,
        groups=(),
        aggregate_phospho=True,
        aggregate_total_protein=True,
    )

    resolved = ProteinAwareDifferentialInputResolver().run(
        replace(validated, technical_replicate_aggregation_plan=plan)
    )

    assert resolved.tested_site_ids == resolved.full_site_ids
    assert resolved.status_counts == ((DIFFERENTIAL_RESULT_STATUS_TESTED, 3),)


def test_protein_aware_input_resolver_marks_preparation_fallback_rows() -> None:
    dataset, preparation = _dataset_with_protein_aware_preparation()
    fallback_site = str(dataset.phospho.index[1])
    fallback_rows = _replace_eligibility_row(
        preparation.report,
        position=1,
        eligibility=ProteinAwarePreparationEligibility.FALLBACK_TO_PHOSPHO_ONLY,
        mapping_status=ProteinMappingStatus.MISSING_TOTAL_PROTEIN_ROW,
        total_protein_row_key=None,
        reasons=(PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW,),
    )
    report = _preparation_report_with(
        preparation.report,
        site_eligibility=fallback_rows,
    )
    matched_pairs = preparation.matched_pairs_dataframe()
    matched_pairs = (
        matched_pairs.loc[matched_pairs.loc[:, "site_key"] != fallback_site, :]
        .reset_index(drop=True)
        .copy(deep=True)
    )
    protein_covariates = preparation.protein_covariate_matrix_dataframe()
    protein_covariates = protein_covariates.drop(index="GSK3B")
    modified_preparation = _preparation_with(
        preparation,
        matched_pairs=matched_pairs,
        protein_covariate_matrix=protein_covariates,
        report=report,
    )
    modified_dataset = _dataset_replacing_protein_aware_preparation(
        dataset,
        modified_preparation,
    )

    resolved = ProteinAwareDifferentialInputResolver().run(
        _validated_protein_aware_request(modified_dataset)
    )

    metadata = resolved.site_eligibility_metadata
    assert metadata.loc[fallback_site, DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE
    )
    assert metadata.loc[fallback_site, "result_status_reason"] == (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK
    )
    assert fallback_site in resolved.full_site_ids
    assert fallback_site not in resolved.tested_site_ids
    assert tuple(
        resolved.computation_request.matched_pairs.loc[:, "site_key"]
        .astype(str)
        .tolist()
    ) == tuple(site for site in resolved.full_site_ids if site != fallback_site)


def test_protein_aware_input_resolver_marks_unmatched_and_ambiguous_rows() -> None:
    dataset, preparation = _dataset_with_protein_aware_preparation()
    unmatched_site = str(dataset.phospho.index[1])
    ambiguous_site = str(dataset.phospho.index[2])
    rows = list(preparation.report.site_eligibility)
    rows[1] = replace(
        rows[1],
        eligibility=ProteinAwarePreparationEligibility.FALLBACK_TO_PHOSPHO_ONLY,
        mapping_status=ProteinMappingStatus.MISSING_TOTAL_PROTEIN_ROW,
        total_protein_row_key=None,
        reasons=(PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW,),
    )
    rows[2] = replace(
        rows[2],
        eligibility=ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION,
        mapping_status=ProteinMappingStatus.AMBIGUOUS_TOTAL_PROTEIN_MAPPING,
        total_protein_row_key=None,
        reasons=(PROTEIN_AWARE_REASON_AMBIGUOUS_PROTEIN_MAPPING,),
    )
    matched_pairs = preparation.matched_pairs_dataframe()
    matched_pairs = matched_pairs.loc[
        ~matched_pairs.loc[:, "site_key"].isin({unmatched_site, ambiguous_site}),
        :,
    ].reset_index(drop=True)
    protein_covariates = preparation.protein_covariate_matrix_dataframe()
    protein_covariates = protein_covariates.drop(index=["GSK3B", "AKT1"])
    modified_preparation = _preparation_with(
        preparation,
        matched_pairs=matched_pairs,
        protein_covariate_matrix=protein_covariates,
        report=_preparation_report_with(
            preparation.report,
            site_eligibility=tuple(rows),
        ),
    )
    modified_dataset = _dataset_replacing_protein_aware_preparation(
        dataset,
        modified_preparation,
    )

    resolved = ProteinAwareDifferentialInputResolver().run(
        _validated_protein_aware_request(modified_dataset)
    )

    metadata = resolved.site_eligibility_metadata
    assert metadata.loc[unmatched_site, DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE
    )
    assert metadata.loc[unmatched_site, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] == (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK
    )
    assert metadata.loc[ambiguous_site, DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE
    )
    assert metadata.loc[ambiguous_site, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] == (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED
    )
    assert metadata.loc[unmatched_site, "protein_aware_preparation_reasons"] == (
        PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW,
    )
    assert metadata.loc[ambiguous_site, "protein_aware_preparation_reasons"] == (
        PROTEIN_AWARE_REASON_AMBIGUOUS_PROTEIN_MAPPING,
    )
    assert resolved.tested_site_ids == (str(dataset.phospho.index[0]),)
    assert (
        tuple(
            resolved.computation_request.matched_pairs.loc[:, "site_key"]
            .astype(str)
            .tolist()
        )
        == resolved.tested_site_ids
    )


def test_protein_aware_input_resolver_preserves_ordinary_numeric_precedence() -> None:
    base_dataset = _dataset()
    phospho = base_dataset.phospho
    ordinary_withheld_site = str(phospho.index[0])
    phospho.loc[ordinary_withheld_site, :] = 42.0
    dataset, _ = _dataset_with_custom_protein_aware_preparation(phospho=phospho)

    resolved = ProteinAwareDifferentialInputResolver().run(
        _validated_protein_aware_request(dataset)
    )

    metadata = resolved.site_eligibility_metadata
    assert (
        metadata.loc[ordinary_withheld_site, DIFFERENTIAL_RESULT_STATUS_COLUMN]
        == DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT
    )
    assert (
        bool(metadata.loc[ordinary_withheld_site, "protein_aware_candidate"]) is False
    )
    assert ordinary_withheld_site in resolved.full_site_ids
    assert ordinary_withheld_site not in resolved.tested_site_ids


def test_protein_aware_input_resolver_marks_non_finite_protein_covariates() -> None:
    dataset, preparation = _six_sample_dataset_with_protein_aware_preparation()
    non_finite_covariates = preparation.protein_covariate_matrix_dataframe()
    non_finite_covariates.loc["MAPK14", "A_2"] = np.nan
    non_finite_preparation = _preparation_with(
        preparation,
        protein_covariate_matrix=non_finite_covariates,
    )
    non_finite_site = str(dataset.phospho.index[0])

    class _TrustedMalformedDatasetView:
        def __init__(
            self,
            dataset: AnalysisReadyPhosphoDataset,
            preparation: ProteinAwarePreparationResult,
        ) -> None:
            self._delegate = DatasetInternalView(dataset)
            self._preparation = preparation

        @property
        def phospho(self) -> pd.DataFrame:
            return self._delegate.phospho

        @property
        def total(self) -> pd.DataFrame | None:
            return self._delegate.total

        @property
        def protein_aware_preparation(self) -> ProteinAwarePreparationInternalView:
            return ProteinAwarePreparationInternalView(self._preparation)

        def validate_protein_aware_preparation_binding(self) -> None:
            return None

    validated = _validated_protein_aware_request(
        dataset,
        sample_ids=("A_1", "A_2", "A_3", "B_1", "B_2", "B_3"),
    )

    resolved = ProteinAwareDifferentialInputResolver().run(
        replace(
            validated,
            dataset_view=cast(
                Any,
                _TrustedMalformedDatasetView(dataset, non_finite_preparation),
            ),
        )
    )

    metadata = resolved.site_eligibility_metadata
    assert metadata.loc[non_finite_site, DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID
    )
    assert metadata.loc[non_finite_site, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] == (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE
    )
    assert non_finite_site not in resolved.tested_site_ids
    assert len(resolved.tested_site_ids) == 2


def test_protein_aware_input_resolver_marks_rank_deficient_augmented_designs() -> None:
    dataset, _ = _six_sample_dataset_with_protein_aware_preparation()
    total = dataset.total
    if total is None:
        raise AssertionError("protein-aware test dataset requires total")
    rank_deficient_total = total.copy(deep=True)
    rank_deficient_total.loc["MAPK14", :] = [10.0, 10.0, 10.0, 20.0, 20.0, 20.0]
    rank_deficient_dataset, _ = _dataset_with_custom_protein_aware_preparation(
        phospho=dataset.phospho, total=rank_deficient_total
    )
    rank_deficient_site = str(rank_deficient_dataset.phospho.index[0])

    resolved = ProteinAwareDifferentialInputResolver().run(
        _validated_protein_aware_request(
            rank_deficient_dataset,
            sample_ids=("A_1", "A_2", "A_3", "B_1", "B_2", "B_3"),
        )
    )

    metadata = resolved.site_eligibility_metadata
    assert metadata.loc[rank_deficient_site, DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID
    )
    assert metadata.loc[
        rank_deficient_site,
        DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    ] == (DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT)
    assert rank_deficient_site not in resolved.tested_site_ids
    assert len(resolved.tested_site_ids) == 2


def test_protein_aware_input_resolver_marks_ill_conditioned_augmented_designs() -> None:
    dataset, _ = _six_sample_dataset_with_protein_aware_preparation()
    total = dataset.total
    if total is None:
        raise AssertionError("protein-aware test dataset requires total")
    ill_conditioned_total = total.copy(deep=True)
    ill_conditioned_total.loc["MAPK14", :] = [
        10.0,
        10.0 + 1.0e-12,
        10.0 + 2.0e-12,
        20.0,
        20.0 + 1.0e-12,
        20.0 + 2.0e-12,
    ]
    ill_conditioned_dataset, _ = _dataset_with_custom_protein_aware_preparation(
        phospho=dataset.phospho, total=ill_conditioned_total
    )
    ill_conditioned_site = str(ill_conditioned_dataset.phospho.index[0])

    resolved = ProteinAwareDifferentialInputResolver().run(
        _validated_protein_aware_request(
            ill_conditioned_dataset,
            sample_ids=("A_1", "A_2", "A_3", "B_1", "B_2", "B_3"),
        )
    )

    metadata = resolved.site_eligibility_metadata
    assert metadata.loc[ill_conditioned_site, DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID
    )
    assert metadata.loc[
        ill_conditioned_site,
        DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    ] == (DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED)
    assert ill_conditioned_site not in resolved.tested_site_ids
    assert len(resolved.tested_site_ids) == 2


def test_protein_aware_input_resolver_allows_fixed_covariate_when_augmented_design_valid() -> (
    None
):
    sample_ids = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
    dataset, _ = _six_sample_dataset_with_protein_aware_preparation()
    sample_covariates = {
        "A_1": {"dose": 0.0},
        "A_2": {"dose": 1.0},
        "A_3": {"dose": 2.0},
        "B_1": {"dose": 0.5},
        "B_2": {"dose": 1.5},
        "B_3": {"dose": 2.5},
    }

    resolved = ProteinAwareDifferentialInputResolver().run(
        _validated_protein_aware_request(
            dataset,
            sample_ids=sample_ids,
            fixed_effects=(ContinuousCovariate("dose"),),
            sample_covariates=sample_covariates,
        )
    )

    assert resolved.tested_site_ids == tuple(dataset.phospho.index.astype(str).tolist())
    assert "dose" in resolved.base_design.frame.columns
    assert (
        resolved.site_eligibility_metadata.loc[
            :, "protein_augmented_design_residual_degrees_of_freedom"
        ]
        > 0.0
    ).all()


def test_protein_aware_input_resolver_does_not_mutate_sidecar_or_input_frames() -> None:
    dataset, preparation = _dataset_with_protein_aware_preparation()
    phospho_before = dataset.phospho
    total_before = dataset.total
    if total_before is None:
        raise AssertionError("protein-aware test dataset requires total")
    matched_pairs_before = preparation.matched_pairs_dataframe()
    covariates_before = preparation.protein_covariate_matrix_dataframe()
    site_eligibility_before = preparation.site_eligibility_table

    resolved = ProteinAwareDifferentialInputResolver().run(
        _validated_protein_aware_request(dataset)
    )

    pd.testing.assert_frame_equal(dataset.phospho, phospho_before)
    pd.testing.assert_frame_equal(dataset.total, total_before)
    pd.testing.assert_frame_equal(
        preparation.matched_pairs_dataframe(),
        matched_pairs_before,
    )
    pd.testing.assert_frame_equal(
        preparation.protein_covariate_matrix_dataframe(),
        covariates_before,
    )
    pd.testing.assert_frame_equal(
        preparation.site_eligibility_table,
        site_eligibility_before,
    )

    resolved.resolved_protein_covariates.iloc[0, 0] = -999.0
    resolved.matched_pairs.iloc[0, 0] = "mutated"
    assert (
        resolved.computation_request.resolved_protein_covariates.iloc[0, 0]
        == covariates_before.iloc[0, 0]
    )
    assert (
        resolved.computation_request.matched_pairs.iloc[0, 0]
        == (matched_pairs_before.iloc[0, 0])
    )
    pd.testing.assert_frame_equal(
        preparation.protein_covariate_matrix_dataframe(),
        covariates_before,
    )


def test_protein_aware_resolved_inputs_reject_trusted_construction_mismatch() -> None:
    dataset, _ = _dataset_with_protein_aware_preparation()
    resolved = ProteinAwareDifferentialInputResolver().run(
        _validated_protein_aware_request(dataset)
    )

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.tested_site_order",
    ):
        ProteinAwareDifferentialResolvedInputs(
            computation_request=resolved.computation_request,
            feature_eligibility_inputs=resolved.feature_eligibility_inputs,
            matched_pairs=resolved.matched_pairs,
            candidate_matched_pairs=resolved.candidate_matched_pairs,
            resolved_protein_covariates=resolved.resolved_protein_covariates,
            site_eligibility_metadata=resolved.site_eligibility_metadata,
            full_site_ids=resolved.full_site_ids,
            tested_site_ids=tuple(reversed(resolved.tested_site_ids)),
            sample_order=resolved.sample_order,
            base_design=resolved.base_design,
            base_contrasts=resolved.base_contrasts,
            method_id=resolved.method_id,
            preparation_policy=resolved.preparation_policy,
            protein_mapping_policy=resolved.protein_mapping_policy,
            eligibility_counts=resolved.eligibility_counts,
            status_counts=resolved.status_counts,
            reason_counts=resolved.reason_counts,
        )


def test_protein_aware_input_resolver_all_withheld_fails_with_counts() -> None:
    dataset, _ = _dataset_with_protein_aware_preparation()
    total = dataset.total
    if total is None:
        raise AssertionError("protein-aware test dataset requires total")
    constant_total = total.copy(deep=True)
    for row_key, value in zip(constant_total.index, (10.0, 20.0, 30.0), strict=True):
        constant_total.loc[row_key, :] = value
    constant_dataset, _ = _dataset_with_custom_protein_aware_preparation(
        total=constant_total
    )

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.all_sites_withheld",
    ) as exc_info:
        ProteinAwareDifferentialInputResolver().run(
            _validated_protein_aware_request(constant_dataset)
        )

    assert exc_info.value.details["status_counts"] == {
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID: 3
    }
    assert exc_info.value.details["reason_counts"] == {
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE: 3
    }
    assert exc_info.value.details["eligibility_counts"] == {
        "total_site_count": 3,
        "ordinary_testable_site_count": 3,
        "protein_preparation_candidate_site_count": 3,
        "protein_aware_tested_site_count": 0,
    }


def test_protein_aware_input_resolver_allows_fixed_block_when_augmented_design_valid() -> (
    None
):
    sample_ids = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
    dataset, _ = _six_sample_dataset_with_protein_aware_preparation()

    resolved = ProteinAwareDifferentialInputResolver().run(
        _validated_protein_aware_request(
            dataset,
            sample_ids=sample_ids,
            paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
        )
    )

    assert resolved.tested_site_ids == tuple(dataset.phospho.index.astype(str).tolist())
    assert resolved.computation_request.base_design.frame.shape[1] == 4
    assert (
        resolved.site_eligibility_metadata.loc[
            :, "protein_augmented_design_residual_degrees_of_freedom"
        ]
        > 0.0
    ).all()


def test_differential_result_references_input_dataset_preprocessing_report() -> None:
    base_dataset = _dataset()
    preprocessing_report = DatasetPreprocessingReport.from_rows()
    dataset_with_report = trusted_analysis_ready_dataset_from_tables(
        phospho=base_dataset.phospho,
        site_metadata=base_dataset.site_metadata,
        sample_metadata=base_dataset.sample_metadata,
        total=base_dataset.total,
        comparisons=base_dataset.comparisons,
        organism=base_dataset.organism,
        intensity_scale_state=base_dataset.intensity_scale_state,
        processing_state=base_dataset.processing_state,
        preprocessing_report=preprocessing_report,
        provenance=base_dataset.provenance,
    )
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset_with_report,
            design=_request().design,
            contrasts=_request().contrasts,
        )
    )
    assert result.input_dataset_preprocessing_report is preprocessing_report


def test_differential_public_stages_expose_run_only() -> None:
    assert _public_methods(DifferentialAnalysisWorkflow) == {"run"}
    assert _public_methods(DifferentialAnalysisValidator) == {"run"}
    assert _public_methods(DifferentialAnalysisInterpreter) == {"run"}
    assert _public_methods(DifferentialAnalysisExecutor) == {"run"}
    assert not hasattr(DifferentialAnalysisWorkflow, "validate")
    assert not hasattr(DifferentialAnalysisWorkflow, "interpret")
    assert not hasattr(DifferentialAnalysisWorkflow, "execute")


def test_differential_validator_rejects_invalid_raw_request_type() -> None:
    validator = DifferentialAnalysisValidator()
    with pytest.raises(
        WorkflowValidationError,
        match="differential workflow input must be a DifferentialAnalysisRequest",
    ):
        validator.run(object())


def test_differential_validator_rejects_unknown_contrast_term_before_interpretation() -> (
    None
):
    request = _request()
    bad_contrasts = (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A_bad",
        ),
    )
    bad_request = DifferentialAnalysisRequest(
        dataset=request.dataset,
        design=request.design,
        contrasts=bad_contrasts,
    )
    with pytest.raises(
        WorkflowValidationError,
        match="unknown denominator condition",
    ):
        DifferentialAnalysisValidator().run(bad_request)


def test_differential_validator_does_not_run_statistical_execution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_if_called(*args: object, **kwargs: object) -> object:
        raise AssertionError("statistical execution must not run during validation")

    monkeypatch.setattr(DifferentialComputationExecutor, "run", fail_if_called)

    request = _request()
    validated = DifferentialAnalysisValidator().run(request)

    assert isinstance(validated, ValidatedDifferentialAnalysisRequest)
    assert validated.dataset is request.dataset
    assert not hasattr(validated, "computation_request")
    assert not hasattr(validated, "result_identity_metadata")


def test_differential_interpreter_checks_sample_to_design_alignment() -> None:
    request = _request()
    validated = DifferentialAnalysisValidator().run(request)
    misaligned_design = validated.design_matrix.to_dataframe()
    misaligned_design.index = pd.Index(["x1", "x2", "x3", "x4"], name="sample")
    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.interpreter.sample_label_alignment",
    ):
        DifferentialAnalysisInterpreter().run(
            ValidatedDifferentialAnalysisRequest(
                dataset=validated.dataset,
                design=validated.design,
                contrasts=validated.contrasts,
                analysis_sample_ids=validated.analysis_sample_ids,
                design_matrix=type(validated.design_matrix)(misaligned_design),
                contrast_matrix=validated.contrast_matrix,
                design_decomposition=validated.design_decomposition,
                config=validated.config,
            )
        )


def test_differential_executor_accepts_only_interpreted_requests() -> None:
    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.executor.interpreted_request_type",
    ):
        DifferentialAnalysisExecutor().run(object())  # type: ignore[arg-type]


def test_differential_executor_consumes_interpreter_resolved_design_inputs() -> None:
    interpreted = DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(_request())
    )
    assert interpreted.execution_design is not None
    assert isinstance(
        interpreted.execution_config,
        ResolvedDifferentialExecutionConfig,
    )
    assert interpreted.execution_config.minimum_condition_replicates == (
        interpreted.config.minimum_condition_replicates
    )
    assert interpreted.computation_request.empirical_bayes is (
        interpreted.execution_config.empirical_bayes
    )
    assert interpreted.computation_request.multiple_testing_method == (
        interpreted.execution_config.multiple_testing_method
    )

    class _ComputationExecutorSpy:
        def __init__(self) -> None:
            self.received_request = None
            self._real_executor = DifferentialComputationExecutor()

        def run(self, request):
            self.received_request = request
            assert request is interpreted.computation_request
            return self._real_executor.run(request)

    computation_executor = _ComputationExecutorSpy()
    result = DifferentialAnalysisExecutor(
        computation_executor=computation_executor,  # type: ignore[arg-type]
    ).run(interpreted)

    assert computation_executor.received_request is interpreted.computation_request
    pd.testing.assert_frame_equal(
        interpreted.computation_request.design.to_dataframe(),
        interpreted.execution_design.design_matrix.to_dataframe(),
    )
    pd.testing.assert_frame_equal(
        interpreted.computation_request.contrasts.to_dataframe(),
        interpreted.execution_design.contrast_matrix.to_dataframe(),
    )
    assert "B_vs_A" in result.contrast_tables


def test_differential_policy_provenance_rejects_different_decomposition_object() -> (
    None
):
    validated = DifferentialAnalysisValidator().run(_request())
    rebuilt_decomposition = decompose_differential_design(
        validated.design_matrix.frame.to_numpy(dtype=float)
    )

    assert rebuilt_decomposition is not validated.design_decomposition
    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.provenance.design_decomposition_identity",
    ):
        build_differential_policy_provenance(
            request=validated,
            design_decomposition=rebuilt_decomposition,
        )


def test_differential_workflow_executor_does_not_mutate_computation_result_tables() -> (
    None
):
    interpreted = DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(_request())
    )

    class _ComputationResultWithoutPrivateContrastTables:
        def __init__(self, result: Any) -> None:
            self._result = result

        def __getattr__(self, name: str) -> object:
            if name == "_contrast_tables":
                raise AssertionError(
                    "workflow executor must use the differential internal view"
                )
            return getattr(self._result, name)

        def _borrow_contrast_tables(self) -> Mapping[str, pd.DataFrame]:
            return DifferentialComputationResultInternalView(
                self._result
            ).contrast_tables

    class _ComputationExecutorSpy:
        def __init__(self) -> None:
            self.result = None
            self.before_tables = None
            self._real_executor = DifferentialComputationExecutor()

        def run(self, request):
            self.result = self._real_executor.run(request)
            self.before_tables = {
                contrast_name: table.copy(deep=True)
                for contrast_name, table in DifferentialComputationResultInternalView(
                    self.result
                ).contrast_tables.items()
            }
            return _ComputationResultWithoutPrivateContrastTables(self.result)

    computation_executor = _ComputationExecutorSpy()
    workflow_result = DifferentialAnalysisExecutor(
        computation_executor=computation_executor,  # type: ignore[arg-type]
    ).run(interpreted)

    assert computation_executor.result is not None
    assert computation_executor.before_tables is not None
    for contrast_name, before_table in computation_executor.before_tables.items():
        pd.testing.assert_frame_equal(
            computation_executor.result.table_for(contrast_name),
            before_table,
        )
    assert "site_key" in workflow_result.table_for("B_vs_A").columns


def test_differential_invalid_contrast_fails_before_executor() -> None:
    calls = {"executor": 0}

    class _ExecutorSpy:
        def run(self, request: InterpretedDifferentialAnalysisRequest):
            calls["executor"] += 1
            raise AssertionError("executor should not be called")

    request = _request()
    bad_request = DifferentialAnalysisRequest(
        dataset=request.dataset,
        design=request.design,
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A_bad",
            ),
        ),
    )

    workflow = DifferentialAnalysisWorkflow._with_components(executor=_ExecutorSpy())  # type: ignore[arg-type]
    with pytest.raises(WorkflowValidationError):
        workflow.run(bad_request)
    assert calls["executor"] == 0


def test_differential_misaligned_design_fails_before_executor() -> None:
    calls = {"executor": 0}

    class _ExecutorSpy:
        def run(self, request: InterpretedDifferentialAnalysisRequest):
            calls["executor"] += 1
            raise AssertionError("executor should not be called")

    request = _request()
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="x1", condition="A"),
            SampleDesignRecord(sample_id="x2", condition="A"),
            SampleDesignRecord(sample_id="x3", condition="B"),
            SampleDesignRecord(sample_id="x4", condition="B"),
        )
    )

    workflow = DifferentialAnalysisWorkflow._with_components(executor=_ExecutorSpy())  # type: ignore[arg-type]
    with pytest.raises(WorkflowValidationError):
        workflow.run(
            DifferentialAnalysisRequest(
                dataset=request.dataset,
                design=design,
                contrasts=request.contrasts,
            )
        )
    assert calls["executor"] == 0


def test_differential_validator_rejects_non_differential_config_type() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="differential workflow request config must be DifferentialAnalysisConfig",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_request().design,
                contrasts=_request().contrasts,
                config=object(),  # type: ignore[arg-type]
            )
        )


def test_differential_validator_rejects_established_linear_scale() -> None:
    valid_dataset = _dataset()
    linear_dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=valid_dataset.phospho,
        site_metadata=valid_dataset.site_metadata,
        organism=valid_dataset.organism,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="requires established log2-scale phospho intensities",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=linear_dataset,
                design=_request().design,
                contrasts=_request().contrasts,
            )
        )


def test_differential_validator_rejects_declared_but_unestablished_log2_scale() -> None:
    dataset = _dataset()
    unsafe_replace_dataset_intensity_scale_state(
        dataset,
        IntensityScaleState(
            phospho=MatrixIntensityScaleState.log2(established_by="test.declaration"),
            total=None,
        ),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="requires established log2-scale phospho intensities",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=dataset,
                design=_request().design,
                contrasts=_request().contrasts,
            )
        )


def test_differential_invalid_scale_fails_before_executor() -> None:
    calls = {"executor": 0}

    class _ExecutorSpy:
        def run(self, request: InterpretedDifferentialAnalysisRequest):
            calls["executor"] += 1
            raise AssertionError("executor should not be called")

    valid_dataset = _dataset()
    linear_dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=valid_dataset.phospho,
        site_metadata=valid_dataset.site_metadata,
        organism=valid_dataset.organism,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    with pytest.raises(WorkflowValidationError):
        DifferentialAnalysisWorkflow._with_components(executor=_ExecutorSpy()).run(  # type: ignore[arg-type]
            DifferentialAnalysisRequest(
                dataset=linear_dataset,
                design=_request().design,
                contrasts=_request().contrasts,
            )
        )
    assert calls["executor"] == 0


def test_differential_unknown_scale_fails_before_executor() -> None:
    calls = {"executor": 0}

    class _ExecutorSpy:
        def run(self, request: InterpretedDifferentialAnalysisRequest):
            calls["executor"] += 1
            raise AssertionError("executor should not be called")

    dataset = _dataset()
    unknown_state = IntensityScaleState.raw(has_total_matrix=False)
    unsafe_replace_dataset_intensity_scale_state(dataset, unknown_state)

    with pytest.raises(
        WorkflowValidationError,
        match="requires established log2-scale phospho intensities",
    ):
        DifferentialAnalysisWorkflow._with_components(executor=_ExecutorSpy()).run(  # type: ignore[arg-type]
            DifferentialAnalysisRequest(
                dataset=dataset,
                design=_request().design,
                contrasts=_request().contrasts,
            )
        )
    assert calls["executor"] == 0


def test_differential_workflow_rejects_imputed_dataset_before_executor() -> None:
    calls = {"executor": 0}

    class _ExecutorSpy:
        def run(self, request: InterpretedDifferentialAnalysisRequest):
            calls["executor"] += 1
            raise AssertionError("executor should not be called")

    request = _request()
    processing_state = request.dataset.processing_state
    imputed_dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=request.dataset.phospho,
        site_metadata=request.dataset.site_metadata,
        organism=request.dataset.organism,
        intensity_scale_state=request.dataset.intensity_scale_state,
        processing_state=valid_imputed_processing_state(processing_state),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="imputed cells as observed measurements",
    ):
        DifferentialAnalysisWorkflow._with_components(executor=_ExecutorSpy()).run(  # type: ignore[arg-type]
            DifferentialAnalysisRequest(
                dataset=imputed_dataset,
                design=request.design,
                contrasts=request.contrasts,
            )
        )
    assert calls["executor"] == 0


def test_differential_validator_rejects_raw_string_technical_replicate_policy() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="technical_replicate_policy must be TechnicalReplicatePolicy",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_request().design,
                contrasts=_request().contrasts,
                config=DifferentialAnalysisConfig(
                    technical_replicate_policy="mean",  # type: ignore[arg-type]
                ),
            )
        )


def test_differential_validator_rejects_non_bool_allow_design_subset() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="allow_design_subset must be a bool",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_request().design,
                contrasts=_request().contrasts,
                config=DifferentialAnalysisConfig(
                    allow_design_subset=1,  # type: ignore[arg-type]
                ),
            )
        )


def test_differential_validator_rejects_non_bool_declared_scale_override() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="allow_suspicious_declared_input_scale must be a bool",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_request().design,
                contrasts=_request().contrasts,
                config=DifferentialAnalysisConfig(
                    allow_suspicious_declared_input_scale=1,  # type: ignore[arg-type]
                ),
            )
        )


@pytest.mark.parametrize("value", [0, -1, 1.5, True])
def test_differential_validator_rejects_invalid_minimum_condition_replicates(
    value: object,
) -> None:
    pattern = (
        "minimum_condition_replicates must be >= 1"
        if isinstance(value, int) and not isinstance(value, bool)
        else "minimum_condition_replicates must be an int"
    )
    with pytest.raises(WorkflowValidationError, match=pattern):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_request().design,
                contrasts=_request().contrasts,
                config=DifferentialAnalysisConfig(
                    minimum_condition_replicates=value,  # type: ignore[arg-type]
                ),
            )
        )


def test_differential_validator_rejects_non_empirical_bayes_config() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="empirical_bayes must be EmpiricalBayesConfig",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_request().design,
                contrasts=_request().contrasts,
                config=DifferentialAnalysisConfig(
                    empirical_bayes=object(),  # type: ignore[arg-type]
                ),
            )
        )


def test_differential_validator_rejects_non_multiple_testing_config() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="multiple_testing must be MultipleTestingConfig",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_request().design,
                contrasts=_request().contrasts,
                config=DifferentialAnalysisConfig(
                    multiple_testing=object(),  # type: ignore[arg-type]
                ),
            )
        )


def test_differential_validator_passes_config_values_to_collaborators() -> None:
    calls: dict[str, object] = {}
    call_order: list[str] = []
    collaborator_views: list[DatasetInternalView] = []
    base_request = _request()

    class _FakeDatasetEligibilityValidator:
        def run(
            self,
            *,
            dataset: AnalysisReadyPhosphoDataset,
            imputed_value_policy: DifferentialImputedValuePolicy,
            allow_suspicious_declared_input_scale: bool,
            dataset_view: DatasetInternalView,
        ) -> None:
            call_order.append("dataset_eligibility_validator")
            collaborator_views.append(dataset_view)
            calls["imputed_value_policy"] = imputed_value_policy
            calls["allow_suspicious_declared_input_scale"] = (
                allow_suspicious_declared_input_scale
            )

    class _FakeTechnicalReplicatePlanner:
        def run(
            self,
            *,
            dataset: AnalysisReadyPhosphoDataset,
            design: ExperimentalDesign,
            technical_replicate_policy: TechnicalReplicatePolicy,
            dataset_view: DatasetInternalView,
        ) -> TechnicalReplicateAggregationPlan:
            call_order.append("technical_replicate_planner")
            collaborator_views.append(dataset_view)
            calls["technical_replicate_policy"] = technical_replicate_policy

            return TechnicalReplicateAggregationPlan(
                technical_replicate_policy=technical_replicate_policy,
                groups=(),
                aggregate_phospho=False,
                aggregate_total_protein=False,
            )

    class _FakeDesignValidator:
        def run(
            self,
            *,
            dataset: AnalysisReadyPhosphoDataset,
            design: ExperimentalDesign,
            contrasts: tuple[Contrast, ...],
            allow_design_subset: bool,
            minimum_condition_replicates: int,
            paired_design_policy: PairedDesignPolicy,
            dataset_view: DatasetInternalView,
        ) -> ValidatedExperimentalDesignContract:
            call_order.append("design_validator")
            collaborator_views.append(dataset_view)
            calls["allow_design_subset"] = allow_design_subset
            calls["minimum_condition_replicates"] = minimum_condition_replicates
            calls["paired_design_policy"] = paired_design_policy

            design_frame = pd.DataFrame(
                {
                    "A": [1.0, 1.0, 0.0, 0.0],
                    "B": [0.0, 0.0, 1.0, 1.0],
                },
                index=pd.Index(["A_1", "A_2", "B_1", "B_2"], name="sample"),
            )
            return ValidatedExperimentalDesignContract(
                design=design,
                contrasts=contrasts,
                analysis_sample_ids=tuple(
                    record.sample_id for record in design.samples
                ),
                condition_labels=("A", "B"),
                design_frame=design_frame,
                contrast_frame=pd.DataFrame(
                    {"B_vs_A": [-1.0, 1.0]},
                    index=pd.Index(["A", "B"], name="coefficient"),
                ),
                design_decomposition=decompose_differential_design(
                    design_frame.to_numpy(dtype=float)
                ),
            )

    request = DifferentialAnalysisRequest(
        dataset=base_request.dataset,
        design=base_request.design,
        contrasts=base_request.contrasts,
        config=DifferentialAnalysisConfig(
            reliability_profile=(
                DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE
            ),
            technical_replicate_policy=TechnicalReplicatePolicy.MEAN,
            allow_design_subset=True,
            minimum_condition_replicates=1,
            empirical_bayes=EmpiricalBayesConfig(method="standard"),
            multiple_testing=MultipleTestingConfig(),
        ),
    )
    validated = DifferentialAnalysisValidator(
        dataset_eligibility_validator=_FakeDatasetEligibilityValidator(),
        technical_replicate_planner=_FakeTechnicalReplicatePlanner(),
        design_validator=_FakeDesignValidator(),
    ).run(request)

    assert call_order == [
        "dataset_eligibility_validator",
        "technical_replicate_planner",
        "design_validator",
    ]
    assert all(view is validated.dataset_view for view in collaborator_views)
    assert calls["imputed_value_policy"] == "reject"
    assert calls["allow_suspicious_declared_input_scale"] is False
    assert calls["technical_replicate_policy"] is TechnicalReplicatePolicy.MEAN
    assert calls["allow_design_subset"] is True
    assert calls["minimum_condition_replicates"] == 1
    assert calls["paired_design_policy"] == "reject"

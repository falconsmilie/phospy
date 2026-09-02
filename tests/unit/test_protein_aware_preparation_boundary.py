from __future__ import annotations

import ast
from dataclasses import replace
from pathlib import Path
from typing import Any, cast

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy.advanced import (
    DatasetProteinAwarePreparationConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    Organism,
)
from phospy.errors import DatasetValidationError, PhosPyInputError
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import DatasetPreprocessingReport
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.datasets.preprocessing.policy_models import (
    TotalProteinCorrectionPolicy,
)
from phospy.science.datasets.preprocessing.protein_aware_alignment import (
    ProteinAwareSampleAlignmentDiagnostics,
)
from phospy.science.datasets.preprocessing.protein_aware_preparation import (
    ProteinAwarePreparationReport,
    ProteinAwarePreparationResult,
    ProteinAwareSiteEligibility,
)
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)

ROOT = Path(__file__).resolve().parents[2]
_USE_DATASET_VALUE = object()
_DIFFERENTIAL_CODE_DIRS = (
    ROOT / "src" / "phospy" / "workflows" / "differential",
    ROOT / "src" / "phospy" / "science" / "differential",
)
_FORBIDDEN_DIFFERENTIAL_OWNER_MODULE_PREFIXES = (
    "phospy.science.datasets.construction",
    "phospy.science.datasets.preprocessing.pipeline",
    "phospy.science.datasets.preprocessing.plan",
    "phospy.science.datasets.preprocessing.plan_",
    "phospy.science.datasets.preprocessing.stage_",
    "phospy.science.datasets.preprocessing.stages.protein_aware_preparation",
    "phospy.science.datasets.preprocessing.state_builder",
    "phospy.science.datasets.preprocessing.trace_builder",
    "phospy.science.datasets.preprocessing.protein_aware_preparation",
    "phospy.science.datasets.preprocessing.protein_mapping",
)
_FORBIDDEN_DIFFERENTIAL_OWNER_AGGREGATE_MODULES = (
    "phospy.science.datasets.preprocessing",
)
_FORBIDDEN_DIFFERENTIAL_OWNER_NAMES = (
    "AnalysisReadyDatasetBuilder",
    "DATASET_PREPROCESSING_STAGE_PROTEIN_AWARE_PREPARATION",
    "DatasetBuildRequest",
    "DatasetProteinAwarePreparationConfig",
    "PreprocessingPlan",
    "PreprocessingPipeline",
    "PreprocessingStage",
    "PreprocessingStageResult",
    "PreprocessingState",
    "ProteinAwareAlignmentConfig",
    "ProteinAwareAlignmentEligibilityDiagnostics",
    "ProteinAwareAlignmentEligibilityResolver",
    "ProteinAwareSiteEligibilityDiagnostic",
    "ProteinAwarePreparationStage",
    "ProteinMappingConfig",
    "ProteinMappingResolver",
    "ProteinMappingRecord",
    "ProteinMappingResult",
)
_FORBIDDEN_DIFFERENTIAL_PREPARATION_CONSTRUCTORS = (
    "ProteinAwarePreparationReport",
    "ProteinAwarePreparationResult",
)
_FORBIDDEN_DIFFERENTIAL_PREPARATION_PRIVATE_FIELDS = (
    "_matched_pairs",
    "_protein_covariate_matrix",
)
_FORBIDDEN_DIFFERENTIAL_PREPARATION_SETATTR_FIELDS = (
    *_FORBIDDEN_DIFFERENTIAL_PREPARATION_PRIVATE_FIELDS,
    "report",
)


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "A_1": [1.0, 1.1],
            "A_2": [1.2, 1.3],
            "B_1": [2.0, 2.1],
            "B_2": [2.2, 2.3],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )


def _site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "protein_id": ["P53778", "P31749"],
            "site_sequence": [
                "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            ],
            "localisation_confidence": [0.95, 0.96],
        },
        index=_phospho().index.copy(),
    )


def _total(*, include_akt1: bool = True) -> pd.DataFrame:
    index = ["P53778", "P31749"] if include_akt1 else ["P53778"]
    rows = {
        "P53778": [10.0, 11.0, 12.0, 13.0],
        "P31749": [20.0, 21.0, 22.0, 23.0],
    }
    return pd.DataFrame(
        {
            sample_id: [rows[row_key][position] for row_key in index]
            for position, sample_id in enumerate(("A_1", "A_2", "B_1", "B_2"))
        },
        index=pd.Index(index, name="protein_id"),
    )


def _build_dataset(
    *,
    preprocessing_config: DatasetPreprocessingConfig,
    total: pd.DataFrame | None = None,
):
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            total=_total() if total is None else total,
            organism=Organism.HUMAN,
            input_intensity_scale="log2",
            preprocessing_config=preprocessing_config,
        )
    )


def _prepared_dataset() -> AnalysisReadyPhosphoDataset:
    return _build_dataset(
        preprocessing_config=DatasetPreprocessingConfig(
            protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                policy="prepare_model_inputs",
                protein_mapping_policy="require_unambiguous",
            )
        )
    )


def _replay_dataset(
    dataset: AnalysisReadyPhosphoDataset,
    *,
    phospho: pd.DataFrame | None = None,
    site_metadata: pd.DataFrame | None = None,
    total: object = _USE_DATASET_VALUE,
    preprocessing_report: object = _USE_DATASET_VALUE,
    protein_aware_preparation: object = _USE_DATASET_VALUE,
    intensity_scale_state: object | None = None,
    processing_state: object | None = None,
) -> AnalysisReadyPhosphoDataset:
    total_value = (
        dataset.total
        if total is _USE_DATASET_VALUE
        else cast(pd.DataFrame | None, total)
    )
    report_value = (
        dataset.preprocessing_report
        if preprocessing_report is _USE_DATASET_VALUE
        else cast(DatasetPreprocessingReport | None, preprocessing_report)
    )
    preparation_value = (
        dataset.protein_aware_preparation
        if protein_aware_preparation is _USE_DATASET_VALUE
        else cast(ProteinAwarePreparationResult | None, protein_aware_preparation)
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=dataset.phospho if phospho is None else phospho,
        site_metadata=dataset.site_metadata if site_metadata is None else site_metadata,
        sample_metadata=dataset.sample_metadata,
        total=total_value,
        organism=dataset.organism,
        intensity_scale_state=(
            dataset.intensity_scale_state
            if intensity_scale_state is None
            else cast(Any, intensity_scale_state)
        ),
        processing_state=(
            dataset.processing_state
            if processing_state is None
            else cast(Any, processing_state)
        ),
        preprocessing_report=report_value,
        protein_aware_preparation=preparation_value,
    )


def _preparation_report_with(
    report: ProteinAwarePreparationReport,
    *,
    site_eligibility: tuple[ProteinAwareSiteEligibility, ...] | None = None,
    sample_alignment: ProteinAwareSampleAlignmentDiagnostics | None = None,
    policy_parameters: dict[str, object] | None = None,
    preparation_policy: str | None = None,
    schema_version: int | None = None,
) -> ProteinAwarePreparationReport:
    return ProteinAwarePreparationReport(
        site_eligibility=(
            report.site_eligibility if site_eligibility is None else site_eligibility
        ),
        mapping_diagnostics=report.mapping_diagnostics,
        sample_alignment=(
            report.sample_alignment if sample_alignment is None else sample_alignment
        ),
        transformation_state=report.transformation_state,
        preparation_policy=(
            report.preparation_policy
            if preparation_policy is None
            else preparation_policy
        ),
        protein_mapping_policy=report.protein_mapping_policy,
        policy_parameters=(
            dict(report.policy_parameters)
            if policy_parameters is None
            else policy_parameters
        ),
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


def _preprocessing_report_for(
    preparation: ProteinAwarePreparationResult,
) -> DatasetPreprocessingReport:
    return DatasetPreprocessingReport.from_rows(
        protein_aware_preparation=preparation.report
    )


def test_protein_aware_preparation_output_is_separate_from_analysis_ready_dataset() -> (
    None
):
    disabled = _build_dataset(preprocessing_config=DatasetPreprocessingConfig())
    prepared = _build_dataset(
        preprocessing_config=DatasetPreprocessingConfig(
            protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                policy="prepare_model_inputs",
                protein_mapping_policy="require_unambiguous",
            )
        )
    )

    assert prepared.protein_aware_preparation is not None
    assert prepared.preprocessing_report is not None
    assert (
        prepared.preprocessing_report.protein_aware_preparation
        is prepared.protein_aware_preparation.report
    )
    pdt.assert_frame_equal(prepared.phospho, disabled.phospho)
    pdt.assert_frame_equal(prepared.site_metadata, disabled.site_metadata)
    assert (
        prepared.protein_aware_preparation.matched_pairs_dataframe()
        .loc[:, "site_key"]
        .tolist()
        == prepared.phospho.index.astype(str).tolist()
    )
    covariate_index = prepared.protein_aware_preparation.protein_covariate_matrix_dataframe().index.tolist()
    assert covariate_index == [
        "P53778",
        "P31749",
    ]


def test_valid_builder_owned_sidecar_can_be_replayed_through_trusted_construction() -> (
    None
):
    prepared = _prepared_dataset()

    replayed = _replay_dataset(prepared)

    assert replayed.protein_aware_preparation is not None
    assert prepared.protein_aware_preparation is not None
    assert replayed.protein_aware_preparation.scientifically_equals(
        prepared.protein_aware_preparation
    )


def test_dataset_preprocessing_rejects_subtraction_and_model_input_preparation() -> (
    None
):
    with pytest.raises(
        PhosPyInputError,
        match="ambiguous double protein adjustment",
    ):
        DatasetPreprocessingConfig(
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy="subtract_log_total"
            ),
            protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                policy="prepare_model_inputs"
            ),
        )

    with pytest.raises(
        PhosPyInputError,
        match="ambiguous double protein adjustment",
    ):
        PreprocessingPlan(
            total_protein_correction_policy=(
                TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
            ),
            protein_aware_preparation_policy="prepare_model_inputs",
        )


def test_trusted_construction_rejects_sidecar_from_other_dataset_same_dimensions() -> (
    None
):
    prepared = _prepared_dataset()
    total = prepared.total
    assert total is not None
    changed_total = total.copy(deep=True)
    changed_total.iloc[0, 0] = float(changed_total.iloc[0, 0]) + 100.0

    with pytest.raises(DatasetValidationError, match="protein_covariate_matrix"):
        _replay_dataset(prepared, total=changed_total)


def test_trusted_construction_rejects_copied_sidecar_from_other_dataset() -> None:
    prepared = _prepared_dataset()
    other_total = _total()
    other_total.iloc[0, 0] = float(other_total.iloc[0, 0]) + 0.5
    other = _build_dataset(
        preprocessing_config=DatasetPreprocessingConfig(
            protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                policy="prepare_model_inputs",
                protein_mapping_policy="require_unambiguous",
            )
        ),
        total=other_total,
    )

    with pytest.raises(DatasetValidationError, match="protein_covariate_matrix"):
        _replay_dataset(
            prepared,
            preprocessing_report=other.preprocessing_report,
            protein_aware_preparation=other.protein_aware_preparation,
        )


def test_trusted_construction_rejects_changed_phospho_values_for_bound_sidecar() -> (
    None
):
    prepared = _prepared_dataset()
    changed_phospho = prepared.phospho
    changed_phospho.iloc[0, 0] = float(changed_phospho.iloc[0, 0]) + 1.0

    with pytest.raises(
        DatasetValidationError,
        match="does not match actual analysis-ready dataset tables",
    ):
        _replay_dataset(prepared, phospho=changed_phospho)


def test_trusted_construction_rejects_unsupported_sidecar_schema_version() -> None:
    prepared = _prepared_dataset()
    preparation = prepared.protein_aware_preparation
    assert preparation is not None
    report = _preparation_report_with(preparation.report, schema_version=999)
    altered_preparation = _preparation_with(preparation, report=report)

    with pytest.raises(DatasetValidationError, match="schema_version"):
        _replay_dataset(
            prepared,
            preprocessing_report=_preprocessing_report_for(altered_preparation),
            protein_aware_preparation=altered_preparation,
        )


def test_trusted_construction_rejects_unsupported_sidecar_preparation_policy() -> None:
    prepared = _prepared_dataset()
    preparation = prepared.protein_aware_preparation
    assert preparation is not None
    report = _preparation_report_with(
        preparation.report,
        preparation_policy="disabled",
    )
    altered_preparation = _preparation_with(preparation, report=report)

    with pytest.raises(DatasetValidationError, match="preparation_policy"):
        _replay_dataset(
            prepared,
            preprocessing_report=_preprocessing_report_for(altered_preparation),
            protein_aware_preparation=altered_preparation,
        )


def test_trusted_construction_rejects_sidecar_without_binding_fingerprints() -> None:
    prepared = _prepared_dataset()
    preparation = prepared.protein_aware_preparation
    assert preparation is not None
    policy_parameters = {
        key: value
        for key, value in dict(preparation.report.policy_parameters).items()
        if key != "dataset_binding_table_fingerprints"
    }
    report = _preparation_report_with(
        preparation.report,
        policy_parameters=policy_parameters,
    )
    altered_preparation = _preparation_with(preparation, report=report)

    with pytest.raises(
        DatasetValidationError,
        match="bound to the owning analysis-ready dataset tables",
    ):
        _replay_dataset(
            prepared,
            preprocessing_report=_preprocessing_report_for(altered_preparation),
            protein_aware_preparation=altered_preparation,
        )


def test_trusted_construction_rejects_altered_phospho_site_key_for_bound_sidecar() -> (
    None
):
    prepared = _prepared_dataset()
    changed_phospho = prepared.phospho
    changed_site_metadata = prepared.site_metadata
    changed_index = changed_phospho.index.astype(str).tolist()
    changed_index[0] = (
        "phospy:v1|organism=human|protein_namespace=protein_id|"
        "protein_identifier=P53778|residue=Y|position=183"
    )
    changed_phospho.index = pd.Index(changed_index, name=changed_phospho.index.name)
    changed_site_metadata.index = pd.Index(
        changed_index,
        name=changed_site_metadata.index.name,
    )
    changed_site_metadata.iloc[
        0,
        changed_site_metadata.columns.get_loc("site"),
    ] = "Y183"
    changed_site_metadata.iloc[
        0,
        changed_site_metadata.columns.get_loc("display_id"),
    ] = "MAPK14;Y183;"
    changed_site_metadata.iloc[
        0,
        changed_site_metadata.columns.get_loc("site_key"),
    ] = changed_index[0]

    with pytest.raises(DatasetValidationError, match="site_eligibility"):
        _replay_dataset(
            prepared,
            phospho=changed_phospho,
            site_metadata=changed_site_metadata,
        )


def test_trusted_construction_rejects_altered_sidecar_site_key() -> None:
    prepared = _prepared_dataset()
    preparation = prepared.protein_aware_preparation
    assert preparation is not None
    altered_site_key = "not-current-site-key"
    matched_pairs = preparation.matched_pairs_dataframe()
    matched_pairs.loc[0, "site_key"] = altered_site_key
    report = _preparation_report_with(
        preparation.report,
        site_eligibility=_replace_eligibility_row(
            preparation.report,
            site_key=altered_site_key,
        ),
    )
    altered_preparation = _preparation_with(
        preparation,
        matched_pairs=matched_pairs,
        report=report,
    )

    with pytest.raises(DatasetValidationError, match="site_eligibility"):
        _replay_dataset(
            prepared,
            preprocessing_report=_preprocessing_report_for(altered_preparation),
            protein_aware_preparation=altered_preparation,
        )


def test_trusted_construction_rejects_report_protein_identifier_change() -> None:
    prepared = _prepared_dataset()
    preparation = prepared.protein_aware_preparation
    assert preparation is not None
    report = _preparation_report_with(
        preparation.report,
        site_eligibility=_replace_eligibility_row(
            preparation.report,
            protein_identifier="Q99999",
        ),
    )
    altered_preparation = _preparation_with(preparation, report=report)

    with pytest.raises(DatasetValidationError, match="must agree"):
        _replay_dataset(
            prepared,
            preprocessing_report=_preprocessing_report_for(altered_preparation),
            protein_aware_preparation=altered_preparation,
        )


def test_trusted_construction_rejects_altered_site_metadata_protein_identifier() -> (
    None
):
    prepared = _prepared_dataset()
    changed_site_metadata = prepared.site_metadata
    changed_site_metadata.loc[
        prepared.phospho.index[0],
        "protein_identifier",
    ] = "Q99999"

    with pytest.raises(DatasetValidationError, match="protein_identifier"):
        _replay_dataset(prepared, site_metadata=changed_site_metadata)


def test_trusted_construction_rejects_altered_sidecar_protein_identifier() -> None:
    prepared = _prepared_dataset()
    preparation = prepared.protein_aware_preparation
    assert preparation is not None
    matched_pairs = preparation.matched_pairs_dataframe()
    matched_pairs.loc[0, "protein_identifier"] = "Q99999"
    report = _preparation_report_with(
        preparation.report,
        site_eligibility=_replace_eligibility_row(
            preparation.report,
            protein_identifier="Q99999",
        ),
    )
    altered_preparation = _preparation_with(
        preparation,
        matched_pairs=matched_pairs,
        report=report,
    )

    with pytest.raises(DatasetValidationError, match="protein_identifier"):
        _replay_dataset(
            prepared,
            preprocessing_report=_preprocessing_report_for(altered_preparation),
            protein_aware_preparation=altered_preparation,
        )


def test_trusted_construction_rejects_altered_sidecar_total_row_key() -> None:
    prepared = _prepared_dataset()
    preparation = prepared.protein_aware_preparation
    assert preparation is not None
    altered_total_row_key = "not-current-total-row"
    matched_pairs = preparation.matched_pairs_dataframe()
    matched_pairs.loc[0, "total_protein_row_key"] = altered_total_row_key
    protein_covariates = preparation.protein_covariate_matrix_dataframe()
    protein_covariates.index = pd.Index(
        [altered_total_row_key, *protein_covariates.index.tolist()[1:]],
        name=protein_covariates.index.name,
    )
    report = _preparation_report_with(
        preparation.report,
        site_eligibility=_replace_eligibility_row(
            preparation.report,
            total_protein_row_key=altered_total_row_key,
        ),
    )
    altered_preparation = _preparation_with(
        preparation,
        matched_pairs=matched_pairs,
        protein_covariate_matrix=protein_covariates,
        report=report,
    )

    with pytest.raises(DatasetValidationError, match="total_protein_row_key"):
        _replay_dataset(
            prepared,
            preprocessing_report=_preprocessing_report_for(altered_preparation),
            protein_aware_preparation=altered_preparation,
        )


def test_trusted_construction_rejects_report_total_row_key_change() -> None:
    prepared = _prepared_dataset()
    preparation = prepared.protein_aware_preparation
    assert preparation is not None
    report = _preparation_report_with(
        preparation.report,
        site_eligibility=_replace_eligibility_row(
            preparation.report,
            total_protein_row_key="P31749",
        ),
    )
    altered_preparation = _preparation_with(preparation, report=report)

    with pytest.raises(DatasetValidationError, match="must agree"):
        _replay_dataset(
            prepared,
            preprocessing_report=_preprocessing_report_for(altered_preparation),
            protein_aware_preparation=altered_preparation,
        )


def test_trusted_construction_rejects_altered_sidecar_total_protein_value() -> None:
    prepared = _prepared_dataset()
    preparation = prepared.protein_aware_preparation
    assert preparation is not None
    protein_covariates = preparation.protein_covariate_matrix_dataframe()
    protein_covariates.iloc[0, 0] = float(protein_covariates.iloc[0, 0]) + 1.0
    altered_preparation = _preparation_with(
        preparation,
        protein_covariate_matrix=protein_covariates,
    )

    with pytest.raises(DatasetValidationError, match="protein_covariate_matrix"):
        _replay_dataset(
            prepared,
            preprocessing_report=_preprocessing_report_for(altered_preparation),
            protein_aware_preparation=altered_preparation,
        )


def test_trusted_construction_rejects_reordered_sidecar_samples() -> None:
    prepared = _prepared_dataset()
    preparation = prepared.protein_aware_preparation
    assert preparation is not None
    reversed_columns = tuple(
        reversed(preparation.protein_covariate_matrix_dataframe().columns.tolist())
    )
    protein_covariates = preparation.protein_covariate_matrix_dataframe().loc[
        :,
        list(reversed_columns),
    ]
    report = _preparation_report_with(
        preparation.report,
        sample_alignment=ProteinAwareSampleAlignmentDiagnostics(
            phospho_sample_columns=reversed_columns,
            total_protein_sample_columns=reversed_columns,
            exact_sample_order_match=True,
            sample_order_compatible=True,
            reordered_sample_columns=False,
            allow_reordered_samples=False,
            missing_total_protein_samples=(),
            extra_total_protein_samples=(),
        ),
    )
    altered_preparation = _preparation_with(
        preparation,
        protein_covariate_matrix=protein_covariates,
        report=report,
    )

    with pytest.raises(DatasetValidationError, match="sample_alignment"):
        _replay_dataset(
            prepared,
            preprocessing_report=_preprocessing_report_for(altered_preparation),
            protein_aware_preparation=altered_preparation,
        )


def test_trusted_construction_rejects_sidecar_report_mismatch() -> None:
    prepared = _prepared_dataset()
    preparation = prepared.protein_aware_preparation
    assert preparation is not None
    report_parameters = dict(preparation.report.policy_parameters)
    report_parameters["test_mismatch"] = True
    mismatched_report = _preparation_report_with(
        preparation.report,
        policy_parameters=report_parameters,
    )

    with pytest.raises(DatasetValidationError, match="preprocessing_report"):
        _replay_dataset(
            prepared,
            preprocessing_report=DatasetPreprocessingReport.from_rows(
                protein_aware_preparation=mismatched_report
            ),
            protein_aware_preparation=preparation,
        )


def test_trusted_construction_rejects_sidecar_without_total_table() -> None:
    prepared = _prepared_dataset()
    assert prepared.protein_aware_preparation is not None

    with pytest.raises(DatasetValidationError, match="requires dataset.total"):
        _replay_dataset(
            prepared,
            total=None,
            intensity_scale_state=supported_log2_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_log2_processing_state(has_total_matrix=False),
        )


def test_internal_view_protein_aware_preparation_is_defensive() -> None:
    prepared = _prepared_dataset()
    preparation = DatasetInternalView(prepared).protein_aware_preparation
    assert preparation is not None

    matched_pairs = preparation.matched_pairs
    protein_covariates = preparation.protein_covariate_matrix
    site_eligibility = preparation.site_eligibility
    matched_pairs.loc[0, "protein_identifier"] = "CHANGED"
    protein_covariates.iloc[0, 0] = -999.0
    site_eligibility.loc[0, "eligibility"] = "changed"

    fresh = DatasetInternalView(prepared).protein_aware_preparation
    assert fresh is not None
    assert fresh.matched_pairs.loc[0, "protein_identifier"] == "P53778"
    assert float(fresh.protein_covariate_matrix.iloc[0, 0]) == 10.0
    assert (
        fresh.site_eligibility.loc[0, "eligibility"]
        == "eligible_for_protein_aware_preparation"
    )


def test_internal_view_has_no_protein_aware_sidecar_when_dataset_has_none() -> None:
    dataset = _build_dataset(preprocessing_config=DatasetPreprocessingConfig())

    assert DatasetInternalView(dataset).protein_aware_preparation is None


def test_protein_aware_diagnostics_are_retained_apart_from_phospho_matrix() -> None:
    disabled = _build_dataset(
        preprocessing_config=DatasetPreprocessingConfig(),
        total=_total(include_akt1=False),
    )
    prepared = _build_dataset(
        preprocessing_config=DatasetPreprocessingConfig(
            protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                policy="prepare_model_inputs",
                protein_mapping_policy="allow_missing_with_report",
            )
        ),
        total=_total(include_akt1=False),
    )

    assert prepared.protein_aware_preparation is not None
    diagnostics = (
        prepared.protein_aware_preparation.missing_protein_abundance_diagnostics
    )
    assert diagnostics.to_dict(orient="records") == [
        {
            "site_key": prepared.protein_aware_preparation.report.fallback_site_keys[0],
            "protein_identifier": "P31749",
            "mapping_status": "missing_total_protein_row",
            "reason": "missing_total_protein_row",
        }
    ]
    assert prepared.preprocessing_report is not None
    assert (
        prepared.preprocessing_report.protein_aware_preparation
        is prepared.protein_aware_preparation.report
    )
    covariate_index = prepared.protein_aware_preparation.protein_covariate_matrix_dataframe().index.tolist()
    assert covariate_index == ["P53778"]
    pdt.assert_frame_equal(prepared.phospho, disabled.phospho)
    for diagnostic_column in ("mapping_status", "reason", "protein_covariate_matrix"):
        assert diagnostic_column not in prepared.phospho.columns


def test_differential_domains_only_consume_typed_protein_aware_sidecar_view() -> None:
    violations: list[str] = []
    for path in _differential_python_sources():
        source = path.read_text(encoding="utf-8")
        relative_path = path.relative_to(ROOT).as_posix()
        tree = ast.parse(source, filename=str(path))
        forbidden_call_names = _forbidden_local_call_names(tree)

        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    module = alias.name
                    if _is_forbidden_differential_owner_module(module):
                        violations.append(
                            f"{relative_path}:{node.lineno}: imports {module!r}"
                        )
                    if module in _FORBIDDEN_DIFFERENTIAL_OWNER_AGGREGATE_MODULES:
                        violations.append(
                            f"{relative_path}:{node.lineno}: imports aggregate "
                            f"owner module {module!r}"
                        )
            elif isinstance(node, ast.ImportFrom):
                module = node.module or ""
                if _is_forbidden_differential_owner_module(module):
                    violations.append(
                        f"{relative_path}:{node.lineno}: imports from {module!r}"
                    )
                if module in _FORBIDDEN_DIFFERENTIAL_OWNER_AGGREGATE_MODULES:
                    violations.append(
                        f"{relative_path}:{node.lineno}: imports from aggregate "
                        f"owner module {module!r}"
                    )
                forbidden_names = sorted(
                    alias.name
                    for alias in node.names
                    if alias.name in _FORBIDDEN_DIFFERENTIAL_OWNER_NAMES
                )
                if forbidden_names:
                    violations.append(
                        f"{relative_path}:{node.lineno}: imports owner names "
                        f"{forbidden_names!r}"
                    )
            elif isinstance(node, ast.Call):
                called_name = _qualified_name(node.func)
                called_leaf = called_name.rsplit(".", maxsplit=1)[-1]
                if (
                    called_name in forbidden_call_names
                    or called_leaf in forbidden_call_names
                ):
                    violations.append(
                        f"{relative_path}:{node.lineno}: constructs {called_name}"
                    )
                setattr_field = _forbidden_setattr_field(node)
                if setattr_field is not None:
                    violations.append(
                        f"{relative_path}:{node.lineno}: mutates sidecar field "
                        f"{setattr_field!r}"
                    )
            elif isinstance(node, ast.Assign | ast.AnnAssign | ast.AugAssign):
                for field_name in _assigned_forbidden_private_fields(node):
                    violations.append(
                        f"{relative_path}:{node.lineno}: assigns sidecar private "
                        f"field {field_name!r}"
                    )
            elif isinstance(node, ast.Delete):
                for field_name in _deleted_forbidden_private_fields(node):
                    violations.append(
                        f"{relative_path}:{node.lineno}: deletes sidecar private "
                        f"field {field_name!r}"
                    )

    assert not violations, (
        "protein-aware preparation and mapping implementation must stay in dataset "
        "preprocessing/building domains. Differential may import the dataset "
        "internal view or typed passive sidecar models, but it must not import "
        "preparation stages, mapping resolvers, preprocessing planners, builder "
        "orchestration, or create/mutate preparation sidecars:\n"
        + "\n".join(violations)
    )


def _differential_python_sources() -> tuple[Path, ...]:
    paths: list[Path] = []
    for directory in _DIFFERENTIAL_CODE_DIRS:
        paths.extend(sorted(directory.rglob("*.py")))
    return tuple(paths)


def _is_forbidden_differential_owner_module(module: str) -> bool:
    return any(
        module == prefix
        or module.startswith(f"{prefix}.")
        or (prefix.endswith("_") and module.startswith(prefix))
        for prefix in _FORBIDDEN_DIFFERENTIAL_OWNER_MODULE_PREFIXES
    )


def _forbidden_local_call_names(tree: ast.AST) -> frozenset[str]:
    forbidden_names = set(_FORBIDDEN_DIFFERENTIAL_PREPARATION_CONSTRUCTORS)
    forbidden_names.update(_FORBIDDEN_DIFFERENTIAL_OWNER_NAMES)
    for node in ast.walk(tree):
        if not isinstance(node, ast.ImportFrom):
            continue
        for alias in node.names:
            if alias.name in forbidden_names and alias.asname is not None:
                forbidden_names.add(alias.asname)
    return frozenset(forbidden_names)


def _forbidden_setattr_field(node: ast.Call) -> str | None:
    if _qualified_name(node.func) not in {"setattr", "object.__setattr__"}:
        return None
    if len(node.args) < 2:
        return None
    field_name = _literal_string(node.args[1])
    if field_name in _FORBIDDEN_DIFFERENTIAL_PREPARATION_SETATTR_FIELDS:
        return field_name
    return None


def _assigned_forbidden_private_fields(
    node: ast.Assign | ast.AnnAssign | ast.AugAssign,
) -> tuple[str, ...]:
    targets: list[ast.expr] = []
    if isinstance(node, ast.Assign):
        targets.extend(node.targets)
    else:
        targets.append(node.target)
    return _target_forbidden_private_fields(tuple(targets))


def _deleted_forbidden_private_fields(node: ast.Delete) -> tuple[str, ...]:
    return _target_forbidden_private_fields(tuple(node.targets))


def _target_forbidden_private_fields(targets: tuple[ast.expr, ...]) -> tuple[str, ...]:
    fields: list[str] = []
    for target in targets:
        for nested in ast.walk(target):
            if (
                isinstance(nested, ast.Attribute)
                and nested.attr in _FORBIDDEN_DIFFERENTIAL_PREPARATION_PRIVATE_FIELDS
            ):
                fields.append(nested.attr)
    return tuple(fields)


def _qualified_name(node: ast.expr) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        owner = _qualified_name(node.value)
        if owner:
            return f"{owner}.{node.attr}"
        return node.attr
    return ""


def _literal_string(node: ast.expr) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None

from __future__ import annotations

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
_FORBIDDEN_DIFFERENTIAL_PREPARATION_TOKENS = (
    "ProteinAwarePreparationResult",
    "ProteinAwarePreparationReport",
    "ProteinAwarePreparationStage",
    "ProteinAwareAlignmentEligibilityResolver",
    "ProteinMappingResolver",
    "ProteinMappingConfig",
    "ProteinMappingResult",
    "ProteinMappingRecord",
)
_FORBIDDEN_DIFFERENTIAL_PREPARATION_MODULE_FRAGMENTS = (
    "phospy.science.datasets.preprocessing.stages.protein_aware_preparation",
    "phospy.science.datasets.preprocessing.protein_mapping import",
)
_FORBIDDEN_DIFFERENTIAL_PREPARATION_MUTATION_TOKENS = (
    'object.__setattr__(preparation, "_protein_covariate_matrix"',
    'object.__setattr__(preparation, "_matched_pairs"',
    "._protein_covariate_matrix",
    "._matched_pairs",
    "protein_aware_preparation=",
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


def test_differential_domains_do_not_import_or_own_protein_aware_preparation() -> None:
    violations: list[str] = []
    for directory in _DIFFERENTIAL_CODE_DIRS:
        for path in directory.rglob("*.py"):
            source = path.read_text(encoding="utf-8")
            for token in (
                _FORBIDDEN_DIFFERENTIAL_PREPARATION_TOKENS
                + _FORBIDDEN_DIFFERENTIAL_PREPARATION_MODULE_FRAGMENTS
                + _FORBIDDEN_DIFFERENTIAL_PREPARATION_MUTATION_TOKENS
            ):
                if token not in source:
                    continue
                relative_path = path.relative_to(ROOT).as_posix()
                violations.append(f"{relative_path}: contains {token!r}")

    assert not violations, (
        "protein-aware preparation and mapping implementation must stay in dataset "
        "preprocessing/building domains; differential may only consume the "
        "dataset-owned typed sidecar view:\n" + "\n".join(violations)
    )

from __future__ import annotations

import json
from typing import Literal, cast

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
)
from phospy.advanced import (
    DatasetIntensityTransformConfig,
    DatasetProteinAwarePreparationConfig,
    DifferentialAnalysisConfig,
    DifferentialProteinAwareModelConfig,
    EmpiricalBayesConfig,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    PhosphositeImportResult,
    SampleDesignRecord,
)
from phospy.api.requests import DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED
from phospy.api.results import DifferentialAnalysisResult
from phospy.errors import WorkflowBoundaryError
from phospy.io.readers import (
    FragPipeColumnMapping,
    FragPipePTMProphetImporter,
    FragPipePTMProphetImportRequest,
    MaxQuantColumnMapping,
    MaxQuantPhosphositeImporter,
    MaxQuantPhosphositeImportRequest,
)
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
    EMPIRICAL_BAYES_TREND_COVARIATE_QUANTIFICATION_DEPTH,
    QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT,
    QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
)

pytestmark = pytest.mark.integration

_SAMPLE_IDS = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
_DISPLAY_IDS = (
    "MAPK14;Y182;",
    "AKT1;T308;",
    "GSK3B;S9;",
    "MTOR;S2448;",
    "RPS6KB1;T389;",
)
_GENES = ("MAPK14", "AKT1", "GSK3B", "MTOR", "RPS6KB1")
_PROTEIN_IDS = ("P53778", "P31749", "P49841", "P42345", "P23443")
_RESIDUES = ("Y", "T", "S", "S", "T")
_POSITIONS = (182, 308, 9, 2448, 389)


def _phospho_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "A_1": [1.00, 2.05, 1.48, 0.80, 3.00],
            "A_2": [1.15, 2.10, 1.50, 0.78, 2.90],
            "A_3": [0.95, 1.92, 1.46, 0.82, 3.05],
            "B_1": [1.75, 2.48, 1.55, 1.10, 2.70],
            "B_2": [1.83, 2.57, 1.58, 1.07, 2.65],
            "B_3": [1.69, 2.41, 1.62, 1.09, 2.72],
        },
        index=pd.Index(_DISPLAY_IDS, name="site_id"),
    )


def _total_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "A_1": [10.0, 8.0, 12.5, 9.1, 6.4],
            "A_2": [10.4, 8.5, 11.9, 9.4, 6.8],
            "A_3": [9.7, 7.8, 12.8, 8.8, 6.3],
            "B_1": [11.8, 8.9, 13.0, 9.9, 7.2],
            "B_2": [12.2, 9.7, 12.7, 10.1, 7.4],
            "B_3": [11.4, 8.8, 13.5, 9.6, 7.0],
        },
        index=pd.Index(_PROTEIN_IDS, name="protein_id"),
    )


def _manual_dataset(
    *,
    quantification_depth: tuple[object, ...] | pd.Series | None = None,
    phospho: pd.DataFrame | None = None,
    protein_aware_preparation: bool = False,
) -> AnalysisReadyPhosphoDataset:
    phospho_matrix = _phospho_matrix() if phospho is None else phospho.copy(deep=True)
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": _GENES,
            "site": [
                f"{residue}{position}"
                for residue, position in zip(_RESIDUES, _POSITIONS, strict=True)
            ],
            "protein_id": _PROTEIN_IDS,
            "site_sequence": [
                ("A" * 15) + residue + ("A" * 15) for residue in _RESIDUES
            ],
            "localisation_confidence": [0.95, 0.96, 0.97, 0.94, 0.93],
        },
        index=phospho_matrix.index.copy(),
    )
    if quantification_depth is not None:
        site_metadata["quantification_depth"] = pd.Series(
            quantification_depth,
            index=phospho_matrix.index.copy(),
        )

    preprocessing_config = (
        DatasetPreprocessingConfig(
            protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                policy="prepare_model_inputs"
            )
        )
        if protein_aware_preparation
        else DatasetPreprocessingConfig()
    )
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho_matrix,
            site_metadata=site_metadata,
            total=_total_matrix() if protein_aware_preparation else None,
            organism=Organism.HUMAN,
            input_intensity_scale="log2",
            preprocessing_config=preprocessing_config,
        )
    )


def _design(*, paired: bool = False) -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=sample_id.split("_", maxsplit=1)[0],
                biological_replicate_id=f"{sample_id}_bio",
                block_id=(
                    f"subject_{sample_id.split('_', maxsplit=1)[1]}" if paired else None
                ),
            )
            for sample_id in _SAMPLE_IDS
        )
    )


def _request(
    dataset: AnalysisReadyPhosphoDataset,
    *,
    empirical_bayes: EmpiricalBayesConfig,
    paired_design_policy: Literal["reject", "duplicate_correlation"] = "reject",
    protein_aware: bool = False,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=_design(paired=paired_design_policy != "reject"),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
        config=DifferentialAnalysisConfig(
            paired_design_policy=paired_design_policy,
            empirical_bayes=empirical_bayes,
            protein_aware_model=(
                DifferentialProteinAwareModelConfig() if protein_aware else None
            ),
        ),
    )


def _depth_config(
    *,
    method: Literal["standard", "robust"] = "standard",
    kind: Literal["psm_count", "peptide_count"] = "psm_count",
) -> EmpiricalBayesConfig:
    return EmpiricalBayesConfig(
        method=method,
        trend=True,
        trend_covariate=EMPIRICAL_BAYES_TREND_COVARIATE_QUANTIFICATION_DEPTH,
        quantification_depth_kind=kind,
    )


def _assert_depth_strategy(
    result: DifferentialAnalysisResult,
    *,
    expected_depth: tuple[float, ...],
    expected_kind: str,
) -> None:
    assert result.mean_variance_trend_diagnostics is None
    diagnostics = result.quantification_depth_trend_diagnostics
    assert diagnostics is not None
    expected_index = result.residual_variance_series().index
    expected_depth_series = pd.Series(
        np.asarray(expected_depth, dtype=float),
        index=expected_index.copy(),
        name="quantification_depth",
    )
    expected_log_depth_series = pd.Series(
        np.log2(np.asarray(expected_depth, dtype=float)),
        index=expected_index.copy(),
        name="log2_quantification_depth",
    )
    assert diagnostics.trend_covariate_name == "quantification_depth"
    assert diagnostics.trend_covariate_transformation == "log2"
    assert diagnostics.quantification_depth_kind == expected_kind
    pdt.assert_series_equal(
        diagnostics.quantification_depth,
        expected_depth_series,
        check_dtype=False,
    )
    pdt.assert_series_equal(
        diagnostics.trend_covariate,
        expected_log_depth_series,
        check_dtype=False,
    )
    assert result.policy_provenance is not None
    assert result.policy_provenance.empirical_bayes.trend_covariate == (
        "quantification_depth"
    )
    assert (
        result.policy_provenance.empirical_bayes.trend_covariate_transformation
        == "log2"
    )
    assert (
        result.policy_provenance.empirical_bayes.quantification_depth_kind
        == expected_kind
    )
    payload = result.to_payload()
    policy_payload = cast(dict[str, object], payload["policy_provenance"])
    empirical_bayes_payload = cast(dict[str, object], policy_payload["empirical_bayes"])
    assert empirical_bayes_payload["trend_covariate"] == "quantification_depth"
    assert empirical_bayes_payload["trend_covariate_transformation"] == "log2"
    assert empirical_bayes_payload["quantification_depth_kind"] == expected_kind
    depth_payload = cast(
        dict[str, object],
        payload["quantification_depth_trend_diagnostics"],
    )
    assert depth_payload["trend_covariate_name"] == "quantification_depth"
    assert depth_payload["trend_covariate_transformation"] == "log2"
    assert depth_payload["quantification_depth_kind"] == expected_kind


def _assert_result_statistics_equal(
    left: DifferentialAnalysisResult,
    right: DifferentialAnalysisResult,
) -> None:
    pdt.assert_series_equal(left.residual_variance_series(), right.residual_variance)
    pdt.assert_series_equal(
        left.posterior_residual_variance_series(),
        right.posterior_residual_variance,
    )
    pdt.assert_series_equal(
        left.prior_residual_variance_series(),
        right.prior_residual_variance,
    )
    pdt.assert_series_equal(
        left.prior_degrees_of_freedom_series(),
        right.prior_degrees_of_freedom_series_value,
    )
    pdt.assert_frame_equal(left.table_for("B_vs_A"), right.table_for("B_vs_A"))


def test_ordinary_workflow_uses_manual_psm_depth_and_serializes_strategy() -> None:
    depth = (1.0, 2.0, 4.0, 8.0, 16.0)
    dataset = _manual_dataset(quantification_depth=depth)
    request = _request(dataset, empirical_bayes=_depth_config())

    result = DifferentialAnalysisWorkflow().run(request)
    repeated = DifferentialAnalysisWorkflow().run(request)

    assert result.diagnostics.variance_method == (
        "ordinary_least_squares_residual_variance"
    )
    assert result.diagnostics.moderation_method == (
        "empirical_bayes_standard_quantification_depth_trend"
    )
    _assert_depth_strategy(
        result,
        expected_depth=depth,
        expected_kind=QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
    )
    assert result.scientifically_equals(repeated)
    assert json.dumps(result.to_payload(), sort_keys=True) == json.dumps(
        repeated.to_payload(),
        sort_keys=True,
    )


def test_robust_workflow_uses_manual_peptide_depth() -> None:
    depth = (3.0, 3.0, 6.0, 12.0, 24.0)
    result = DifferentialAnalysisWorkflow().run(
        _request(
            _manual_dataset(quantification_depth=depth),
            empirical_bayes=_depth_config(
                method="robust",
                kind=QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT,
            ),
        )
    )

    assert result.empirical_bayes_method == "robust"
    assert result.empirical_bayes_robust is True
    assert result.diagnostics.moderation_method == (
        "empirical_bayes_robust_quantification_depth_trend"
    )
    assert np.isfinite(result.table_for("B_vs_A").loc[:, "P.Value"]).all()
    _assert_depth_strategy(
        result,
        expected_depth=depth,
        expected_kind=QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT,
    )


def test_duplicate_correlation_workflow_uses_manual_depth() -> None:
    depth = (1.0, 2.0, 4.0, 8.0, 16.0)
    result = DifferentialAnalysisWorkflow().run(
        _request(
            _manual_dataset(quantification_depth=depth),
            empirical_bayes=_depth_config(),
            paired_design_policy="duplicate_correlation",
        )
    )

    assert result.diagnostics.model_type == "moderated_gls_duplicate_correlation"
    assert (
        result.diagnostics.variance_method == "compound_symmetry_gls_residual_variance"
    )
    assert result.policy_provenance is not None
    assert result.policy_provenance.duplicate_correlation is not None
    _assert_depth_strategy(
        result,
        expected_depth=depth,
        expected_kind=QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
    )


def test_protein_aware_workflow_uses_manual_depth() -> None:
    depth = (2.0, 4.0, 8.0, 16.0, 32.0)
    dataset = _manual_dataset(
        quantification_depth=depth,
        protein_aware_preparation=True,
    )
    result = DifferentialAnalysisWorkflow().run(
        _request(
            dataset,
            empirical_bayes=_depth_config(),
            protein_aware=True,
        )
    )

    assert result.protein_aware_diagnostics is not None
    assert result.diagnostics.variance_method == (
        "protein_covariate_adjusted_ordinary_least_squares_residual_variance"
    )
    assert result.policy_provenance is not None
    assert result.policy_provenance.protein_aware is not None
    _assert_depth_strategy(
        result,
        expected_depth=depth,
        expected_kind=QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
    )


def test_depth_diagnostics_reexpand_withheld_features() -> None:
    phospho = _phospho_matrix()
    phospho.iloc[0, :] = 5.0
    depth = (99.0, 2.0, 4.0, 8.0, 16.0)

    result = DifferentialAnalysisWorkflow().run(
        _request(
            _manual_dataset(quantification_depth=depth, phospho=phospho),
            empirical_bayes=_depth_config(),
        )
    )

    table = result.table_for("B_vs_A")
    status = table.loc[:, DIFFERENTIAL_RESULT_STATUS_COLUMN].astype(str)
    withheld_site_key = str(result.residual_variance_series().index[0])
    assert (
        status.loc[withheld_site_key]
        == DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT
    )
    assert int((status == DIFFERENTIAL_RESULT_STATUS_TESTED).sum()) == 4
    withheld_values = table.loc[withheld_site_key]
    missing_statistics = withheld_values.loc[
        ["logFC", "t", "P.Value", "adj.P.Val"]
    ].isna()
    assert bool(np.asarray(missing_statistics, dtype=bool).all())
    diagnostics = result.quantification_depth_trend_diagnostics
    assert diagnostics is not None
    assert diagnostics.quantification_depth.index.equals(result.residual_variance.index)
    assert np.isnan(float(diagnostics.quantification_depth.loc[withheld_site_key]))
    assert diagnostics.quantification_depth.drop(index=withheld_site_key).notna().all()
    assert diagnostics.trend_covariate.drop(index=withheld_site_key).notna().all()


def test_ordinary_depth_ignores_invalid_depth_on_pre_excluded_site() -> None:
    phospho = _phospho_matrix()
    phospho.iloc[0, :] = 5.0
    depth_with_invalid_excluded_site = (np.nan, 2.0, 4.0, 8.0, 16.0)
    depth_with_valid_excluded_site = (99.0, 2.0, 4.0, 8.0, 16.0)

    result = DifferentialAnalysisWorkflow().run(
        _request(
            _manual_dataset(
                quantification_depth=depth_with_invalid_excluded_site,
                phospho=phospho,
            ),
            empirical_bayes=_depth_config(),
        )
    )
    clean_baseline = DifferentialAnalysisWorkflow().run(
        _request(
            _manual_dataset(
                quantification_depth=depth_with_valid_excluded_site,
                phospho=phospho,
            ),
            empirical_bayes=_depth_config(),
        )
    )

    _assert_result_statistics_equal(result, clean_baseline)
    table = result.table_for("B_vs_A")
    status = table.loc[:, DIFFERENTIAL_RESULT_STATUS_COLUMN].astype(str)
    withheld_site_key = str(result.residual_variance_series().index[0])
    tested = status == DIFFERENTIAL_RESULT_STATUS_TESTED

    assert (
        status.loc[withheld_site_key]
        == DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT
    )
    assert "all-constant" in str(
        table.loc[withheld_site_key, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN]
    )
    assert int(tested.sum()) == 4
    assert (
        np.isfinite(table.loc[tested, ["logFC", "t", "P.Value", "adj.P.Val"]])
        .all()
        .all()
    )
    assert (
        table.loc[withheld_site_key, ["logFC", "t", "P.Value", "adj.P.Val"]]
        .isna()
        .all()
    )
    for series in (
        result.residual_variance_series(),
        result.posterior_residual_variance_series(),
        result.prior_residual_variance_series(),
        result.prior_degrees_of_freedom_series(),
    ):
        assert np.isnan(float(series.loc[withheld_site_key]))

    assert result.mean_variance_trend_diagnostics is None
    diagnostics = result.quantification_depth_trend_diagnostics
    assert diagnostics is not None
    expected_tested_depth = np.asarray([2.0, 4.0, 8.0, 16.0], dtype=float)
    np.testing.assert_allclose(
        diagnostics.quantification_depth.loc[tested].to_numpy(dtype=float),
        expected_tested_depth,
    )
    np.testing.assert_allclose(
        diagnostics.trend_covariate.loc[tested].to_numpy(dtype=float),
        np.log2(expected_tested_depth),
    )
    assert np.isnan(float(diagnostics.quantification_depth.loc[withheld_site_key]))
    assert np.isnan(float(diagnostics.trend_covariate.loc[withheld_site_key]))
    assert np.isnan(float(diagnostics.log_residual_variance.loc[withheld_site_key]))
    assert np.isnan(float(diagnostics.fitted_log_prior_variance.loc[withheld_site_key]))
    assert diagnostics.quantification_depth.drop(index=withheld_site_key).notna().all()
    assert diagnostics.trend_covariate.drop(index=withheld_site_key).notna().all()
    assert result.diagnostics.moderation_method == (
        "empirical_bayes_standard_quantification_depth_trend"
    )
    assert result.policy_provenance is not None
    assert result.policy_provenance.empirical_bayes.trend_covariate == (
        "quantification_depth"
    )
    assert (
        result.policy_provenance.empirical_bayes.trend_covariate_transformation
        == "log2"
    )
    assert result.policy_provenance.empirical_bayes.quantification_depth_kind == (
        QUANTIFICATION_DEPTH_KIND_PSM_COUNT
    )


def test_ordinary_depth_rejects_invalid_depth_on_tested_site() -> None:
    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.interpreter.quantification_depth",
    ) as exc_info:
        DifferentialAnalysisWorkflow().run(
            _request(
                _manual_dataset(
                    quantification_depth=(1.0, np.nan, 4.0, 8.0, 16.0),
                ),
                empirical_bayes=_depth_config(),
            )
        )

    error = str(exc_info.value)
    assert "quantification_depth" in error
    assert "missing values" in error


def test_existing_global_and_mean_intensity_modes_ignore_depth_metadata() -> None:
    depth = (1.0, 2.0, 4.0, 8.0, 16.0)
    dataset_without_depth = _manual_dataset()
    dataset_with_depth = _manual_dataset(quantification_depth=depth)

    global_without_depth = DifferentialAnalysisWorkflow().run(
        _request(dataset_without_depth, empirical_bayes=EmpiricalBayesConfig())
    )
    global_with_depth = DifferentialAnalysisWorkflow().run(
        _request(dataset_with_depth, empirical_bayes=EmpiricalBayesConfig())
    )
    assert global_with_depth.mean_variance_trend_diagnostics is None
    assert global_with_depth.quantification_depth_trend_diagnostics is None
    assert global_with_depth.policy_provenance is not None
    assert global_with_depth.policy_provenance.empirical_bayes.trend_covariate is None
    _assert_result_statistics_equal(global_without_depth, global_with_depth)

    trend_without_depth = DifferentialAnalysisWorkflow().run(
        _request(
            dataset_without_depth, empirical_bayes=EmpiricalBayesConfig(trend=True)
        )
    )
    trend_with_depth = DifferentialAnalysisWorkflow().run(
        _request(dataset_with_depth, empirical_bayes=EmpiricalBayesConfig(trend=True))
    )
    assert trend_with_depth.mean_variance_trend_diagnostics is not None
    assert trend_with_depth.quantification_depth_trend_diagnostics is None
    assert trend_with_depth.policy_provenance is not None
    assert trend_with_depth.policy_provenance.empirical_bayes.trend_covariate == (
        "mean_intensity"
    )
    assert (
        trend_with_depth.policy_provenance.empirical_bayes.trend_covariate_transformation
        == "identity"
    )
    assert (
        trend_with_depth.policy_provenance.empirical_bayes.quantification_depth_kind
        is None
    )
    _assert_result_statistics_equal(trend_without_depth, trend_with_depth)


def test_maxquant_imported_depth_feeds_depth_aware_workflow() -> None:
    depth = (3.0, 4.0, 5.0, 6.0, 7.0)
    unmapped = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(source=_maxquant_depth_source(depth))
    )
    assert "quantification_depth" not in unmapped.site_metadata_candidate.columns

    imported = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=_maxquant_depth_source(depth),
            column_mapping=MaxQuantColumnMapping(
                quantification_depth="Depth",
                quantification_depth_kind=QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
            ),
        )
    )
    dataset = _dataset_from_imported_depth(imported)
    result = DifferentialAnalysisWorkflow().run(
        _request(dataset, empirical_bayes=_depth_config())
    )

    _assert_depth_strategy(
        result,
        expected_depth=depth,
        expected_kind=QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
    )


def test_fragpipe_imported_depth_feeds_depth_aware_workflow() -> None:
    depth = (3.0, 4.0, 5.0, 6.0, 7.0)
    unmapped = FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=_fragpipe_depth_source(depth),
            ptmprophet_position_reference="protein",
        )
    )
    assert "quantification_depth" not in unmapped.site_metadata_candidate.columns

    imported = FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=_fragpipe_depth_source(depth),
            column_mapping=FragPipeColumnMapping(
                quantification_depth="Depth",
                quantification_depth_kind=QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
            ),
            ptmprophet_position_reference="protein",
        )
    )
    dataset = _dataset_from_imported_depth(imported)
    result = DifferentialAnalysisWorkflow().run(
        _request(dataset, empirical_bayes=_depth_config())
    )

    _assert_depth_strategy(
        result,
        expected_depth=depth,
        expected_kind=QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
    )


def _dataset_from_imported_depth(
    import_result: PhosphositeImportResult,
) -> AnalysisReadyPhosphoDataset:
    return AnalysisReadyDatasetBuilder().run(
        import_result.to_dataset_build_request(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            ),
        )
    )


def _maxquant_depth_source(depth: tuple[float, ...]) -> pd.DataFrame:
    rows = len(depth)
    frame = _importer_base_frame(depth)
    frame["Proteins"] = list(_PROTEIN_IDS[:rows])
    frame["Gene names"] = list(_GENES[:rows])
    frame["Amino acid"] = list(_RESIDUES[:rows])
    frame["Positions within proteins"] = [str(value) for value in _POSITIONS[:rows]]
    frame["Localization prob"] = ["0.95"] * rows
    frame["Sequence"] = [f"AAAAA{residue}AAAA" for residue in _RESIDUES[:rows]]
    frame["Modified sequence"] = [
        f"AAAAA(ph){residue}AAAA" for residue in _RESIDUES[:rows]
    ]
    frame["Sequence window"] = [
        ("A" * 15) + residue + ("A" * 15) for residue in _RESIDUES[:rows]
    ]
    frame["Potential contaminant"] = [""] * rows
    frame["Reverse"] = [""] * rows
    return frame


def _fragpipe_depth_source(depth: tuple[float, ...]) -> pd.DataFrame:
    rows = len(depth)
    frame = _importer_base_frame(depth)
    frame["Protein"] = list(_PROTEIN_IDS[:rows])
    frame["Gene"] = list(_GENES[:rows])
    frame["Peptide"] = [f"AAAAA{residue}AAAA" for residue in _RESIDUES[:rows]]
    frame["Modified Peptide"] = [
        f"AAAAA[p{residue}]AAAA" for residue in _RESIDUES[:rows]
    ]
    frame["PTMProphet Probability"] = [
        f"{residue}{position}(0.95)"
        for residue, position in zip(
            _RESIDUES[:rows],
            _POSITIONS[:rows],
            strict=True,
        )
    ]
    frame["Site"] = [
        f"{residue}{position}"
        for residue, position in zip(
            _RESIDUES[:rows],
            _POSITIONS[:rows],
            strict=True,
        )
    ]
    frame["Sequence Window"] = [
        ("A" * 15) + residue + ("A" * 15) for residue in _RESIDUES[:rows]
    ]
    frame["Spectrum"] = [f"scan.{row + 1}.{row + 1}.2" for row in range(rows)]
    frame["PSM ID"] = [f"psm-{row + 1}" for row in range(rows)]
    return frame


def _importer_base_frame(depth: tuple[float, ...]) -> pd.DataFrame:
    rows = len(depth)
    values: dict[str, list[float]] = {
        "Depth": [float(value) for value in depth],
    }
    for sample_position, sample_id in enumerate(_SAMPLE_IDS):
        values[f"Intensity {sample_id}"] = [
            10.0 + float(row) + float(sample_position) * 2.0 for row in range(rows)
        ]
    return pd.DataFrame(values)

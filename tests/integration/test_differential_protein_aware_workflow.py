from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
)
from phospy.advanced import (
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    DatasetProteinAwarePreparationConfig,
    DifferentialAnalysisConfig,
    DifferentialProteinAwareModelConfig,
    TechnicalReplicatePolicy,
)
from phospy.api import (
    ContinuousCovariate,
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.api.results import DifferentialAnalysisResult
from phospy.errors import WorkflowBoundaryError
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
)
from phospy.science.statistics.multiple_testing import adjust_p_values
from tests.support.unsafe_dataset_states import (
    unsafe_mark_dataset_total_protein_correction_applied,
)

pytestmark = pytest.mark.integration


_SAMPLE_IDS = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
_DISPLAY_IDS = ("MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;")
_PROTEIN_IDS = ("P53778", "P31749", "P49841")
_TOTAL_VALUES = {
    "P53778": {
        "A_1": 10.0,
        "A_2": 10.4,
        "A_3": 9.7,
        "B_1": 11.8,
        "B_2": 12.2,
        "B_3": 11.4,
    },
    "P31749": {
        "A_1": 8.0,
        "A_2": 8.5,
        "A_3": 7.8,
        "B_1": 8.9,
        "B_2": 9.7,
        "B_3": 8.8,
    },
    "P49841": {
        "A_1": 12.5,
        "A_2": 11.9,
        "A_3": 12.8,
        "B_1": 13.0,
        "B_2": 12.7,
        "B_3": 13.5,
    },
}


def _build_dataset(
    *,
    protein_aware_preparation: bool = True,
    constant_total_proteins: frozenset[str] = frozenset(),
) -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {
            "A_1": [1.00, 2.05, 1.48],
            "A_2": [1.15, 2.10, 1.50],
            "A_3": [0.95, 1.92, 1.46],
            "B_1": [1.75, 2.48, 1.55],
            "B_2": [1.83, 2.57, 1.58],
            "B_3": [1.69, 2.41, 1.62],
        },
        index=pd.Index(_DISPLAY_IDS, name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ("MAPK14", "AKT1", "GSK3B"),
            "site": ("Y182", "T308", "S9"),
            "protein_id": _PROTEIN_IDS,
            "site_sequence": (
                "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            ),
            "localisation_confidence": (0.95, 0.96, 0.97),
        },
        index=phospho.index.copy(),
    )
    total = pd.DataFrame(
        {
            sample_id: [
                (
                    14.0 + float(position) * 0.3
                    if protein_id in constant_total_proteins
                    else _TOTAL_VALUES[protein_id][sample_id]
                )
                for position, protein_id in enumerate(_PROTEIN_IDS)
            ]
            for sample_id in _SAMPLE_IDS
        },
        index=pd.Index(_PROTEIN_IDS, name="protein_id"),
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
            phospho=phospho,
            site_metadata=site_metadata,
            total=total,
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=preprocessing_config,
        )
    )


def _contrasts(*, reciprocal: bool = False) -> tuple[Contrast, ...]:
    contrasts = [
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        )
    ]
    if reciprocal:
        contrasts.append(
            Contrast(
                name="A_vs_B",
                numerator_condition="A",
                denominator_condition="B",
            )
        )
    return tuple(contrasts)


def _design(
    *,
    sample_ids: tuple[str, ...] = _SAMPLE_IDS,
    paired_design_policy: str = "reject",
    fixed_effects: tuple[ContinuousCovariate, ...] = (),
    sample_covariates: dict[str, dict[str, float]] | None = None,
    technical_replicates: bool = False,
) -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=sample_id.split("_", maxsplit=1)[0],
                biological_replicate_id=_biological_replicate_id(
                    sample_id,
                    technical_replicates=technical_replicates,
                ),
                technical_replicate_id=(
                    _technical_replicate_id(sample_id) if technical_replicates else None
                ),
                block_id=(
                    f"pair_{sample_id.split('_', maxsplit=1)[1]}"
                    if paired_design_policy
                    in {
                        PAIRED_DESIGN_POLICY_FIXED_BLOCK,
                        PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
                    }
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
    )


def _biological_replicate_id(
    sample_id: str,
    *,
    technical_replicates: bool,
) -> str:
    if not technical_replicates:
        return f"{sample_id}_bio"
    condition, replicate = sample_id.split("_", maxsplit=1)
    if replicate in {"1", "2"}:
        return f"{condition}_bio_technical_pair"
    return f"{sample_id}_bio"


def _technical_replicate_id(sample_id: str) -> str:
    replicate = sample_id.split("_", maxsplit=1)[1]
    return f"tech_{replicate}"


def _request(
    dataset: AnalysisReadyPhosphoDataset,
    *,
    sample_ids: tuple[str, ...] = _SAMPLE_IDS,
    paired_design_policy: str = "reject",
    allow_design_subset: bool = False,
    fixed_effects: tuple[ContinuousCovariate, ...] = (),
    sample_covariates: dict[str, dict[str, float]] | None = None,
    technical_replicate_policy: TechnicalReplicatePolicy = TechnicalReplicatePolicy.REJECT,
    technical_replicates: bool = False,
    reciprocal_contrast: bool = False,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=_design(
            sample_ids=sample_ids,
            paired_design_policy=paired_design_policy,
            fixed_effects=fixed_effects,
            sample_covariates=sample_covariates,
            technical_replicates=technical_replicates,
        ),
        contrasts=_contrasts(reciprocal=reciprocal_contrast),
        config=DifferentialAnalysisConfig(
            paired_design_policy=paired_design_policy,  # type: ignore[arg-type]
            allow_design_subset=allow_design_subset,
            technical_replicate_policy=technical_replicate_policy,
            protein_aware_model=DifferentialProteinAwareModelConfig(),
        ),
    )


def _ordinary_request(
    dataset: AnalysisReadyPhosphoDataset,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=_design(),
        contrasts=_contrasts(),
    )


def _assert_adjusted_public_result(
    result: DifferentialAnalysisResult,
    *,
    contrast_names: tuple[str, ...] = ("B_vs_A",),
    tested_site_count: int = 3,
    withheld_site_count: int = 0,
) -> None:
    if not hasattr(result, "protein_aware_diagnostics"):
        raise AssertionError("result must expose protein-aware diagnostics")
    diagnostics = result.protein_aware_diagnostics
    assert diagnostics is not None
    assert diagnostics.method_id == DifferentialProteinAwareModelConfig().method
    assert diagnostics.tested_site_count == tested_site_count
    assert diagnostics.withheld_site_count == withheld_site_count
    assert diagnostics.protein_covariate_centered is True
    assert diagnostics.protein_covariate_standardized is False
    assert diagnostics.fallback_policy == "no_fallback_to_ordinary_differential_lane"
    assert result.policy_provenance is not None
    assert result.policy_provenance.protein_aware is not None
    assert result.workflow_provenance is not None
    row_attrition = result.workflow_provenance["row_attrition_metrics"]
    assert row_attrition == {
        "input_sites": tested_site_count + withheld_site_count,
        "sites_retained_for_model_fitting": tested_site_count,
        "sites_excluded_before_testing": withheld_site_count,
        "sites_with_failed_model_fit": 0,
        "sites_included_in_multiple_testing_family": tested_site_count,
    }

    for contrast_name in contrast_names:
        table = result.table_for(contrast_name)
        tested = table[DIFFERENTIAL_RESULT_STATUS_COLUMN].astype(str) == (
            DIFFERENTIAL_RESULT_STATUS_TESTED
        )
        assert int(tested.sum()) == tested_site_count
        assert (
            table.loc[tested, ["logFC", "t", "P.Value", "adj.P.Val"]]
            .notna()
            .all()
            .all()
        )
        assert np.isfinite(
            table.loc[tested, ["logFC", "t", "P.Value", "adj.P.Val"]].to_numpy(
                dtype=float
            )
        ).all()


def test_differential_protein_aware_public_workflow_runs_adjusted_multiple_contrasts() -> (
    None
):
    dataset = _build_dataset()

    result = DifferentialAnalysisWorkflow().run(
        _request(dataset, reciprocal_contrast=True)
    )

    _assert_adjusted_public_result(
        result,
        contrast_names=("B_vs_A", "A_vs_B"),
    )
    assert set(result.contrast_tables) == {"B_vs_A", "A_vs_B"}


def test_differential_protein_aware_public_workflow_supports_covariate_and_block() -> (
    None
):
    dataset = _build_dataset()
    dose_covariates = {
        "A_1": {"dose": 0.0},
        "A_2": {"dose": 1.5},
        "A_3": {"dose": 3.0},
        "B_1": {"dose": 0.2},
        "B_2": {"dose": 1.7},
        "B_3": {"dose": 3.2},
    }

    covariate_result = DifferentialAnalysisWorkflow().run(
        _request(
            dataset,
            fixed_effects=(ContinuousCovariate("dose"),),
            sample_covariates=dose_covariates,
        )
    )
    block_result = DifferentialAnalysisWorkflow().run(
        _request(
            dataset,
            paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
        )
    )

    _assert_adjusted_public_result(covariate_result)
    _assert_adjusted_public_result(block_result)


def test_differential_protein_aware_public_workflow_respects_subset_reorder() -> None:
    dataset = _build_dataset()
    sample_ids = ("B_3", "B_1", "A_2", "A_1")

    result = DifferentialAnalysisWorkflow().run(
        _request(
            dataset,
            sample_ids=sample_ids,
            allow_design_subset=True,
        )
    )
    repeated = DifferentialAnalysisWorkflow().run(
        _request(
            dataset,
            sample_ids=sample_ids,
            allow_design_subset=True,
        )
    )

    _assert_adjusted_public_result(result)
    assert result.protein_aware_diagnostics is not None
    assert result.protein_aware_diagnostics.execution_sample_order == sample_ids
    pd.testing.assert_frame_equal(
        result.table_for("B_vs_A"),
        repeated.table_for("B_vs_A"),
    )


def test_differential_protein_aware_public_workflow_retains_withheld_rows() -> None:
    dataset = _build_dataset(constant_total_proteins=frozenset({"P49841"}))

    result = DifferentialAnalysisWorkflow().run(_request(dataset))

    _assert_adjusted_public_result(
        result,
        tested_site_count=2,
        withheld_site_count=1,
    )
    table = result.table_for("B_vs_A")
    withheld = table[DIFFERENTIAL_RESULT_STATUS_COLUMN].astype(str) == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID
    )
    assert int(withheld.sum()) == 1
    assert (
        table.loc[withheld, ["logFC", "t", "P.Value", "adj.P.Val"]].isna().all().all()
    )
    assert (
        table.loc[withheld, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN]
        .astype(str)
        .str.contains("protein_covariate_zero_variance")
        .all()
    )
    expected_adjusted = adjust_p_values(
        table.loc[~withheld, "P.Value"].to_numpy(dtype=float),
        method="benjamini_hochberg",
    )
    np.testing.assert_allclose(
        table.loc[~withheld, "adj.P.Val"].to_numpy(dtype=float),
        expected_adjusted,
        rtol=1e-12,
        atol=1e-12,
    )


def test_differential_protein_aware_public_workflow_rejects_missing_sidecar() -> None:
    dataset = _build_dataset(protein_aware_preparation=False)

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.sidecar_missing",
    ):
        DifferentialAnalysisWorkflow().run(_request(dataset))


def test_differential_protein_aware_public_workflow_rejects_prior_subtraction() -> None:
    dataset = _build_dataset()
    unsafe_mark_dataset_total_protein_correction_applied(dataset)

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.prior_total_protein_subtraction",
    ):
        DifferentialAnalysisWorkflow().run(_request(dataset))


def test_differential_protein_aware_public_workflow_rejects_duplicate_correlation() -> (
    None
):
    dataset = _build_dataset()

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.duplicate_correlation",
    ):
        DifferentialAnalysisWorkflow().run(
            _request(
                dataset,
                paired_design_policy=PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
            )
        )


def test_differential_protein_aware_public_workflow_rejects_technical_aggregation() -> (
    None
):
    dataset = _build_dataset()

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.technical_replicate_aggregation",
    ):
        DifferentialAnalysisWorkflow().run(
            _request(
                dataset,
                technical_replicate_policy=TechnicalReplicatePolicy.MEAN,
                technical_replicates=True,
            )
        )


def test_differential_protein_aware_public_workflow_rejects_all_withheld() -> None:
    dataset = _build_dataset(constant_total_proteins=frozenset(_PROTEIN_IDS))

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.all_sites_withheld",
    ):
        DifferentialAnalysisWorkflow().run(_request(dataset))


def test_differential_protein_aware_ordinary_workflow_ignores_sidecar() -> None:
    dataset = _build_dataset()
    dataset_without_sidecar = _build_dataset(protein_aware_preparation=False)

    result = DifferentialAnalysisWorkflow().run(_ordinary_request(dataset))
    result_without_sidecar = DifferentialAnalysisWorkflow().run(
        _ordinary_request(dataset_without_sidecar)
    )

    assert result.protein_aware_diagnostics is None
    assert result.policy_provenance is not None
    assert result.policy_provenance.protein_aware is None
    pd.testing.assert_frame_equal(
        result.table_for("B_vs_A"),
        result_without_sidecar.table_for("B_vs_A"),
    )
    pd.testing.assert_series_equal(
        result.residual_variance_series(),
        result_without_sidecar.residual_variance_series(),
    )


def test_differential_protein_aware_request_has_no_caller_sidecar_injection() -> None:
    request_type: Any = DifferentialAnalysisRequest
    dataset = _build_dataset()
    design = _design()
    contrasts = _contrasts()

    with pytest.raises(TypeError):
        request_type(
            dataset=dataset,
            design=design,
            contrasts=contrasts,
            protein_aware_preparation=object(),
        )
    with pytest.raises(TypeError):
        request_type(
            dataset=dataset,
            design=design,
            contrasts=contrasts,
            protein_covariate_matrix=pd.DataFrame(),
        )

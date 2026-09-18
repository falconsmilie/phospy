from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.advanced import (
    CorrectionMissingnessPolicy,
    DifferentialAnalysisConfig,
    SpsDiscoveryConfig,
    SpsDiscoveryRequest,
    SpsDiscoveryWorkflow,
    SpsReferenceDataset,
    SpsRuvBatchCorrectionConfig,
    TechnicalReplicatePolicy,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from tests.support.site_keys import protein_site_key_index

pytestmark = pytest.mark.integration

_STABLE_SITE_COUNT = 8
_CHANGING_SITE_COUNT = 12
_SITE_COUNT = _STABLE_SITE_COUNT + _CHANGING_SITE_COUNT
_SAMPLES = (
    "A_bio1_tech1",
    "A_bio1_tech2",
    "A_bio2_tech1",
    "A_bio2_tech2",
    "B_bio1_tech1",
    "B_bio1_tech2",
    "B_bio2_tech1",
    "B_bio2_tech2",
)
_CONDITIONS = ("A", "A", "A", "A", "B", "B", "B", "B")
_REPLICATES = (
    "A_bio1",
    "A_bio1",
    "A_bio2",
    "A_bio2",
    "B_bio1",
    "B_bio1",
    "B_bio2",
    "B_bio2",
)


@dataclass(frozen=True, slots=True)
class _SyntheticTarget:
    absolute_log2_abundance: pd.DataFrame
    biological_truth: pd.DataFrame
    condition_effect_by_site: pd.Series


def _site_keys() -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=[
            f"P{position:04d}" for position in range(1, _SITE_COUNT + 1)
        ],
        sites=[f"S{position}" for position in range(1, _SITE_COUNT + 1)],
    )


def _synthetic_sps_reference_evidence() -> tuple[SpsReferenceDataset, ...]:
    """Build governed SPS evidence from absolute log2 reference abundance.

    The explicit baseline subtraction in this helper is the scientific boundary
    between processed reference abundance and condition-relative SPS evidence.
    It is deliberately unrelated to ``_synthetic_target_experiment``.
    """

    references: list[SpsReferenceDataset] = []
    site_position = np.arange(_SITE_COUNT, dtype=float)
    stable_change = 0.01 + 0.002 * (site_position[:_STABLE_SITE_COUNT] % 3.0)
    changing_change = 1.5 + 0.15 * (site_position[_STABLE_SITE_COUNT:] % 5.0)
    known_change = np.concatenate((stable_change, changing_change))

    for reference_number in range(3):
        prefix = f"reference_{reference_number + 1}"
        samples = (
            f"{prefix}_baseline_1",
            f"{prefix}_baseline_2",
            f"{prefix}_stimulated_1",
            f"{prefix}_stimulated_2",
        )
        baseline = 18.0 + 0.2 * site_position + 0.1 * reference_number
        replicate_noise = np.vstack(
            (
                0.004 * ((site_position + reference_number) % 3.0 - 1.0),
                -0.003 * ((site_position + 2.0 * reference_number) % 4.0 - 1.5),
                0.005 * ((site_position + reference_number) % 5.0 - 2.0),
                -0.004 * ((site_position + reference_number) % 6.0 - 2.5),
            )
        ).T
        source_absolute_log2 = pd.DataFrame(
            baseline[:, None]
            + replicate_noise
            + np.asarray((0.0, 0.0, 1.0, 1.0))[None, :]
            * (known_change + 0.01 * reference_number)[:, None],
            index=_site_keys(),
            columns=samples,
        )
        if reference_number == 1:
            # Controlled reference missingness: the site remains observed in both
            # conditions and therefore still has scientifically usable evidence.
            source_absolute_log2.iloc[2, 3] = np.nan

        baseline_mean = source_absolute_log2.loc[:, samples[:2]].mean(axis=1)
        condition_relative_log2 = source_absolute_log2.subtract(
            baseline_mean,
            axis="index",
        )
        np.testing.assert_allclose(
            condition_relative_log2.loc[:, samples[:2]].mean(axis=1),
            0.0,
            atol=1e-12,
        )
        references.append(
            SpsReferenceDataset.from_condition_relative_log2(
                dataset_id=prefix,
                intensities=condition_relative_log2,
                condition_by_sample={
                    samples[0]: "baseline",
                    samples[1]: "baseline",
                    samples[2]: "stimulated",
                    samples[3]: "stimulated",
                },
                log2_scale_established_by=(
                    f"{prefix} synthetic processed log2 abundance generator"
                ),
                baseline_centering_established_by=(
                    f"{prefix} explicit baseline-condition mean subtraction"
                ),
                organism="rat",
                baseline_context="synthetic baseline condition",
                reference_context="synthetic rat SPS/RUV benchmark",
                source_name="deterministic synthetic SPS reference evidence",
                source_version="PHOSPY-RUV-07-v1",
                source_uri=f"https://example.test/sps-ruv/{prefix}",
            )
        )
    return tuple(references)


def _discover_controls():
    return SpsDiscoveryWorkflow().run(
        SpsDiscoveryRequest(
            reference_datasets=_synthetic_sps_reference_evidence(),
            config=SpsDiscoveryConfig(
                top_n=_STABLE_SITE_COUNT,
                minimum_reference_datasets=3,
                minimum_datasets_per_site=3,
                minimum_shared_sites=_SITE_COUNT,
            ),
        )
    )


def _synthetic_target_experiment() -> _SyntheticTarget:
    """Build absolute target abundance without SPS-reference semantics."""

    site_position = np.arange(_SITE_COUNT, dtype=float)
    condition_effect = np.concatenate(
        (
            np.zeros(_STABLE_SITE_COUNT),
            np.where(
                np.arange(_CHANGING_SITE_COUNT) % 2 == 0,
                2.0,
                -1.5,
            ),
        )
    )
    condition_indicator = np.asarray((0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0))
    biological_replicate = np.asarray((-0.2, -0.2, 0.2, 0.2, -0.15, -0.15, 0.15, 0.15))
    replicate_loading = 0.2 + 0.05 * (site_position % 4.0)
    base_abundance = 14.0 + 0.12 * site_position
    representative_noise = 0.006 * np.sin(
        (site_position[:, None] + 1.0)
        * (np.arange(len(_SAMPLES), dtype=float)[None, :] + 1.0)
    )
    biological_truth = (
        base_abundance[:, None]
        + condition_effect[:, None] * condition_indicator[None, :]
        + replicate_loading[:, None] * biological_replicate[None, :]
        + representative_noise
    )

    unwanted_factor_1 = np.asarray((-2.0, 2.0, -1.0, 1.0, -3.0, 3.0, -1.5, 1.5))
    unwanted_factor_2 = np.asarray((-0.4, 0.4, 1.5, -1.5, -0.8, 0.8, 2.2, -2.2))
    loading_1 = 0.55 + 0.11 * (site_position % 7.0)
    loading_2 = -0.35 + 0.13 * (site_position % 6.0)
    unwanted = (
        loading_1[:, None] * unwanted_factor_1[None, :]
        + loading_2[:, None] * unwanted_factor_2[None, :]
    )

    keys = _site_keys()
    return _SyntheticTarget(
        absolute_log2_abundance=pd.DataFrame(
            biological_truth + unwanted,
            index=keys,
            columns=_SAMPLES,
        ),
        biological_truth=pd.DataFrame(
            biological_truth,
            index=keys,
            columns=_SAMPLES,
        ),
        condition_effect_by_site=pd.Series(
            condition_effect,
            index=keys,
            name="known_B_minus_A",
        ),
    )


def _site_metadata() -> pd.DataFrame:
    display_ids = [
        f"P{position:04d};S{position};" for position in range(1, _SITE_COUNT + 1)
    ]
    return pd.DataFrame(
        {
            "gene_symbol": [
                f"G{position:04d}" for position in range(1, _SITE_COUNT + 1)
            ],
            "protein_id": [
                f"P{position:04d}" for position in range(1, _SITE_COUNT + 1)
            ],
            "site": [f"S{position}" for position in range(1, _SITE_COUNT + 1)],
            "site_sequence": [("A" * 15) + "S" + ("A" * 15)] * _SITE_COUNT,
            "localisation_confidence": [0.99] * _SITE_COUNT,
        },
        index=pd.Index(display_ids, name="site_id"),
    )


def _sample_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ("tech1", "tech2") * 4,
            "condition": _CONDITIONS,
            "replicate": _REPLICATES,
        },
        index=pd.Index(_SAMPLES, name="sample_id"),
    )


def _build_target(
    *,
    method: str | None,
    k: int = 2,
) -> AnalysisReadyPhosphoDataset:
    target = _synthetic_target_experiment().absolute_log2_abundance
    target_for_builder = target.copy(deep=True)
    target_for_builder.index = _site_metadata().index.copy()
    preprocessing = DatasetPreprocessingConfig()
    if method is not None:
        preprocessing = DatasetPreprocessingConfig(
            batch_correction=SpsRuvBatchCorrectionConfig(
                control_site_set=_discover_controls().control_site_set,
                batch_column="batch",
                condition_columns=("condition",),
                replicate_column=("replicate" if method == "ruv_iii_style" else None),
                missingness_policy=CorrectionMissingnessPolicy(),
                n_unwanted_factors=k,
                method=method,
            )
        )
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=target_for_builder,
            site_metadata=_site_metadata(),
            sample_metadata=_sample_metadata(),
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=preprocessing,
        )
    )


def _rmse_from_known_biology(matrix: pd.DataFrame) -> float:
    truth = _synthetic_target_experiment().biological_truth
    # RUV factors are identifiable only up to feature-wise location shifts. The
    # scientific target is therefore the sample-relative signal for each site.
    centered_matrix = matrix.subtract(matrix.mean(axis=1), axis="index")
    centered_truth = truth.subtract(truth.mean(axis=1), axis="index")
    return float(
        np.sqrt(
            np.mean(np.square(centered_matrix.to_numpy() - centered_truth.to_numpy()))
        )
    )


def _condition_effect(matrix: pd.DataFrame) -> pd.Series:
    return matrix.loc[:, list(_SAMPLES[4:])].mean(axis=1) - matrix.loc[
        :, list(_SAMPLES[:4])
    ].mean(axis=1)


def _within_replicate_disagreement(matrix: pd.DataFrame) -> float:
    pair_differences = [
        matrix.iloc[:, right].to_numpy() - matrix.iloc[:, left].to_numpy()
        for left, right in ((0, 1), (2, 3), (4, 5), (6, 7))
    ]
    return float(np.sqrt(np.mean(np.square(np.concatenate(pair_differences)))))


def test_separate_reference_evidence_discovers_the_known_stable_population() -> None:
    discovery = _discover_controls()
    stable_keys = tuple(str(value) for value in _site_keys()[:_STABLE_SITE_COUNT])
    changing_keys = set(str(value) for value in _site_keys()[_STABLE_SITE_COUNT:])

    assert set(discovery.selected_site_keys) == set(stable_keys)
    assert not changing_keys.intersection(discovery.selected_site_keys)
    assert set(
        record.site_key for record in discovery.site_ranking[:_STABLE_SITE_COUNT]
    ) == set(stable_keys)
    assert min(
        record.consensus_stability_score
        for record in discovery.site_ranking[:_STABLE_SITE_COUNT]
    ) > max(
        record.consensus_stability_score
        for record in discovery.site_ranking[_STABLE_SITE_COUNT:]
    )
    assert all(
        source.quantitative_meaning.value == "contrast_log2_fold_change"
        for source in discovery.provenance.source_datasets
    )
    assert discovery == _discover_controls()


def test_ruv_iii_reduces_planted_technical_variation_and_preserves_biology() -> None:
    uncorrected = _build_target(method=None)
    native = _build_target(method="sps_ruv_style")
    ruv_k1 = _build_target(method="ruv_iii_style", k=1)
    ruv_k2 = _build_target(method="ruv_iii_style", k=2)

    raw_rmse = _rmse_from_known_biology(uncorrected.phospho)
    k1_rmse = _rmse_from_known_biology(ruv_k1.phospho)
    k2_rmse = _rmse_from_known_biology(ruv_k2.phospho)
    assert k2_rmse < k1_rmse < raw_rmse
    assert k2_rmse < 0.1 * raw_rmse
    assert _within_replicate_disagreement(ruv_k2.phospho) < (
        0.1 * _within_replicate_disagreement(uncorrected.phospho)
    )

    known_effect = _synthetic_target_experiment().condition_effect_by_site
    for corrected in (native.phospho, ruv_k2.phospho):
        observed_effect = _condition_effect(corrected)
        np.testing.assert_allclose(
            observed_effect.iloc[_STABLE_SITE_COUNT:],
            known_effect.iloc[_STABLE_SITE_COUNT:],
            atol=0.02,
        )
        assert float(observed_effect.iloc[:_STABLE_SITE_COUNT].abs().max()) < 0.02

    # The historical method remains a distinct, useful and deterministic
    # scientific baseline; this scenario does not require it to equal RUV-III.
    assert _rmse_from_known_biology(native.phospho) < raw_rmse
    pdt.assert_frame_equal(
        native.phospho,
        _build_target(method="sps_ruv_style").phospho,
    )
    pdt.assert_frame_equal(
        ruv_k2.phospho,
        _build_target(method="ruv_iii_style", k=2).phospho,
    )


def test_corrected_target_keeps_dataset_state_and_runs_differential_workflow() -> None:
    corrected = _build_target(method="ruv_iii_style", k=2)
    repeated = _build_target(method="ruv_iii_style", k=2)

    assert corrected.phospho.shape == (_SITE_COUNT, len(_SAMPLES))
    assert corrected.phospho.index.equals(_site_keys())
    assert corrected.phospho.columns.tolist() == list(_SAMPLES)
    assert corrected.intensity_scale_state.quantity is not None
    assert corrected.intensity_scale_state.quantity.value == "phosphosite_log_abundance"
    assert corrected.intensity_scale_state.quantity.value != "contrast_log2_fold_change"
    assert corrected.processing_state.missing_data.complete_matrix is True
    assert corrected.processing_state.missing_data.imputed is False
    assert corrected.processing_state.missing_data.missing_value_count == 0
    assert corrected.scientifically_equals(repeated)
    assert corrected.provenance == repeated.provenance

    design = ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=condition,
                biological_replicate_id=replicate,
                technical_replicate_id=sample_id,
            )
            for sample_id, condition, replicate in zip(
                _SAMPLES,
                _CONDITIONS,
                _REPLICATES,
                strict=True,
            )
        )
    )
    differential = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=corrected,
            design=design,
            contrasts=(
                Contrast(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A",
                ),
            ),
            config=DifferentialAnalysisConfig(
                technical_replicate_policy=TechnicalReplicatePolicy.MEAN
            ),
        )
    )
    table = differential.table_for("B_vs_A")
    observed = table.set_index("site_key").loc[
        list(_site_keys()[_STABLE_SITE_COUNT:]),
        "logFC",
    ]
    expected = _synthetic_target_experiment().condition_effect_by_site.iloc[
        _STABLE_SITE_COUNT:
    ]
    np.testing.assert_allclose(observed, expected, atol=0.02)

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence
from itertools import permutations

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.errors.input import PhosPyInputError
from phospy.science.batch_correction import (
    RUV_III_ALGORITHM_ID,
    RUV_III_METHOD,
    RuvIIIKernel,
    RuvIIIReplicateStructure,
    run_ruv_iii,
)
from tests.support.site_keys import protein_site_key_index

SAMPLES = ("sample_1", "sample_2", "sample_3", "sample_4", "sample_5", "sample_6")
ASSIGNMENTS = {
    "sample_1": "replicate_a",
    "sample_2": "replicate_a",
    "sample_3": "replicate_b",
    "sample_4": "replicate_b",
    "sample_5": "replicate_c",
    "sample_6": "replicate_c",
}


def _site_keys(count: int = 5) -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=[f"P{position:04d}" for position in range(count)],
        sites=[f"S{position + 1}" for position in range(count)],
    )


def _matrix() -> pd.DataFrame:
    unwanted = np.asarray((-2.0, 2.0, -1.0, 1.0, -3.0, 3.0))
    biological = np.asarray((-4.0, -4.0, 1.0, 1.0, 5.0, 5.0))
    values = np.vstack(
        (
            40.0 + unwanted,
            20.0 - (2.0 * unwanted),
            30.0 + biological + (0.5 * unwanted),
            4.0 + (2.0 * biological),
            np.repeat(7.0, len(SAMPLES)),
        )
    )
    return pd.DataFrame(values, index=_site_keys(), columns=SAMPLES)


def _replicates(
    *,
    sample_order: Sequence[str] = SAMPLES,
    assignments: Mapping[str, str] = ASSIGNMENTS,
) -> RuvIIIReplicateStructure:
    return RuvIIIReplicateStructure.from_assignments(
        sample_order=sample_order,
        replicate_by_sample=assignments,
    )


def _run(
    matrix: pd.DataFrame | None = None,
    *,
    controls: Sequence[str] | None = None,
    replicates: RuvIIIReplicateStructure | None = None,
    k: int = 1,
):
    resolved = _matrix() if matrix is None else matrix
    return run_ruv_iii(
        resolved,
        control_site_keys=(
            tuple(str(value) for value in resolved.index[:2])
            if controls is None
            else controls
        ),
        replicate_structure=_replicates() if replicates is None else replicates,
        k=k,
    )


def test_ruv_iii_k_zero_is_validated_no_op() -> None:
    matrix = _matrix()
    result = _run(matrix, k=0)

    pdt.assert_frame_equal(result.corrected_matrix, matrix)
    assert result.estimated_unwanted_factors.shape == (len(SAMPLES), 0)
    assert result.diagnostics.requested_k == 0
    assert result.diagnostics.effective_k == 0
    assert result.diagnostics.control_loading_rank == 0


@pytest.mark.parametrize("invalid_k", (-1, 1.5, True))
def test_ruv_iii_rejects_invalid_k(invalid_k: object) -> None:
    with pytest.raises(PhosPyInputError, match="RUV-III k must"):
        RuvIIIKernel().run(
            phospho=_matrix(),
            control_site_keys=tuple(str(value) for value in _matrix().index[:2]),
            replicate_structure=_replicates(),
            k=invalid_k,  # type: ignore[arg-type]
        )


def test_ruv_iii_removes_known_unwanted_component_and_reports_ranks() -> None:
    matrix = _matrix()
    result = _run(matrix)
    corrected = result.corrected_matrix

    np.testing.assert_allclose(corrected.iloc[0].to_numpy(), 40.0, atol=1e-12)
    np.testing.assert_allclose(corrected.iloc[1].to_numpy(), 20.0, atol=1e-12)
    np.testing.assert_allclose(
        corrected.iloc[2].to_numpy(),
        np.asarray((26.0, 26.0, 31.0, 31.0, 35.0, 35.0)),
        atol=1e-12,
    )
    assert np.isfinite(corrected.to_numpy()).all()

    diagnostics = result.diagnostics
    assert diagnostics.method == RUV_III_METHOD
    assert diagnostics.algorithm_id == RUV_III_ALGORITHM_ID
    assert diagnostics.sample_count == 6
    assert diagnostics.feature_count == 5
    assert diagnostics.control_count == 2
    assert diagnostics.replicate_set_count == 3
    assert diagnostics.replicate_mapping_rank == 3
    assert diagnostics.replicate_residual_degrees_of_freedom == 3
    assert diagnostics.replicate_residual_rank == 1
    assert diagnostics.control_loading_rank == 1
    assert diagnostics.estimated_unwanted_factor_dimensions == (6, 1)
    assert diagnostics == result.diagnostics
    assert diagnostics.to_payload()["effective_k"] == 1


def test_ruv_iii_retains_signal_not_explained_by_unwanted_variation() -> None:
    matrix = _matrix()
    result = _run(matrix)

    pdt.assert_series_equal(
        result.corrected_matrix.iloc[3],
        matrix.iloc[3],
        check_exact=False,
        atol=1e-12,
        rtol=1e-12,
    )
    original_group_means = matrix.iloc[2].groupby(list(ASSIGNMENTS.values())).mean()
    corrected_group_means = (
        result.corrected_matrix.iloc[2].groupby(list(ASSIGNMENTS.values())).mean()
    )
    np.testing.assert_allclose(corrected_group_means, original_group_means, atol=1e-12)


def test_ruv_iii_is_invariant_to_valid_feature_reordering() -> None:
    matrix = _matrix()
    controls = tuple(str(value) for value in matrix.index[:2])
    expected = _run(matrix, controls=controls).corrected_matrix
    reordered = matrix.iloc[[3, 1, 4, 0, 2], :]

    observed = _run(reordered, controls=controls).corrected_matrix.loc[matrix.index, :]

    pdt.assert_frame_equal(
        observed, expected, check_exact=False, atol=1e-11, rtol=1e-11
    )


def test_ruv_iii_remaps_replicates_under_sample_reordering() -> None:
    matrix = _matrix()
    expected = _run(matrix)
    reordered_samples = (
        "sample_5",
        "sample_2",
        "sample_3",
        "sample_6",
        "sample_1",
        "sample_4",
    )
    reordered = matrix.loc[:, list(reordered_samples)]
    replicates = _replicates(sample_order=reordered_samples)

    observed = _run(reordered, replicates=replicates)

    pdt.assert_frame_equal(
        observed.corrected_matrix.loc[:, list(SAMPLES)],
        expected.corrected_matrix,
        check_exact=False,
        atol=1e-11,
        rtol=1e-11,
    )
    pdt.assert_frame_equal(
        observed.estimated_unwanted_factors.loc[list(SAMPLES), :],
        expected.estimated_unwanted_factors,
        check_exact=False,
        atol=1e-11,
        rtol=1e-11,
    )


def test_ruv_iii_rejects_every_sample_order_when_k_splits_singular_value_tie() -> None:
    samples = ("a", "b", "c", "d")
    assignments = {"a": "pair_1", "b": "pair_1", "c": "pair_2", "d": "pair_2"}
    matrix = pd.DataFrame(
        ((1.0, -1.0, 0.0, 0.0), (0.0, 0.0, 1.0, -1.0)),
        index=_site_keys(2),
        columns=samples,
    )
    controls = tuple(str(value) for value in matrix.index)
    messages: set[str] = set()

    for sample_order in permutations(samples):
        replicates = RuvIIIReplicateStructure.from_assignments(
            sample_order=sample_order,
            replicate_by_sample=assignments,
        )
        with pytest.raises(
            PhosPyInputError,
            match="splits a numerically tied singular-value block",
        ) as error:
            run_ruv_iii(
                matrix.loc[:, list(sample_order)],
                control_site_keys=controls,
                replicate_structure=replicates,
                k=1,
            )
        messages.add(str(error.value))

    assert len(messages) == 1


def test_ruv_iii_remaps_factors_when_k_includes_complete_singular_value_tie() -> None:
    samples = ("a", "b", "c", "d")
    reordered_samples = ("a", "c", "b", "d")
    assignments = {"a": "pair_1", "b": "pair_1", "c": "pair_2", "d": "pair_2"}
    matrix = pd.DataFrame(
        ((1.0, -1.0, 0.0, 0.0), (0.0, 0.0, 1.0, -1.0)),
        index=_site_keys(2),
        columns=samples,
    )
    controls = tuple(str(value) for value in matrix.index)
    expected = run_ruv_iii(
        matrix,
        control_site_keys=controls,
        replicate_structure=RuvIIIReplicateStructure.from_assignments(
            sample_order=samples,
            replicate_by_sample=assignments,
        ),
        k=2,
    )
    observed = run_ruv_iii(
        matrix.loc[:, list(reordered_samples)],
        control_site_keys=controls,
        replicate_structure=RuvIIIReplicateStructure.from_assignments(
            sample_order=reordered_samples,
            replicate_by_sample=assignments,
        ),
        k=2,
    )

    pdt.assert_frame_equal(
        observed.corrected_matrix.loc[:, list(samples)],
        expected.corrected_matrix,
        check_exact=False,
        atol=1e-12,
        rtol=1e-12,
    )
    pdt.assert_frame_equal(
        observed.estimated_unwanted_factors.loc[list(samples), :],
        expected.estimated_unwanted_factors,
        check_exact=False,
        atol=1e-12,
        rtol=1e-12,
    )


def test_ruv_iii_rejects_rank_deficient_replicate_residuals() -> None:
    matrix = _matrix()
    matrix.loc[:, :] = np.asarray(
        (
            (1.0, 1.0, 2.0, 2.0, 3.0, 3.0),
            (2.0, 2.0, 4.0, 4.0, 6.0, 6.0),
            (3.0, 3.0, 4.0, 4.0, 5.0, 5.0),
            (4.0, 4.0, 4.0, 4.0, 4.0, 4.0),
            (5.0, 5.0, 5.0, 5.0, 5.0, 5.0),
        )
    )

    with pytest.raises(PhosPyInputError, match="replicate-residual matrix rank=0"):
        _run(matrix)


def test_ruv_iii_rejects_insufficient_controls_for_k() -> None:
    matrix = _matrix()
    with pytest.raises(PhosPyInputError, match="requires at least 2 negative controls"):
        _run(matrix, controls=(str(matrix.index[0]),), k=2)


def test_ruv_iii_rejects_missing_replicate_assignment() -> None:
    assignments = dict(ASSIGNMENTS)
    assignments.pop("sample_6")

    with pytest.raises(PhosPyInputError, match="missing assignments='sample_6'"):
        _replicates(assignments=assignments)


def test_ruv_iii_rejects_singleton_replicate_sets() -> None:
    assignments = dict(ASSIGNMENTS)
    assignments["sample_6"] = "singleton"

    with pytest.raises(PhosPyInputError, match="singleton set.*replicate_c.*singleton"):
        _replicates(assignments=assignments)


def test_ruv_iii_rejects_inconsistent_sample_order() -> None:
    matrix = _matrix().loc[:, list(reversed(SAMPLES))]

    with pytest.raises(PhosPyInputError, match="sample_order must exactly match"):
        _run(matrix)


@pytest.mark.parametrize("non_finite", (np.nan, np.inf, -np.inf))
def test_ruv_iii_rejects_non_finite_input_without_imputation(non_finite: float) -> None:
    matrix = _matrix()
    matrix.iloc[0, 0] = non_finite

    with pytest.raises(
        PhosPyInputError,
        match="complete finite phosphosite matrix.*does not impute missing values",
    ):
        _run(matrix)


@pytest.mark.parametrize("k", (0, 1))
def test_ruv_iii_rejects_complex_intensities_before_float_conversion(k: int) -> None:
    matrix = _matrix().astype("complex128")
    matrix.iloc[0, 0] += 9j

    with pytest.raises(
        PhosPyInputError,
        match="real-valued phospho intensities.*sample_1",
    ):
        _run(matrix, k=k)


def test_ruv_iii_handles_estimable_high_magnitude_finite_input_stably() -> None:
    samples = ("a", "b", "c", "d")
    assignments = {"a": "pair_1", "b": "pair_1", "c": "pair_2", "d": "pair_2"}
    unwanted = 1e200 * np.asarray((-1.0, 1.0, -2.0, 2.0))
    matrix = pd.DataFrame(
        np.vstack((unwanted, 2.0 * unwanted, 3.0 * unwanted)),
        index=_site_keys(3),
        columns=samples,
    )
    replicates = RuvIIIReplicateStructure.from_assignments(
        sample_order=samples,
        replicate_by_sample=assignments,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = run_ruv_iii(
            matrix,
            control_site_keys=tuple(str(value) for value in matrix.index[:2]),
            replicate_structure=replicates,
            k=1,
        )

    assert np.isfinite(result.corrected_matrix.to_numpy()).all()
    assert np.isfinite(result.estimated_unwanted_factors.to_numpy()).all()


def test_ruv_iii_requires_governed_site_key_identity() -> None:
    matrix = _matrix()
    matrix.index = pd.Index([f"display_{position}" for position in range(5)])

    with pytest.raises(PhosPyInputError, match="index.name='site_key'"):
        _run(matrix)

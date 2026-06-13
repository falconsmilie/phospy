from __future__ import annotations

import pytest

from phospy.science.enrichment import (
    ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    MULTIPLE_TESTING_CORRECTION_NONE,
    EnrichmentSet,
    EnrichmentSetCollection,
)
from phospy.science.enrichment.ora import OraConfig
from phospy.science.enrichment.ora import run as run_ora
from phospy.science.statistics.multiple_testing import (
    run as run_multiple_testing_correction,
)


def _collection(sets: dict[str, tuple[str, ...]]) -> EnrichmentSetCollection:
    return EnrichmentSetCollection(
        sets=tuple(
            EnrichmentSet(
                set_id=set_id,
                name=set_id.replace("_", " "),
                identifiers=identifiers,
                identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            )
            for set_id, identifiers in sets.items()
        )
    )


def test_enrichment_multiple_testing_correction_bh_known_example() -> None:
    adjusted = run_multiple_testing_correction(
        (0.01, 0.04, 0.03, 0.002, 0.05),
        method=MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    )

    assert adjusted == pytest.approx((0.025, 0.05, 0.05, 0.01, 0.05))


def test_enrichment_multiple_testing_correction_none_preserves_raw_values() -> None:
    raw = (0.2, 0.01, 1.0)

    adjusted = run_multiple_testing_correction(
        raw,
        method=MULTIPLE_TESTING_CORRECTION_NONE,
    )

    assert adjusted == pytest.approx(raw)


def test_enrichment_multiple_testing_correction_preserves_input_order() -> None:
    adjusted = run_multiple_testing_correction(
        (0.04, 0.001, 0.03),
        method=MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    )

    assert adjusted == pytest.approx((0.04, 0.003, 0.04))


def test_enrichment_multiple_testing_correction_empty_result_list() -> None:
    adjusted = run_multiple_testing_correction(
        (),
        method=MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    )

    assert adjusted == ()


def test_enrichment_multiple_testing_correction_all_one_p_values() -> None:
    adjusted = run_multiple_testing_correction(
        (1.0, 1.0, 1.0),
        method=MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    )

    assert adjusted == pytest.approx((1.0, 1.0, 1.0))


def test_enrichment_multiple_testing_correction_zero_p_values() -> None:
    adjusted = run_multiple_testing_correction(
        (0.0, 0.01, 0.0),
        method=MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    )

    assert adjusted == pytest.approx((0.0, 0.01, 0.0))


def test_enrichment_ora_correction_records_bh_adjusted_p_values() -> None:
    result = run_ora(
        selected_identifiers=("A", "B", "C"),
        background_universe=("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
        enrichment_sets=_collection(
            {
                "SET_A": ("A", "B", "C"),
                "SET_B": ("A", "B", "D", "E"),
                "SET_C": ("A", "D", "E", "F"),
            }
        ),
    )

    raw_p_values = tuple(record.p_value for record in result.records)
    expected_adjusted = run_multiple_testing_correction(
        raw_p_values,
        method=MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    )

    assert tuple(record.adjusted_p_value for record in result.records) == (
        pytest.approx(expected_adjusted)
    )
    assert all(
        record.correction_method == MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG
        for record in result.records
    )
    assert tuple(record.p_value for record in result.records) == raw_p_values


def test_enrichment_ora_correction_none_leaves_adjusted_equal_to_raw() -> None:
    result = run_ora(
        selected_identifiers=("A", "B", "C"),
        background_universe=("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
        enrichment_sets=_collection(
            {
                "SET_A": ("A", "B", "C"),
                "SET_B": ("A", "B", "D", "E"),
            }
        ),
        config=OraConfig(multiple_testing_correction=MULTIPLE_TESTING_CORRECTION_NONE),
    )

    assert tuple(record.adjusted_p_value for record in result.records) == pytest.approx(
        tuple(record.p_value for record in result.records)
    )
    assert all(
        record.correction_method == MULTIPLE_TESTING_CORRECTION_NONE
        for record in result.records
    )


def test_enrichment_ora_correction_keeps_existing_stable_row_ordering() -> None:
    result = run_ora(
        selected_identifiers=("A", "B"),
        background_universe=("A", "B", "C", "D", "E", "F"),
        enrichment_sets=_collection(
            {
                "B_TIE": ("A", "C"),
                "Z_BEST": ("A", "B"),
                "A_TIE": ("B", "D"),
            }
        ),
    )

    assert tuple(row.set_id for row in result.records) == (
        "Z_BEST",
        "A_TIE",
        "B_TIE",
    )
    assert all(row.adjusted_p_value is not None for row in result.records)

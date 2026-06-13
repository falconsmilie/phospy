from __future__ import annotations

import pytest

from phospy.science.enrichment import (
    ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    EnrichmentSet,
    EnrichmentSetCollection,
)
from phospy.science.enrichment.ora import (
    ORA_OUTSIDE_BACKGROUND_POLICY_DROP,
    ORA_OUTSIDE_BACKGROUND_POLICY_ERROR,
    OraConfig,
    OraEngine,
    run,
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


def test_enrichment_ora_simple_hypergeometric_calculation_known_p_value() -> None:
    result = run(
        selected_identifiers=("A", "B", "C"),
        background_universe=("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
        enrichment_sets=_collection({"SET": ("A", "B", "D", "E")}),
    )

    row = result.records[0]
    assert row.background_size == 10
    assert row.selected_size == 3
    assert row.set_size == 4
    assert row.overlap_size == 2
    assert row.overlap_identifiers == ("A", "B")
    assert row.p_value == pytest.approx(1.0 / 3.0)
    assert row.enrichment_ratio == pytest.approx(2.0 / (3.0 * 4.0 / 10.0))


def test_enrichment_ora_background_universe_affects_p_value() -> None:
    selected = ("A", "B", "C")
    collection = _collection({"SET": ("A", "B", "D", "E")})

    small_background = run(
        selected_identifiers=selected,
        background_universe=("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
        enrichment_sets=collection,
    )
    large_background = run(
        selected_identifiers=selected,
        background_universe=(
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "J",
            "K",
            "L",
            "M",
            "N",
            "O",
            "P",
            "Q",
            "R",
            "S",
            "T",
        ),
        enrichment_sets=collection,
    )

    assert small_background.records[0].p_value == pytest.approx(1.0 / 3.0)
    assert large_background.records[0].p_value == pytest.approx(100.0 / 1140.0)
    assert large_background.records[0].p_value < small_background.records[0].p_value


def test_enrichment_ora_empty_selected_set_returns_non_significant_rows() -> None:
    result = run(
        selected_identifiers=(),
        background_universe=("A", "B", "C", "D"),
        enrichment_sets=_collection({"SET": ("A", "B")}),
    )

    row = result.records[0]
    assert result.selected_size == 0
    assert row.selected_size == 0
    assert row.overlap_size == 0
    assert row.p_value == pytest.approx(1.0)
    assert row.enrichment_ratio is None


def test_enrichment_ora_empty_enrichment_set_after_background_filtering() -> None:
    result = run(
        selected_identifiers=("A", "B"),
        background_universe=("A", "B", "C", "D"),
        enrichment_sets=_collection({"EMPTY_AFTER_BACKGROUND": ("X", "Y")}),
    )

    row = result.records[0]
    assert row.raw_set_size == 2
    assert row.set_size == 0
    assert row.set_identifiers_outside_background_count == 2
    assert row.overlap_size == 0
    assert row.p_value == pytest.approx(1.0)
    assert row.enrichment_ratio is None


def test_enrichment_ora_no_overlap_returns_p_value_one_and_zero_ratio() -> None:
    result = run(
        selected_identifiers=("A", "B"),
        background_universe=("A", "B", "C", "D"),
        enrichment_sets=_collection({"NO_OVERLAP": ("C", "D")}),
    )

    row = result.records[0]
    assert row.overlap_size == 0
    assert row.overlap_identifiers == ()
    assert row.p_value == pytest.approx(1.0)
    assert row.enrichment_ratio == pytest.approx(0.0)


def test_enrichment_ora_identifiers_outside_background_follow_policy() -> None:
    with pytest.raises(ValueError, match="selected_identifiers"):
        run(
            selected_identifiers=("A", "Z"),
            background_universe=("A", "B", "C"),
            enrichment_sets=_collection({"SET": ("A", "B")}),
        )

    result = run(
        selected_identifiers=("A", "Z"),
        background_universe=("A", "B", "C"),
        enrichment_sets=_collection({"SET": ("A", "B")}),
        config=OraConfig(
            selected_outside_background_policy=ORA_OUTSIDE_BACKGROUND_POLICY_DROP
        ),
    )

    assert result.selected_identifiers == ("A",)
    assert result.dropped_selected_identifiers == ("Z",)
    assert result.records[0].selected_size == 1

    with pytest.raises(ValueError, match="set_id='SET'"):
        run(
            selected_identifiers=("A",),
            background_universe=("A", "B", "C"),
            enrichment_sets=_collection({"SET": ("A", "X")}),
            config=OraConfig(
                set_outside_background_policy=ORA_OUTSIDE_BACKGROUND_POLICY_ERROR
            ),
        )


def test_enrichment_ora_deterministic_row_ordering() -> None:
    result = OraEngine().run(
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

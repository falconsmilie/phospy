from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.protein_aware_alignment import (
    PROTEIN_AWARE_REASON_AMBIGUOUS_PROTEIN_MAPPING,
    PROTEIN_AWARE_REASON_INCOMPATIBLE_TRANSFORMATION_STATE,
    PROTEIN_AWARE_REASON_MATCHED_PROTEIN_AVAILABLE,
    PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW,
    PROTEIN_AWARE_REASON_SAMPLE_MISMATCH,
    ProteinAwareAlignmentConfig,
    ProteinAwareAlignmentEligibilityResolver,
    ProteinAwarePreparationEligibility,
)
from phospy.science.datasets.preprocessing.protein_mapping import (
    ProteinMappingConfig,
    ProteinMappingResolver,
)
from phospy.science.transformations.models import (
    IntensityScaleState,
    MatrixIntensityScaleState,
)


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [3.0, 4.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )


def _site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "protein_accession": ["P53778", "P31749"],
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
        },
        index=_phospho().index.copy(),
    )


def _total(
    *,
    columns: list[str] | None = None,
    index: list[str] | None = None,
) -> pd.DataFrame:
    resolved_columns = ["sample_a", "sample_b"] if columns is None else columns
    resolved_index = ["P53778", "P31749"] if index is None else index
    return pd.DataFrame(
        {
            column: [float(row + position + 1) for row in range(len(resolved_index))]
            for position, column in enumerate(resolved_columns)
        },
        index=pd.Index(resolved_index, name="protein_id"),
    )


def _compatible_state() -> IntensityScaleState:
    return IntensityScaleState(
        phospho=MatrixIntensityScaleState.log2(established_by="test.phospho"),
        total=MatrixIntensityScaleState.log2(established_by="test.total"),
    )


def _incompatible_state() -> IntensityScaleState:
    return IntensityScaleState(
        phospho=MatrixIntensityScaleState.log2(established_by="test.phospho"),
        total=MatrixIntensityScaleState.linear(established_by="test.total"),
    )


def _mapping_result(
    *,
    site_metadata: pd.DataFrame | None = None,
    phospho: pd.DataFrame | None = None,
    total: pd.DataFrame | None = None,
):
    phospho_matrix = _phospho() if phospho is None else phospho
    total_matrix = _total() if total is None else total
    metadata = _site_metadata() if site_metadata is None else site_metadata
    return ProteinMappingResolver().run(
        site_metadata=metadata,
        phospho_matrix_index=phospho_matrix.index,
        total_protein_matrix_index=total_matrix.index,
        config=ProteinMappingConfig(protein_identifier_columns=("protein_accession",)),
    )


def _diagnose(
    *,
    phospho: pd.DataFrame | None = None,
    total: pd.DataFrame | None = None,
    site_metadata: pd.DataFrame | None = None,
    intensity_scale_state: IntensityScaleState | None = None,
    config: ProteinAwareAlignmentConfig | None = None,
):
    phospho_matrix = _phospho() if phospho is None else phospho
    total_matrix = _total() if total is None else total
    mapping = _mapping_result(
        site_metadata=site_metadata,
        phospho=phospho_matrix,
        total=total_matrix,
    )
    return ProteinAwareAlignmentEligibilityResolver().run(
        phospho=phospho_matrix,
        total=total_matrix,
        mapping_result=mapping,
        intensity_scale_state=(
            _compatible_state()
            if intensity_scale_state is None
            else intensity_scale_state
        ),
        config=config,
    )


def test_protein_aware_sample_alignment_matching_sample_columns() -> None:
    diagnostics = _diagnose()

    assert diagnostics.sample_alignment.phospho_sample_columns == (
        "sample_a",
        "sample_b",
    )
    assert diagnostics.sample_alignment.total_protein_sample_columns == (
        "sample_a",
        "sample_b",
    )
    assert diagnostics.sample_alignment.exact_sample_order_match is True
    assert diagnostics.sample_alignment.sample_order_compatible is True
    assert diagnostics.sample_alignment.missing_total_protein_samples == ()
    assert diagnostics.sample_alignment.extra_total_protein_samples == ()
    assert diagnostics.global_reasons == ()


def test_protein_aware_sample_alignment_reordered_columns_when_explicitly_supported() -> (
    None
):
    diagnostics = _diagnose(
        total=_total(columns=["sample_b", "sample_a"]),
        config=ProteinAwareAlignmentConfig(allow_reordered_samples=True),
    )

    assert diagnostics.sample_alignment.exact_sample_order_match is False
    assert diagnostics.sample_alignment.reordered_sample_columns is True
    assert diagnostics.sample_alignment.sample_order_compatible is True
    assert diagnostics.sample_alignment.missing_total_protein_samples == ()
    assert diagnostics.sample_alignment.extra_total_protein_samples == ()
    assert diagnostics.global_reasons == ()


def test_protein_aware_sample_alignment_missing_total_sample_diagnostic() -> None:
    diagnostics = _diagnose(total=_total(columns=["sample_a"]))

    assert diagnostics.sample_alignment.sample_order_compatible is False
    assert diagnostics.sample_alignment.missing_total_protein_samples == ("sample_b",)
    assert diagnostics.sample_alignment.extra_total_protein_samples == ()
    assert PROTEIN_AWARE_REASON_SAMPLE_MISMATCH in diagnostics.global_reasons
    assert diagnostics.excluded_from_preparation == (
        "MAPK14;Y182;",
        "AKT1;T308;",
    )


def test_protein_aware_sample_alignment_extra_total_sample_diagnostic() -> None:
    diagnostics = _diagnose(total=_total(columns=["sample_a", "sample_b", "sample_c"]))

    assert diagnostics.sample_alignment.sample_order_compatible is False
    assert diagnostics.sample_alignment.missing_total_protein_samples == ()
    assert diagnostics.sample_alignment.extra_total_protein_samples == ("sample_c",)
    assert PROTEIN_AWARE_REASON_SAMPLE_MISMATCH in diagnostics.global_reasons


def test_protein_aware_sample_alignment_incompatible_transformation_state() -> None:
    diagnostics = _diagnose(intensity_scale_state=_incompatible_state())

    assert diagnostics.transformation_state.compatible is False
    assert diagnostics.transformation_state.phospho_transformation_state == {
        "kind": "log2",
        "transformed": True,
        "established_by": "test.phospho",
    }
    assert diagnostics.transformation_state.total_protein_transformation_state == {
        "kind": "linear",
        "transformed": False,
        "established_by": "test.total",
    }
    assert (
        PROTEIN_AWARE_REASON_INCOMPATIBLE_TRANSFORMATION_STATE
        in diagnostics.global_reasons
    )
    assert diagnostics.eligibility_by_site["MAPK14;Y182;"] is (
        ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION
    )


def test_protein_aware_sample_alignment_eligible_site() -> None:
    diagnostics = _diagnose()

    assert diagnostics.eligible_for_protein_aware_preparation == (
        "MAPK14;Y182;",
        "AKT1;T308;",
    )
    assert diagnostics.fallback_to_phospho_only == ()
    assert diagnostics.excluded_from_preparation == ()
    assert diagnostics.reasons_by_site["MAPK14;Y182;"] == (
        PROTEIN_AWARE_REASON_MATCHED_PROTEIN_AVAILABLE,
    )


def test_protein_aware_sample_alignment_fallback_site_due_to_missing_protein_row() -> (
    None
):
    diagnostics = _diagnose(
        total=_total(index=["P53778"]),
        config=ProteinAwareAlignmentConfig(
            protein_mapping_policy="allow_missing_with_report"
        ),
    )

    assert diagnostics.eligible_for_protein_aware_preparation == ("MAPK14;Y182;",)
    assert diagnostics.fallback_to_phospho_only == ("AKT1;T308;",)
    assert diagnostics.excluded_from_preparation == ()
    assert diagnostics.reasons_by_site["AKT1;T308;"] == (
        PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW,
    )


def test_protein_aware_sample_alignment_require_unambiguous_excludes_missing_row() -> (
    None
):
    diagnostics = _diagnose(total=_total(index=["P53778"]))

    assert diagnostics.eligible_for_protein_aware_preparation == ("MAPK14;Y182;",)
    assert diagnostics.fallback_to_phospho_only == ()
    assert diagnostics.excluded_from_preparation == ("AKT1;T308;",)
    assert diagnostics.reasons_by_site["AKT1;T308;"] == (
        PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW,
    )


def test_protein_aware_sample_alignment_excluded_site_due_to_ambiguity() -> None:
    site_metadata = pd.DataFrame(
        {
            "protein_accession": ["P53778", "Q16539", "P31749"],
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1"],
            "site": ["Y182", "Y182", "T308"],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "MAPK14;Y182;", "AKT1;T308;"],
            name="site_id",
        ),
    )

    diagnostics = _diagnose(
        total=_total(index=["P53778", "Q16539", "P31749"]),
        site_metadata=site_metadata,
        config=ProteinAwareAlignmentConfig(
            protein_mapping_policy="allow_missing_with_report"
        ),
    )

    assert diagnostics.eligible_for_protein_aware_preparation == ("AKT1;T308;",)
    assert diagnostics.fallback_to_phospho_only == ()
    assert diagnostics.excluded_from_preparation == ("MAPK14;Y182;",)
    assert diagnostics.reasons_by_site["MAPK14;Y182;"] == (
        PROTEIN_AWARE_REASON_AMBIGUOUS_PROTEIN_MAPPING,
    )


def test_protein_aware_sample_alignment_config_rejects_non_bool_reorder_policy() -> (
    None
):
    with pytest.raises(PhosPyInputError, match="allow_reordered_samples"):
        ProteinAwareAlignmentConfig(
            allow_reordered_samples="yes"  # type: ignore[arg-type]
        )

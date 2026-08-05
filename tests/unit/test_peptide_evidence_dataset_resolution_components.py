from __future__ import annotations

from collections.abc import Mapping

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api import PhosPyInputError
from phospy.science.evidence.dataset_resolution import (
    DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
    DATASET_MULTI_SITE_POLICY_SPLIT,
    DATASET_PEPTIDE_ALLOCATION_DOMAIN_DECLARED_SCALE_UNIT_MAPPING_PASSTHROUGH,
    DATASET_PEPTIDE_INPUT_QUANTITATIVE_MEANING_PEPTIDE_LOG2_ABUNDANCE,
    DATASET_PEPTIDE_LOCALISATION_SUMMARY_COLUMN,
    DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS,
    DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS_COLUMN,
    DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL,
    DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT,
    DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL,
    DATASET_PEPTIDE_MISSING_VALUE_POLICY_FINITE_MEAN,
    DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION,
    DATASET_PEPTIDE_SIGNAL_CONSERVATION_POLICY_NOT_CONSERVED,
    DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR,
    DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS,
    DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS,
    DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LINEAR_ALLOCATED_MEAN_V1,
    DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN,
    PeptideEvidenceDatasetResolver,
    PeptideEvidenceResolutionResult,
    PeptideEvidenceResolutionSummary,
    build_multi_site_handling_config_for_dataset_policy,
    normalise_peptide_evidence_resolution_summary_payload,
)
from phospy.science.evidence.dataset_resolution.allocation import (
    AllocatedEvidence,
    SiteSignalSummary,
    allocate_peptide_signals_to_resolved_sites,
    summarise_allocated_site_signals,
)
from phospy.science.evidence.dataset_resolution.contracts import (
    MAPPING_FRACTION_COLUMN,
    PeptideEvidenceResolutionInputMetrics,
    build_peptide_to_site_aggregation_policy,
)
from phospy.science.evidence.dataset_resolution.mapping import (
    JoinedMappingRows,
    ResolvedMappingFractions,
    join_peptide_rows_to_site_mapping,
    resolve_and_validate_mapping_fractions,
)
from phospy.science.evidence.dataset_resolution.site_metadata import (
    SiteMetadataResolution,
    aggregate_site_metadata_and_localisation,
)
from phospy.science.evidence.dataset_resolution.site_sequence import (
    SITE_SEQUENCE_SOURCE_MISSING,
    SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT,
    SITE_SEQUENCE_SOURCE_PROVIDED,
    SiteSequenceResolutionDiagnostics,
    derive_site_sequence_from_peptide_context,
    resolve_site_sequence_for_resolved_site,
)
from phospy.science.evidence.dataset_resolution.summary import (
    build_resolution_summary,
)
from phospy.science.transformations.models import IntensityScaleKind


def _evidence_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "peptide_row_id": "pep_split",
                "protein_accession": "P28482",
                "site_string": "S10,T12",
                "peptide_sequence": "AAASTTAAAA",
                "multi_site": True,
                "site_sequence": "AAASAAA",
                "localisation_confidence": 0.9,
                "sample_a": 10.0,
                "sample_b": 20.0,
            },
            {
                "peptide_row_id": "pep_single",
                "protein_accession": "P31749",
                "site_string": "S473",
                "peptide_sequence": "GGGSGGG",
                "multi_site": False,
                "site_sequence": "GGGSGGG",
                "localisation_confidence": 0.8,
                "sample_a": 7.0,
                "sample_b": 9.0,
            },
        ]
    )


def test_import_path_preserves_dataset_resolution_contract_exports() -> None:
    assert PeptideEvidenceDatasetResolver is not None
    assert PeptideEvidenceResolutionResult is not None
    assert PeptideEvidenceResolutionSummary is not None
    assert callable(build_multi_site_handling_config_for_dataset_policy)


def test_mapping_uses_explicit_mapping_weights_and_preserves_mapped_row_order() -> None:
    mapping = pd.DataFrame(
        {
            "peptide_row_id": ["pep_single", "pep_split", "pep_split"],
            "site_id": ["AKT1;S473;", "MAPK1;T12;", "MAPK1;S10;"],
            "mapping_weight": [1.0, 0.25, 0.75],
        }
    )

    joined = join_peptide_rows_to_site_mapping(
        evidence_frame=_evidence_frame(),
        mapping=mapping,
        sample_columns=("sample_a", "sample_b"),
    )
    resolved = resolve_and_validate_mapping_fractions(joined_mapping=joined)

    assert joined.rows.loc[:, "site_id"].tolist() == mapping.loc[:, "site_id"].tolist()
    assert (
        resolved.mapping_weight_source == DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT
    )
    assert resolved.rows.loc[:, MAPPING_FRACTION_COLUMN].tolist() == [
        1.0,
        0.25,
        0.75,
    ]


def test_mapping_derives_equal_mapping_weights_when_explicit_weights_are_absent() -> (
    None
):
    joined = JoinedMappingRows(
        rows=pd.DataFrame(
            {
                "peptide_row_id": ["pep_split", "pep_split", "pep_single"],
                "site_id": ["MAPK1;S10;", "MAPK1;T12;", "AKT1;S473;"],
            }
        )
    )

    resolved = resolve_and_validate_mapping_fractions(joined_mapping=joined)

    assert (
        resolved.mapping_weight_source
        == DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL
    )
    assert resolved.rows.loc[:, MAPPING_FRACTION_COLUMN].tolist() == [0.5, 0.5, 1.0]


@pytest.mark.parametrize(
    ("weights", "message"),
    [
        (("bad", 1.0), "non-positive or non-numeric"),
        ((0.0, 1.0), "non-positive or non-numeric"),
        ((0.7, 0.4), "must sum to 1.0 per peptide_row_id"),
    ],
)
def test_mapping_rejects_invalid_mapping_weights(
    weights: tuple[object, object],
    message: str,
) -> None:
    joined = JoinedMappingRows(
        rows=pd.DataFrame(
            {
                "peptide_row_id": ["pep_split", "pep_split"],
                "site_id": ["MAPK1;S10;", "MAPK1;T12;"],
                "mapping_weight": list(weights),
            }
        )
    )

    with pytest.raises(PhosPyInputError, match=message):
        resolve_and_validate_mapping_fractions(joined_mapping=joined)


def test_allocation_applies_weights_and_site_summarisation_preserves_matrix_contract() -> (
    None
):
    resolved_mapping = ResolvedMappingFractions(
        rows=pd.DataFrame(
            {
                "peptide_row_id": [
                    "pep_split",
                    "pep_split",
                    "pep_duplicate",
                    "pep_missing",
                ],
                "site_id": [
                    "MAPK1;T12;",
                    "MAPK1;S10;",
                    "MAPK1;S10;",
                    "AKT1;S473;",
                ],
                "peptide_sequence": [
                    "AAASTTAAAA",
                    "AAASTTAAAA",
                    "AAASTTAAAA",
                    "GGGSGGG",
                ],
                MAPPING_FRACTION_COLUMN: [0.25, 0.75, 1.0, 1.0],
                "sample_b": [20.0, 20.0, 24.0, 40.0],
                "sample_a": [10.0, 10.0, 14.0, float("nan")],
            }
        ),
        mapping_weight_source=DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT,
    )

    allocated = allocate_peptide_signals_to_resolved_sites(
        resolved_mapping=resolved_mapping,
        sample_columns=("sample_b", "sample_a"),
        aggregation_policy=build_peptide_to_site_aggregation_policy(
            input_intensity_scale=IntensityScaleKind.LINEAR,
            mapping_rows=resolved_mapping.rows,
        ),
    )
    site_signals = summarise_allocated_site_signals(
        allocated_evidence=allocated,
    )

    assert allocated.rows.loc[0, "sample_b"] == pytest.approx(5.0)
    assert allocated.rows.loc[1, "sample_a"] == pytest.approx(7.5)
    assert pd.isna(allocated.rows.loc[3, "sample_a"])
    assert site_signals.phospho.columns.tolist() == ["sample_b", "sample_a"]
    assert site_signals.phospho.index.name == "site_id"
    assert site_signals.phospho.index.astype(str).tolist() == [
        "AKT1;S473;",
        "MAPK1;S10;",
        "MAPK1;T12;",
    ]
    assert float(site_signals.phospho.loc["MAPK1;S10;", "sample_b"]) == pytest.approx(
        (15.0 + 24.0) / 2.0
    )
    assert float(site_signals.phospho.loc["MAPK1;S10;", "sample_a"]) == pytest.approx(
        (7.5 + 14.0) / 2.0
    )
    assert pd.isna(site_signals.phospho.loc["AKT1;S473;", "sample_a"])


def test_allocation_requires_typed_peptide_to_site_policy() -> None:
    resolved_mapping = ResolvedMappingFractions(
        rows=pd.DataFrame(
            {
                "peptide_row_id": ["pep_single"],
                "site_id": ["AKT1;S473;"],
                MAPPING_FRACTION_COLUMN: [1.0],
                "sample_a": [7.0],
            }
        ),
        mapping_weight_source=DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT,
    )

    with pytest.raises(PhosPyInputError, match="typed PeptideToSiteAggregationPolicy"):
        allocate_peptide_signals_to_resolved_sites(
            resolved_mapping=resolved_mapping,
            sample_columns=("sample_a",),
            aggregation_policy=object(),  # type: ignore[arg-type]
        )


def test_log2_unit_mapping_policy_is_typed_passthrough_not_fractional_inversion() -> (
    None
):
    resolved_mapping = ResolvedMappingFractions(
        rows=pd.DataFrame(
            {
                "peptide_row_id": ["pep_single"],
                "site_id": ["AKT1;S473;"],
                MAPPING_FRACTION_COLUMN: [1.0],
                "sample_a": [3.0],
            }
        ),
        mapping_weight_source=DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT,
    )
    policy = build_peptide_to_site_aggregation_policy(
        input_intensity_scale=IntensityScaleKind.LOG2,
        mapping_rows=resolved_mapping.rows,
    )

    allocated = allocate_peptide_signals_to_resolved_sites(
        resolved_mapping=resolved_mapping,
        sample_columns=("sample_a",),
        aggregation_policy=policy,
    )

    assert allocated.rows.loc[0, "sample_a"] == pytest.approx(3.0)
    payload = policy.to_payload()
    assert payload["input_quantitative_meaning"] == (
        DATASET_PEPTIDE_INPUT_QUANTITATIVE_MEANING_PEPTIDE_LOG2_ABUNDANCE
    )
    assert payload["allocation_domain"] == (
        DATASET_PEPTIDE_ALLOCATION_DOMAIN_DECLARED_SCALE_UNIT_MAPPING_PASSTHROUGH
    )
    assert payload["output_intensity_scale"] == "log2"
    assert payload["output_quantitative_meaning"] == "phosphosite_log_abundance"
    assert payload["fractional_mapping_present"] is False
    assert "2**" not in str(payload["aggregation_formula"])


def test_metadata_aggregates_one_site_repeated_evidence_and_localisation() -> None:
    allocated_evidence = AllocatedEvidence(
        rows=pd.DataFrame(
            {
                "site_id": ["MAPK1;S10;", "MAPK1;S10;", "MAPK1;S10;"],
                "protein_accession": ["P28482", "P28482", "P28482"],
                "site_string": ["S10", "S10", "S10"],
                "peptide_sequence": ["AAASAAA", "CCCSCCC", "GGGSGGG"],
                "multi_site": [False, False, False],
                "site_sequence": ["AAASAAA", "AAASAAA", "AAASAAA"],
                "localisation_confidence": [0.2, pd.NA, 1.0],
                "sample_a": [10.0, 20.0, 30.0],
            }
        ),
        sample_columns=("sample_a",),
    )

    resolved = aggregate_site_metadata_and_localisation(
        allocated_evidence=allocated_evidence,
        site_ids=pd.Index(["MAPK1;S10;"], name="site_id"),
        multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
    )

    assert resolved.site_metadata.index.name == "site_id"
    assert resolved.site_metadata.loc["MAPK1;S10;", "gene_symbol"] == "MAPK1"
    assert resolved.site_metadata.loc["MAPK1;S10;", "site"] == "S10"
    assert resolved.site_metadata.loc["MAPK1;S10;", "protein_accession"] == "P28482"
    assert float(
        resolved.site_metadata.loc["MAPK1;S10;", "localisation_confidence"]
    ) == pytest.approx(0.6)
    assert float(
        resolved.site_metadata.loc[
            "MAPK1;S10;",
            DATASET_PEPTIDE_LOCALISATION_SUMMARY_COLUMN,
        ]
    ) == pytest.approx(0.6)
    assert (
        resolved.site_metadata.loc[
            "MAPK1;S10;",
            DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS_COLUMN,
        ]
        == DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS
    )
    assert resolved.sequence_diagnostics.provided_site_sequence_used_count == 1


def test_metadata_records_missing_localisation_when_no_finite_values_exist() -> None:
    allocated_evidence = AllocatedEvidence(
        rows=pd.DataFrame(
            {
                "site_id": ["AKT1;S473;"],
                "protein_accession": ["P31749"],
                "site_string": ["S473"],
                "peptide_sequence": ["GGGSGGG"],
                "multi_site": [False],
                "site_sequence": ["GGGSGGG"],
                "localisation_confidence": [pd.NA],
                "sample_a": [7.0],
            }
        ),
        sample_columns=("sample_a",),
    )

    resolved = aggregate_site_metadata_and_localisation(
        allocated_evidence=allocated_evidence,
        site_ids=pd.Index(["AKT1;S473;"], name="site_id"),
        multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
    )

    assert pd.isna(resolved.site_metadata.loc["AKT1;S473;", "localisation_confidence"])


def test_metadata_rejects_conflicting_metadata_for_resolved_site() -> None:
    allocated_evidence = AllocatedEvidence(
        rows=pd.DataFrame(
            {
                "site_id": ["MAPK1;S10;", "MAPK1;S10;"],
                "protein_accession": ["P28482", "Q9Y243"],
                "site_string": ["S10", "S10"],
                "peptide_sequence": ["AAASAAA", "CCCSCCC"],
                "multi_site": [False, False],
                "site_sequence": ["AAASAAA", "AAASAAA"],
                "sample_a": [10.0, 20.0],
            }
        ),
        sample_columns=("sample_a",),
    )

    with pytest.raises(PhosPyInputError) as exc_info:
        aggregate_site_metadata_and_localisation(
            allocated_evidence=allocated_evidence,
            site_ids=pd.Index(["MAPK1;S10;"], name="site_id"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
        )

    message = str(exc_info.value)
    assert "protein_accession" in message
    assert "P28482" in message
    assert "Q9Y243" in message


def test_metadata_derives_sequences_for_split_multi_site_rows() -> None:
    allocated_evidence = AllocatedEvidence(
        rows=pd.DataFrame(
            {
                "site_id": ["MAPK1;S10;", "MAPK1;T12;"],
                "protein_accession": ["P28482", "P28482"],
                "site_string": ["S10,T12", "S10,T12"],
                "peptide_sequence": ["AAASTTAAAA", "AAASTTAAAA"],
                "multi_site": [True, True],
                "sample_a": [5.0, 5.0],
            }
        ),
        sample_columns=("sample_a",),
    )

    resolved = aggregate_site_metadata_and_localisation(
        allocated_evidence=allocated_evidence,
        site_ids=pd.Index(["MAPK1;S10;", "MAPK1;T12;"], name="site_id"),
        multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
    )

    assert resolved.site_metadata.loc["MAPK1;S10;", "site_sequence"] == "AAASTTA"
    assert resolved.site_metadata.loc["MAPK1;T12;", "site_sequence"] == "AASTTAAAA"
    assert (
        resolved.sequence_diagnostics.peptide_context_derived_site_sequence_count == 2
    )


def test_site_sequence_accepts_valid_supplied_sequence() -> None:
    resolved = resolve_site_sequence_for_resolved_site(
        group=pd.DataFrame(
            {
                "site_sequence": [" aaAsaaa "],
                "site_string": ["S10"],
                "peptide_sequence": ["AAASAAA"],
                "multi_site": [False],
            }
        ),
        site_id="MAPK1;S10;",
        resolved_site_token="S10",
        multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
    )

    assert resolved.site_sequence == "AAASAAA"
    assert resolved.source == SITE_SEQUENCE_SOURCE_PROVIDED
    assert resolved.rejected_provided_context_count == 0


def test_site_sequence_rejects_invalid_supplied_sequence() -> None:
    with pytest.raises(PhosPyInputError) as exc_info:
        resolve_site_sequence_for_resolved_site(
            group=pd.DataFrame(
                {
                    "site_sequence": ["AAATAAA"],
                    "site_string": ["S10"],
                    "peptide_sequence": ["AAASAAA"],
                    "multi_site": [False],
                }
            ),
            site_id="MAPK1;S10;",
            resolved_site_token="S10",
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
        )

    message = str(exc_info.value)
    assert "invalid_supplied_values" in message
    assert "AAATAAA" in message
    assert "expected='S'" in message
    assert "observed='T'" in message


def test_site_sequence_rejects_conflicting_supplied_sequences_with_stable_preview() -> (
    None
):
    with pytest.raises(PhosPyInputError) as exc_info:
        resolve_site_sequence_for_resolved_site(
            group=pd.DataFrame(
                {
                    "site_sequence": ["CCCSCCC", "AAASAAA"],
                    "site_string": ["S10", "S10"],
                    "peptide_sequence": ["CCCSCCC", "AAASAAA"],
                    "multi_site": [False, False],
                }
            ),
            site_id="MAPK1;S10;",
            resolved_site_token="S10",
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
        )

    message = str(exc_info.value)
    assert "distinct_normalized_value_count=2" in message
    assert "values=['AAASAAA', 'CCCSCCC']" in message
    assert "row order" in message


def test_site_sequence_rejects_mixed_valid_and_invalid_supplied_evidence() -> None:
    with pytest.raises(PhosPyInputError) as exc_info:
        resolve_site_sequence_for_resolved_site(
            group=pd.DataFrame(
                {
                    "site_sequence": ["AAASAAA", "AAATAAA"],
                    "site_string": ["S10", "S10"],
                    "peptide_sequence": ["AAASAAA", "AAASAAA"],
                    "multi_site": [False, False],
                }
            ),
            site_id="MAPK1;S10;",
            resolved_site_token="S10",
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
        )

    message = str(exc_info.value)
    assert "Mixed valid and invalid supplied evidence" in message
    assert "valid_normalized_value='AAASAAA'" in message


def test_site_sequence_derives_from_peptide_context_and_split_multisite_context() -> (
    None
):
    group = pd.DataFrame(
        {
            "site_string": ["S10,T12"],
            "peptide_sequence": ["AAASTTAAAA"],
            "multi_site": [True],
        }
    )

    derivation = derive_site_sequence_from_peptide_context(
        group=group,
        site_id="MAPK1;S10;",
        resolved_site_token="S10",
    )
    resolved = resolve_site_sequence_for_resolved_site(
        group=group,
        site_id="MAPK1;S10;",
        resolved_site_token="S10",
        multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
    )

    assert derivation.site_sequence == "AAASTTA"
    assert derivation.distinct_sequences == ("AAASTTA",)
    assert resolved.site_sequence == "AAASTTA"
    assert resolved.source == SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT


def test_site_sequence_falls_back_to_peptide_context_for_invalid_split_evidence() -> (
    None
):
    resolved = resolve_site_sequence_for_resolved_site(
        group=pd.DataFrame(
            {
                "site_sequence": ["AAATAAA"],
                "site_string": ["S10,T12"],
                "peptide_sequence": ["AAASTTAAAA"],
                "multi_site": [True],
            }
        ),
        site_id="MAPK1;S10;",
        resolved_site_token="S10",
        multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
    )

    assert resolved.site_sequence == "AAASTTA"
    assert resolved.source == SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT
    assert resolved.rejected_provided_context_count == 1


def test_site_sequence_records_missing_sequence_when_no_context_can_resolve() -> None:
    resolved = resolve_site_sequence_for_resolved_site(
        group=pd.DataFrame(
            {
                "site_string": ["S10"],
                "peptide_sequence": ["AAASAAA"],
                "multi_site": [False],
            }
        ),
        site_id="MAPK1;S10;",
        resolved_site_token="S10",
        multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
    )

    assert resolved.site_sequence is None
    assert resolved.source == SITE_SEQUENCE_SOURCE_MISSING


def test_summary_assembles_exact_fields_and_round_trips_payloads() -> None:
    input_metrics = PeptideEvidenceResolutionInputMetrics(
        peptide_observations_received=4,
        ambiguous_observations=2,
        unambiguous_observations=2,
        excluded_observations=0,
        split_observations=2,
        duplicate_peptide_rows=2,
        site_sequence_column_present=True,
        provided_site_sequence_count=3,
    )
    resolved_mapping = ResolvedMappingFractions(
        rows=pd.DataFrame(
            {
                "peptide_row_id": ["pep_split", "pep_split", "pep_single"],
                MAPPING_FRACTION_COLUMN: [0.5, 0.5, 1.0],
            }
        ),
        mapping_weight_source=DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT,
    )
    site_signals = SiteSignalSummary(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0, 3.0]},
            index=pd.Index(["AKT1;S473;", "MAPK1;S10;", "MAPK1;T12;"], name="site_id"),
        )
    )
    site_metadata_resolution = SiteMetadataResolution(
        site_metadata=pd.DataFrame(
            {
                "site_sequence": ["GGGSGGG", "AAASTTA", None],
            },
            index=pd.Index(["AKT1;S473;", "MAPK1;S10;", "MAPK1;T12;"], name="site_id"),
        ),
        sequence_diagnostics=SiteSequenceResolutionDiagnostics(
            rejected_provided_context_count=1,
            provided_site_sequence_used_count=1,
            peptide_context_derived_site_sequence_count=1,
            missing_site_sequence_count=1,
        ),
    )

    summary = build_resolution_summary(
        multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
        input_metrics=input_metrics,
        resolved_mapping=resolved_mapping,
        aggregation_policy=build_peptide_to_site_aggregation_policy(
            input_intensity_scale=IntensityScaleKind.LINEAR,
            mapping_rows=resolved_mapping.rows,
        ),
        site_signals=site_signals,
        site_metadata_resolution=site_metadata_resolution,
    )
    payload = summary.to_payload()

    assert payload == {
        "input_mode": "peptide_evidence",
        "multi_site_policy": DATASET_MULTI_SITE_POLICY_SPLIT,
        "peptide_to_site_aggregation_policy_id": (
            DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LINEAR_ALLOCATED_MEAN_V1
        ),
        "supported_input_scales": ["linear", "log2"],
        "supported_input_quantitative_meanings": [
            "peptide_abundance",
            "peptide_log2_abundance",
        ],
        "input_intensity_scale": "linear",
        "input_quantitative_meaning": "peptide_abundance",
        "output_intensity_scale": "linear",
        "output_quantitative_meaning": "phosphosite_abundance",
        "allocation_domain": "linear_abundance",
        "fractional_mapping_present": True,
        "peptide_observations_received": 4,
        "mapped_peptide_observations": 2,
        "site_mapping_rows": 3,
        "allocated_evidence_rows": 3,
        "unique_site_ids_produced": 3,
        "ambiguous_observations": 2,
        "unambiguous_observations": 2,
        "excluded_observations": 0,
        "split_observations": 2,
        "fractional_mapping_rows": 2,
        "unit_mapping_rows": 1,
        "mapping_weight_source_policy": (
            DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL
        ),
        "mapping_weight_normalization_policy": "sum_to_one_per_peptide_evidence_row",
        "mapping_weight_semantics": (
            "unitless_fraction_of_one_peptide_evidence_row_allocated_to_each_resolved_site"
        ),
        "signal_allocation_policy": (
            DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION
        ),
        "site_summarisation_policy": (
            DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS
        ),
        "missing_value_policy": DATASET_PEPTIDE_MISSING_VALUE_POLICY_FINITE_MEAN,
        "duplicate_evidence_policy": (
            "retain_duplicate_peptide_evidence_rows_as_separate_observations"
        ),
        "mixed_ambiguity_policy": (
            "combine_ambiguous_and_unambiguous_allocated_signals_in_site_mean"
        ),
        "localisation_aggregation_policy": (
            "arithmetic_mean_of_finite_reported_localisation_values"
        ),
        "localisation_summary_policy": (
            "descriptive_mean_of_finite_reported_localisation_confidence_values"
        ),
        "localisation_summary_semantics": (
            "descriptive_arithmetic_mean_not_calibrated_posterior_probability"
        ),
        "localisation_output_column": "localisation_confidence_descriptive_mean",
        "localisation_compatibility_alias_column": "localisation_confidence",
        "signal_conservation_policy": (
            DATASET_PEPTIDE_SIGNAL_CONSERVATION_POLICY_NOT_CONSERVED
        ),
        "uncertainty_limitations": [
            "no_model_based_uncertainty_or_posterior_localisation_combination",
            "peptide_evidence_rows_are_not_modelled_as_independent_replicates",
        ],
        "aggregation_policy": DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS,
        "aggregation_formula": (
            "a[p,s,j] [linear abundance units] = w[p,s] [unitless allocation "
            "fraction] * x[p,j] [linear peptide-abundance units]; y[s,j] "
            "[linear phosphosite-abundance estimate units] = arithmetic_mean("
            "a[p,s,j] over finite retained peptide evidence rows p mapped to "
            "site s for sample j)"
        ),
        "mapping_weight_source": DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT,
        "mapping_weight_normalisation": "sum_to_one_per_peptide_evidence_row",
        "duplicate_peptide_policy": (
            "retain_duplicate_peptide_evidence_rows_as_separate_observations"
        ),
        "duplicate_peptide_rows": 2,
        "site_sequence_column_present": True,
        "provided_site_sequence_count": 3,
        "accepted_site_sequence_count": 2,
        "rejected_site_sequence_count": 1,
        "provided_site_sequence_used_count": 1,
        "peptide_context_derived_site_sequence_count": 1,
        "missing_site_sequence_count": 1,
        "site_sequence_policy": (
            DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR
        ),
    }
    assert PeptideEvidenceResolutionSummary.from_payload(payload) == summary


def test_summary_normalizes_legacy_payload_compatibility_fields() -> None:
    legacy_payload: Mapping[str, object] = {
        "input_mode": "peptide_evidence",
        "multi_site_policy": DATASET_MULTI_SITE_POLICY_SPLIT,
        "peptide_observations_received": 1,
        "unique_site_ids_produced": 2,
        "ambiguous_observations": 1,
        "excluded_observations": 0,
        "split_observations": 1,
        "aggregation_policy": (
            DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN
        ),
        "aggregation_formula": "site_intensity = mean(per_peptide_intensity)",
        "mapping_weight_source": DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT,
        "mapping_weight_normalisation": "sum_to_one_per_peptide_row",
        "duplicate_peptide_policy": "retain_all_peptide_rows_as_independent_observations",
        "duplicate_peptide_rows": 0,
        "mixed_ambiguity_policy": (
            "mixed_ambiguous_and_unambiguous_rows_share_same_weighted_mean_aggregation"
        ),
        "site_sequence_column_present": False,
        "provided_site_sequence_count": 0,
        "accepted_site_sequence_count": 0,
        "rejected_site_sequence_count": 0,
        "provided_site_sequence_used_count": 0,
        "peptide_context_derived_site_sequence_count": 0,
        "missing_site_sequence_count": 2,
        "site_sequence_policy": DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR,
    }

    normalized = normalise_peptide_evidence_resolution_summary_payload(legacy_payload)
    summary = PeptideEvidenceResolutionSummary.from_payload(legacy_payload)

    assert normalized["mapping_weight_source_policy"] == (
        DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL
    )
    assert normalized["mapping_weight_normalization_policy"] == (
        "sum_to_one_per_peptide_evidence_row"
    )
    assert normalized["mapping_weight_normalisation"] == (
        "sum_to_one_per_peptide_evidence_row"
    )
    assert normalized["aggregation_policy"] == (
        DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS
    )
    assert summary.signal_allocation_policy == (
        DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION
    )
    assert summary.site_summarisation_policy == (
        DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS
    )


def test_resolver_outputs_are_isolated_from_input_frame_mutation() -> None:
    evidence_frame = pd.DataFrame(
        [
            {
                "peptide_row_id": "pep_1",
                "site_id": "MAPK1;S10;",
                "unique_feature_id": "feat_1",
                "gene_symbol": "MAPK1",
                "protein_accession": "P28482",
                "site_string": "S10",
                "sample_a": 10.0,
                "peptide_sequence": "AAASAAA",
                "modified_peptide_sequence": "AAASAAA",
                "multi_site": False,
                "provenance_source": "unit-test",
                "site_sequence": "AAASAAA",
            }
        ]
    )
    from phospy.science.evidence import PeptideEvidenceTable

    evidence = PeptideEvidenceTable(
        frame=evidence_frame,
        sample_intensity_columns=("sample_a",),
    )
    resolved = PeptideEvidenceDatasetResolver().run(
        evidence=evidence,
        multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
        input_intensity_scale="linear",
    )
    original_phospho = resolved.phospho.copy(deep=True)
    original_metadata = resolved.site_metadata.copy(deep=True)

    evidence_frame.loc[0, "sample_a"] = 999.0
    evidence_frame.loc[0, "site_sequence"] = "AAATAAA"

    pdt.assert_frame_equal(resolved.phospho, original_phospho)
    pdt.assert_frame_equal(resolved.site_metadata, original_metadata)

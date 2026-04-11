from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import phospy.motifs as motifs
from phospy.motifs import (
    KinaseMotifScorer,
    ValidatedMotifLibrary,
    create_frequency_matrix,
    frequency_scoring,
    score_phosphosite_motifs,
)
from phospy.references import (
    BundledReferenceProvider,
    ReferenceBundle,
    ReferenceBundleProvenance,
    ReferenceBundleSourceMetadata,
    ReferenceProvider,
)
from phospy.validation.errors import InputCompatibilityError, TableSchemaError


def test_create_frequency_matrix_normalizes_counts_and_ignores_gaps() -> None:
    frequency_mat = create_frequency_matrix(["A_A", "ACA"], flank_size=1)

    assert float(frequency_mat.loc["A", "p1"]) == pytest.approx(1.0)
    assert float(frequency_mat.loc["C", "p2"]) == pytest.approx(0.5)
    assert float(frequency_mat.loc["A", "p3"]) == pytest.approx(1.0)


def test_create_frequency_matrix_normalizes_lowercase_sequences() -> None:
    lower = create_frequency_matrix(["aca"], flank_size=1)
    upper = create_frequency_matrix(["ACA"], flank_size=1)

    pd.testing.assert_frame_equal(lower, upper)


def test_frequency_scoring_normalizes_lowercase_sequences() -> None:
    frequency_mat = create_frequency_matrix(["ACA", "AAA"], flank_size=1)

    result = frequency_scoring(
        sequence_list=pd.Series(["ACA", "aca"], index=["SITE_UPPER", "SITE_LOWER"]),
        frequency_mat=frequency_mat,
    )

    assert float(result.loc["SITE_UPPER"]) == pytest.approx(2.5)
    assert float(result.loc["SITE_LOWER"]) == pytest.approx(2.5)


def test_frequency_scoring_rejects_invalid_amino_acid_characters() -> None:
    frequency_mat = create_frequency_matrix(["ACA", "AAA"], flank_size=1)

    with pytest.raises(TableSchemaError, match="invalid amino-acid characters"):
        frequency_scoring(
            sequence_list=pd.Series(["AXA", np.nan], index=["SITE_X", "SITE_NA"]),
            frequency_mat=frequency_mat,
        )


def test_kinase_motif_scorer_extracts_centered_windows_and_scales_scores() -> None:
    scorer = KinaseMotifScorer.from_substrate_sequences(
        {
            "KINASE_A": ["AAAAA"],
            "KINASE_B": ["TTTTT"],
        },
        flank_size=2,
    )

    result = scorer.score_sequences(
        seqs={
            "SITE_A": "QQAAAAAYY",
            "SITE_B": "QQTTTTTYY",
        },
        min_motif_size=1,
    )

    assert list(result.motif_scores.columns) == ["KINASE_A", "KINASE_B"]
    assert float(result.motif_scores.loc["SITE_A", "KINASE_A"]) == pytest.approx(1.0)
    assert float(result.motif_scores.loc["SITE_A", "KINASE_B"]) == pytest.approx(0.0)
    assert result.sequence_windows.loc["SITE_A"] == "AAAAA"


def test_kinase_motif_scorer_encodes_sequence_windows_once_per_batch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    scorer = KinaseMotifScorer(
        motif_frequency_matrices={
            "KINASE_A": create_frequency_matrix(["AAAAA"], flank_size=2),
            "KINASE_B": create_frequency_matrix(["TTTTT"], flank_size=2),
            "KINASE_C": create_frequency_matrix(["SSSSS"], flank_size=2),
        },
        motif_sizes=pd.Series(
            {"KINASE_A": 2, "KINASE_B": 2, "KINASE_C": 2}, dtype=float
        ),
        flank_size=2,
    )

    encode_calls = 0
    original_encode = motifs._encode_sequence_positions

    def counting_encode(sequences: object, width: int) -> np.ndarray:
        nonlocal encode_calls
        encode_calls += 1
        return original_encode(sequences, width)

    monkeypatch.setattr(motifs, "_encode_sequence_positions", counting_encode)

    result = scorer.score_sequences(
        seqs={
            "SITE_A": "QQAAAAAYY",
            "SITE_B": "QQTTTTTYY",
            "SITE_C": "QQSSSSSYY",
        },
        min_motif_size=1,
    )

    assert encode_calls == 1
    assert list(result.motif_scores.columns) == ["KINASE_A", "KINASE_B", "KINASE_C"]


def test_score_phosphosite_motifs_filters_by_minimum_motif_size() -> None:
    motif_frequency_matrices = {
        "KINASE_A": create_frequency_matrix(["AAAAA"], flank_size=2),
        "KINASE_B": create_frequency_matrix(["TTTTT"], flank_size=2),
    }
    motif_sizes = pd.Series({"KINASE_A": 5, "KINASE_B": 1}, dtype=float)

    result = score_phosphosite_motifs(
        seqs={"SITE_A": "QQAAAAAYY"},
        motif_frequency_matrices=motif_frequency_matrices,
        motif_sizes=motif_sizes,
        min_motif_size=2,
        flank_size=2,
    )

    assert list(result.motif_scores.columns) == ["KINASE_A"]
    assert list(result.motif_sizes.index) == ["KINASE_A"]


def test_validated_motif_library_is_not_a_frozen_dataclass() -> None:
    assert ValidatedMotifLibrary.__dataclass_params__.frozen is False


def test_motif_scoring_result_is_detached_from_input_sequences() -> None:
    seqs = pd.Series({"SITE_A": "QQAAAAAYY", "SITE_B": "QQTTTTTYY"})
    original = seqs.copy(deep=True)
    scorer = KinaseMotifScorer.from_substrate_sequences(
        {"KINASE_A": ["AAAAA"], "KINASE_B": ["TTTTT"]},
        flank_size=2,
    )

    result = scorer.score_sequences(seqs=seqs, min_motif_size=1)

    result.motif_scores.loc["SITE_A", "KINASE_A"] = -999.0
    result.sequence_windows.loc["SITE_A"] = "CHANGED"

    pd.testing.assert_series_equal(seqs, original)


def test_reference_bundle_constructs_validated_kinase_prior_contract() -> None:
    bundle = ReferenceBundle(
        substrate_map={"KINASE_A": ["SITE_1", "SITE_2"]},
        motif_sequences={"KINASE_A": ["QQSQQ", "QQTQQ"]},
        species="human",
        source_metadata=ReferenceBundleSourceMetadata(
            source="bundled",
            reference="phosr-like",
            version="2026.04",
        ),
        provenance=ReferenceBundleProvenance(
            provider="BundledReferenceProvider",
            notes=("validated",),
        ),
    )

    assert bundle.species == "human"
    assert bundle.substrate_map == {"KINASE_A": ("SITE_1", "SITE_2")}
    assert bundle.motif_sequences == {"KINASE_A": ("QQSQQ", "QQTQQ")}


def test_reference_bundle_rejects_mismatched_kinase_sets() -> None:
    with pytest.raises(
        InputCompatibilityError,
        match=r"ReferenceBundle kinase sets must match exactly",
    ):
        ReferenceBundle(
            substrate_map={"KINASE_A": ["SITE_1"]},
            motif_sequences={"KINASE_B": ["QQSQQ"]},
            species="human",
            source_metadata=ReferenceBundleSourceMetadata(
                source="bundled",
                reference="phosr-like",
            ),
            provenance=ReferenceBundleProvenance(provider="BundledReferenceProvider"),
        )


def test_reference_bundle_rejects_empty_entries() -> None:
    with pytest.raises(
        InputCompatibilityError,
        match=r"ReferenceBundle substrate_map entries must not be empty: KINASE_A",
    ):
        ReferenceBundle(
            substrate_map={"KINASE_A": []},
            motif_sequences={"KINASE_A": ["QQSQQ"]},
            species="human",
            source_metadata=ReferenceBundleSourceMetadata(
                source="bundled",
                reference="phosr-like",
            ),
            provenance=ReferenceBundleProvenance(provider="BundledReferenceProvider"),
        )


def test_reference_provider_protocol_is_runtime_checkable() -> None:
    class StaticReferenceProvider:
        def resolve(
            self,
            *,
            species: str,
            reference: str = "auto",
        ) -> ReferenceBundle:
            return ReferenceBundle(
                substrate_map={"KINASE_A": ["SITE_1"]},
                motif_sequences={"KINASE_A": ["QQSQQ"]},
                species=species,
                source_metadata=ReferenceBundleSourceMetadata(
                    source="bundled",
                    reference=reference,
                ),
                provenance=ReferenceBundleProvenance(
                    provider="StaticReferenceProvider"
                ),
            )

    assert isinstance(StaticReferenceProvider(), ReferenceProvider)


def test_bundled_reference_provider_resolves_supported_rat_l6_bundle() -> None:
    provider = BundledReferenceProvider()

    bundle = provider.resolve(species="rat")

    assert isinstance(provider, ReferenceProvider)
    assert bundle.species == "rat"
    assert bundle.source_metadata.source == "phospy-bundled"
    assert bundle.source_metadata.reference == "l6_native"
    assert bundle.provenance.provider == "BundledReferenceProvider"
    assert "AKT1" in bundle.substrate_map
    assert bundle.substrate_map["AKT1"]
    assert len(bundle.motif_sequences["AKT1"]) == len(bundle.substrate_map["AKT1"])


def test_bundled_reference_provider_rejects_unsupported_species() -> None:
    provider = BundledReferenceProvider()

    with pytest.raises(
        InputCompatibilityError,
        match=r"Unsupported bundled reference species 'human'\. Supported species: rat",
    ):
        provider.resolve(species="human")


def test_bundled_reference_provider_rejects_unsupported_reference() -> None:
    provider = BundledReferenceProvider()

    with pytest.raises(
        InputCompatibilityError,
        match=(
            r"Unsupported bundled reference 'unknown' for species 'rat'\. "
            r"Supported references: l6_native"
        ),
    ):
        provider.resolve(species="rat", reference="unknown")

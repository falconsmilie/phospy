from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
)
from phospy.api import (
    DatasetBuildRequest,
    Organism,
    ReferenceBundle,
)
from phospy.errors import (
    DatasetValidationError,
    ReferenceValidationError,
    UnsupportedInputFormatError,
)
from phospy.references.models import ReferencePreset
from phospy.references.resolution import ReferenceResolver
from phospy.site_ids import canonicalize_site_identifier
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)


def test_builder_canonicalizes_site_ids_and_reorders_site_metadata() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index([" MAPK14;Y182; ", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "site_id": ["AKT1;T308;", " MAPK14;Y182; "],
            "gene_symbol": ["AKT1", "MAPK14"],
            "site": ["T308", "Y182"],
            "site_sequence": ["A" * 31, "B" * 31],
        }
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )

    assert list(built.phospho.index) == ["MAPK14;Y182;", "AKT1;T308;"]
    assert list(built.site_metadata.index) == ["MAPK14;Y182;", "AKT1;T308;"]
    assert built.site_metadata.loc["MAPK14;Y182;", "gene_symbol"] == "MAPK14"
    assert built.site_metadata.loc["AKT1;T308;", "gene_symbol"] == "AKT1"


def test_builder_canonicalizes_lowercase_site_ids() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0],
            "sample_b": [2.0],
        },
        index=pd.Index(["mapk14;y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "site_id": ["mapk14;y182;"],
            "gene_symbol": ["mapk14"],
            "site": ["y182"],
            "site_sequence": ["A" * 31],
        }
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )

    assert list(built.phospho.index) == ["MAPK14;Y182;"]
    assert list(built.site_metadata.index) == ["MAPK14;Y182;"]


def test_builder_rejects_ambiguous_site_ids_after_canonicalization() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0]},
        index=pd.Index(["mapk14;y182;", "MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["A" * 31, "B" * 31],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        UnsupportedInputFormatError,
        match="duplicate site identifiers after canonicalization",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
            )
        )


def test_reference_bundle_rejects_ambiguous_site_sequence_ids() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="conflicting site_sequence values after canonicalization",
    ) as exc_info:
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {"kinase": ["MAP2K6"], "substrate_site": ["MAPK14;Y182;"]}
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31, "B" * 31]},
                index=pd.Index(
                    ["MAPK14;Y182;", " MAPK14;Y182;"],
                    name="site_id",
                ),
            ),
        )
    report = getattr(exc_info.value, "identifier_normalisation_report", None)
    assert report is not None
    assert report.conflict_count > 0


def test_dataset_boundary_rejects_non_canonical_site_ids() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="site identifiers must use 'GENE;SITE;' format",
    ):
        from phospy.datasets.models import AnalysisReadyPhosphoDataset

        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame(
                {
                    "sample_a": [1.0, 3.0],
                    "sample_b": [2.0, 4.0],
                },
                index=pd.Index([101, 202], name="site_id"),
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "AKT1"],
                    "site": ["Y182", "T308"],
                    "site_sequence": ["A" * 31, "B" * 31],
                },
                index=pd.Index([101, 202], name="site_id"),
            ),
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        )


def test_dataset_boundary_rejects_lowercase_site_ids() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="dataset\\.phospho\\.index must contain canonical site identifiers",
    ):
        from phospy.datasets.models import AnalysisReadyPhosphoDataset

        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame(
                {
                    "sample_a": [1.0],
                    "sample_b": [2.0],
                },
                index=pd.Index(["mapk14;y182;"], name="site_id"),
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["A" * 31],
                },
                index=pd.Index(["mapk14;y182;"], name="site_id"),
            ),
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        )


def test_dataset_boundary_rejects_whitespace_site_ids() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="dataset\\.phospho\\.index must contain canonical site identifiers",
    ):
        from phospy.datasets.models import AnalysisReadyPhosphoDataset

        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame(
                {
                    "sample_a": [1.0],
                    "sample_b": [2.0],
                },
                index=pd.Index([" MAPK14;Y182; "], name="site_id"),
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["A" * 31],
                },
                index=pd.Index([" MAPK14;Y182; "], name="site_id"),
            ),
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        )


def test_dataset_boundary_rejects_missing_trailing_delimiter_site_ids() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="dataset\\.phospho\\.index must contain canonical site identifiers",
    ):
        from phospy.datasets.models import AnalysisReadyPhosphoDataset

        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame(
                {
                    "sample_a": [1.0],
                    "sample_b": [2.0],
                },
                index=pd.Index(["MAPK14;Y182"], name="site_id"),
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["A" * 31],
                },
                index=pd.Index(["MAPK14;Y182"], name="site_id"),
            ),
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        )


def test_dataset_boundary_rejects_duplicates_after_site_id_canonicalization() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="contains duplicate site identifiers after canonicalization",
    ):
        from phospy.datasets.models import AnalysisReadyPhosphoDataset

        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame(
                {
                    "sample_a": [1.0, 3.0],
                    "sample_b": [2.0, 4.0],
                },
                index=pd.Index(
                    ["mapk14;y182;", "MAPK14;Y182;"],
                    name="site_id",
                ),
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "MAPK14"],
                    "site": ["Y182", "Y182"],
                    "site_sequence": ["A" * 31, "B" * 31],
                },
                index=pd.Index(
                    ["mapk14;y182;", "MAPK14;Y182;"],
                    name="site_id",
                ),
            ),
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        )


def test_dataset_boundary_rejects_colliding_dirty_site_ids() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="contains duplicate site identifiers after canonicalization",
    ):
        from phospy.datasets.models import AnalysisReadyPhosphoDataset

        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame(
                {
                    "sample_a": [1.0, 3.0],
                    "sample_b": [2.0, 4.0],
                },
                index=pd.Index(
                    ["MAPK14;Y182;", " MAPK14;Y182;"],
                    name="site_id",
                ),
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "MAPK14"],
                    "site": ["Y182", "Y182"],
                    "site_sequence": ["A" * 31, "B" * 31],
                },
                index=pd.Index(
                    ["MAPK14;Y182;", " MAPK14;Y182;"],
                    name="site_id",
                ),
            ),
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        )


def test_dataset_boundary_accepts_strict_canonical_site_ids() -> None:
    from phospy.datasets.models import AnalysisReadyPhosphoDataset

    dataset = AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {
                "sample_a": [1.0],
                "sample_b": [2.0],
            },
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["A" * 31],
            },
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )

    assert list(dataset.phospho.index) == ["MAPK14;Y182;"]


def test_dataset_boundary_requires_explicit_intensity_and_processing_state() -> None:
    with pytest.raises(TypeError, match="missing .* required positional argument"):
        from phospy.datasets.models import AnalysisReadyPhosphoDataset

        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame(
                {"sample_a": [1.0]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["A" * 31],
                },
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            organism=Organism.RAT,
        )  # type: ignore[call-arg]


def test_reference_bundle_rejects_duplicate_kinase_substrate_pairs() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="contains duplicate \\(kinase, substrate_site\\) pairs",
    ):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {
                    "kinase": ["MAP2K6", "MAP2K6"],
                    "substrate_site": ["MAPK14;Y182;", "MAPK14;Y182;"],
                }
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
        )


def test_reference_bundle_rejects_colliding_dirty_substrate_site_ids() -> None:
    bundle = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K7"],
                "substrate_site": ["MAPK14;Y182;", " mapk14 ; y182 "],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
    )
    assert set(bundle.kinase_substrate_map.loc[:, "substrate_site"]) == {"MAPK14;Y182;"}


def test_reference_bundle_rejects_duplicate_pairs_after_site_id_normalization() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="contains duplicate \\(kinase, substrate_site\\) pairs",
    ):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {
                    "kinase": ["MAP2K6", "MAP2K6"],
                    "substrate_site": ["MAPK14;Y182;", " mapk14 ; y182 "],
                }
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
        )


def test_reference_bundle_normalises_mixed_case_kinase_ids() -> None:
    bundle = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["akt1", "Akt1", "AKT1"],
                "substrate_site": [
                    "MAPK14;Y182;",
                    "MAPK14;T185;",
                    "GSK3B;S9;",
                ],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31, "B" * 31, "C" * 31]},
            index=pd.Index(
                ["MAPK14;Y182;", "MAPK14;T185;", "GSK3B;S9;"],
                name="site_id",
            ),
        ),
    )
    assert set(bundle.kinase_substrate_map.loc[:, "kinase"]) == {"AKT1"}


def test_reference_bundle_rejects_duplicate_pairs_after_kinase_normalization() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="contains duplicate \\(kinase, substrate_site\\) pairs",
    ):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {
                    "kinase": ["akt1", "AKT1"],
                    "substrate_site": ["MAPK14;Y182;", "MAPK14;Y182;"],
                }
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
        )


def test_reference_provider_shapes_bundled_resources_for_strict_boundary() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )

    assert bundle.site_sequences.index.is_unique
    assert (
        bundle.kinase_substrate_map.duplicated(
            subset=["kinase", "substrate_site"]
        ).sum()
        == 0
    )
    assert all(
        isinstance(value, str) and value == value.strip() and value != ""
        for value in bundle.site_sequences.index.tolist()
    )
    assert all(
        isinstance(value, str) and value == value.strip() and value != ""
        for value in bundle.kinase_substrate_map.loc[:, "substrate_site"].tolist()
    )


def test_shared_normalizer_supports_whitespace_case_and_missing_trailing_semicolon() -> (
    None
):
    assert (
        canonicalize_site_identifier(
            " mapk1 ; s123 ",
            field_name="test.site_id",
            error_type=UnsupportedInputFormatError,
        )
        == "MAPK1;S123;"
    )
    assert (
        canonicalize_site_identifier(
            "mapk1;s123",
            field_name="test.site_id",
            error_type=UnsupportedInputFormatError,
        )
        == "MAPK1;S123;"
    )


def test_shared_normalizer_rejects_malformed_site_id_with_clear_error() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="site identifiers must use 'GENE;SITE;' format",
    ):
        canonicalize_site_identifier(
            "MAPK14-Y182",
            field_name="test.site_id",
            error_type=UnsupportedInputFormatError,
        )


def test_dataset_and_reference_ids_align_after_shared_normalization() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=pd.DataFrame(
                {"sample_a": [1.0], "sample_b": [2.0]},
                index=pd.Index([" mapk14 ; y182 "], name="site_id"),
            ),
            site_metadata=pd.DataFrame(
                {
                    "site_id": ["MAPK14;Y182"],
                    "gene_symbol": ["mapk14"],
                    "site": ["y182"],
                    "site_sequence": ["A" * 31],
                }
            ),
            organism=Organism.RAT,
        )
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["MAP2K6"], "substrate_site": [" mapk14 ; y182 "]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31]},
            index=pd.Index(["MAPK14;Y182"], name="site_id"),
        ),
    )
    overlap = KinaseWorkflowInterpreter._summarize_overlap(
        dataset=built.phospho,
        kinase_substrate_map=references.kinase_substrate_map,
    )
    scoring_index = KinaseWorkflowInterpreter._resolve_scoring_site_index(
        dataset=built.phospho,
        site_sequences=references.site_sequences,
    )
    assert overlap["overlap_sites"] == 1
    assert list(scoring_index) == ["MAPK14;Y182;"]

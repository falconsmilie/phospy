from __future__ import annotations

from collections.abc import Mapping
from itertools import permutations

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import DatasetBuildRequest, Organism, PhosPyInputError
from phospy.api.requests import (
    DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
    DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
    DATASET_MULTI_SITE_POLICY_REJECT,
    DATASET_MULTI_SITE_POLICY_SPLIT,
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
    DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
)
from phospy.provenance.serialization import from_payload
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.evidence import (
    DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS,
    DATASET_PEPTIDE_DUPLICATE_POLICY_RETAIN_ALL_ROWS,
    DATASET_PEPTIDE_LOCALISATION_AGGREGATION_POLICY_ARITHMETIC_MEAN_OF_FINITE_VALUES,
    DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE,
    DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW,
    DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT,
    DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL,
    DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS,
    DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION,
    DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR,
    DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS,
    DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS,
    DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN,
    PeptideEvidenceDatasetResolver,
    PeptideEvidenceResolutionResult,
    PeptideEvidenceResolutionSummary,
    PeptideEvidenceTable,
    build_multi_site_handling_config_for_dataset_policy,
)
from phospy.science.evidence.dataset_resolution import (
    _single_non_empty_string_or_error,
)
from phospy.science.sites.site_keys import decode_site_key
from phospy.validation.datasets.builder_request import DatasetBuildRequestValidator


def _site_key_for_display_id(
    dataset: AnalysisReadyPhosphoDataset,
    display_id: str,
) -> str:
    matches = dataset.site_metadata.index[
        dataset.site_metadata.loc[:, "display_id"].astype(str) == display_id
    ].astype(str)
    assert len(matches) == 1
    return str(matches[0])


def _display_ids(dataset: AnalysisReadyPhosphoDataset) -> list[str]:
    return dataset.site_metadata.loc[:, "display_id"].astype(str).tolist()


def _site_level_phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )


def _site_level_site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            "protein_id": ["MAPK14"],
        },
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )


def _peptide_evidence_frame(*, include_single_site: bool = True) -> pd.DataFrame:
    rows = [
        {
            "peptide_row_id": "pep_joint",
            "site_id": "MAPK1;S10;",
            "unique_feature_id": "feat_joint",
            "gene_symbol": "MAPK1",
            "protein_accession": "P28482",
            "site_string": "S10,T12",
            "sample_a": 10.0,
            "sample_b": 12.0,
            "peptide_sequence": "AAAAA",
            "modified_peptide_sequence": "AA[+80]AAA",
            "multi_site": True,
            "provenance_source": "maxquant",
            "site_sequence": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAA",
            "localisation_confidence": 0.95,
        }
    ]
    if include_single_site:
        rows.append(
            {
                "peptide_row_id": "pep_single",
                "site_id": "AKT1;S473;",
                "unique_feature_id": "feat_single",
                "gene_symbol": "AKT1",
                "protein_accession": "P31749",
                "site_string": "S473",
                "sample_a": 7.0,
                "sample_b": 9.0,
                "peptide_sequence": "BBBBB",
                "modified_peptide_sequence": "BB[+80]BBB",
                "multi_site": False,
                "provenance_source": "maxquant",
                "site_sequence": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                "localisation_confidence": 0.9,
            }
        )
    return pd.DataFrame(rows)


def _resolution_policy_evidence_row(
    *,
    peptide_row_id: str,
    sample_a: float,
    site_id: str = "MAPK1;S10;",
    site_string: str = "S10",
    peptide_sequence: str = "AAASAAA",
    multi_site: bool = False,
    localisation_confidence: object = 0.9,
) -> dict[str, object]:
    return {
        "peptide_row_id": peptide_row_id,
        "site_id": site_id,
        "unique_feature_id": f"feat_{peptide_row_id}",
        "gene_symbol": "MAPK1",
        "protein_accession": "P28482",
        "site_string": site_string,
        "sample_a": sample_a,
        "peptide_sequence": peptide_sequence,
        "modified_peptide_sequence": peptide_sequence,
        "multi_site": multi_site,
        "provenance_source": "unit-test",
        "localisation_confidence": localisation_confidence,
    }


def _resolution_policy_evidence_table(
    rows: list[dict[str, object]],
    *,
    split: bool = False,
    site_mapping: pd.DataFrame | None = None,
) -> PeptideEvidenceTable:
    return PeptideEvidenceTable(
        frame=pd.DataFrame(rows),
        sample_intensity_columns=("sample_a",),
        site_mapping=site_mapping,
        multi_site_handling_config=(
            build_multi_site_handling_config_for_dataset_policy(
                multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT
            )
            if split
            else None
        ),
    )


def _run_policy_resolution(
    rows: list[dict[str, object]],
    *,
    split: bool = False,
    site_mapping: pd.DataFrame | None = None,
) -> PeptideEvidenceResolutionResult:
    return PeptideEvidenceDatasetResolver().run(
        evidence=_resolution_policy_evidence_table(
            rows,
            split=split,
            site_mapping=site_mapping,
        ),
        multi_site_policy=(
            DATASET_MULTI_SITE_POLICY_SPLIT
            if split
            else DATASET_MULTI_SITE_POLICY_KEEP_JOINT
        ),
        input_intensity_scale="linear",
    )


def _single_site_peptide_evidence_frame() -> pd.DataFrame:
    evidence = _peptide_evidence_frame(include_single_site=True)
    return evidence.loc[
        evidence.loc[:, "peptide_row_id"] == "pep_single",
        :,
    ].reset_index(drop=True)


def _peptide_evidence_center_residue_mismatch_frame() -> pd.DataFrame:
    evidence = _peptide_evidence_frame(include_single_site=True)
    single_site = evidence.loc[
        evidence.loc[:, "peptide_row_id"] == "pep_single", :
    ].copy(deep=True)
    mismatched_sequence = "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA"
    assert len(mismatched_sequence) % 2 == 1
    assert mismatched_sequence[len(mismatched_sequence) // 2] == "T"
    single_site.loc[:, "site_id"] = "AKT1;S473;"
    single_site.loc[:, "site_string"] = "S473"
    single_site.loc[:, "site_sequence"] = mismatched_sequence
    return single_site.reset_index(drop=True)


def _same_resolved_site_evidence_frame(
    site_sequences: tuple[object, ...],
) -> pd.DataFrame:
    base_row = _single_site_peptide_evidence_frame().iloc[0].to_dict()
    peptide_sequences = (
        "AAASAAA",
        "CCCSCCC",
        "GGGSGGG",
        "HHHSHHH",
        "KKKSKKK",
    )
    rows: list[dict[str, object]] = []
    for index, site_sequence in enumerate(site_sequences, start=1):
        peptide_sequence = peptide_sequences[(index - 1) % len(peptide_sequences)]
        row = dict(base_row)
        row["peptide_row_id"] = f"pep_single_{index}"
        row["unique_feature_id"] = f"feat_single_{index}"
        row["sample_a"] = float(6 + index)
        row["sample_b"] = float(8 + index)
        row["peptide_sequence"] = peptide_sequence
        row["modified_peptide_sequence"] = peptide_sequence
        row["site_sequence"] = site_sequence
        rows.append(row)
    return pd.DataFrame(rows)


def _split_multisite_ambiguous_derivation_frame() -> pd.DataFrame:
    base_row = _peptide_evidence_frame(include_single_site=False).iloc[0].to_dict()
    rows: list[dict[str, object]] = []
    for index, peptide_sequence in enumerate(("AAASTTAAAA", "CCCSATCCCC"), start=1):
        row = dict(base_row)
        row["peptide_row_id"] = f"pep_joint_{index}"
        row["unique_feature_id"] = f"feat_joint_{index}"
        row["site_string"] = "S10,T12"
        row["sample_a"] = float(10 * index)
        row["sample_b"] = float(12 * index)
        row["peptide_sequence"] = peptide_sequence
        row["modified_peptide_sequence"] = peptide_sequence
        row["site_sequence"] = "AAAAASAAAAA"
        rows.append(row)
    return pd.DataFrame(rows)


def _assert_site_sequence_center_mismatch_message(message: str) -> None:
    lower_message = message.lower()
    assert "site_sequence" in message
    assert "centre" in lower_message or "center" in lower_message
    assert "expected='S'" in message
    assert "observed='T'" in message
    assert "AKT1;S473;" in message


def _assert_site_sequence_conflict_message(message: str) -> None:
    assert "site_sequence values conflict" in message
    assert "site_id='AKT1;S473;'" in message
    assert "distinct_normalized_value_count=2" in message
    assert "values=['AAASAAA', 'CCCSCCC']" in message
    assert "row order" in message
    assert "Correct the source evidence" in message
    assert "explicit upstream reference-resolution policy" in message


def _assert_log2_fractional_allocation_message(message: str) -> None:
    assert "peptide-evidence mode" in message
    assert "input_intensity_scale='log2'" in message
    assert "fractional allocation" in message
    assert "non-unit" in message
    assert "Supported corrective action" in message


def test_site_level_input_works_with_safe_default_or_explicit_declaration() -> None:
    safe_default_request = DatasetBuildRequest(
        phospho=_site_level_phospho(),
        site_metadata=_site_level_site_metadata(),
        input_intensity_scale="linear",
    )
    explicit_request = DatasetBuildRequest(
        phospho=_site_level_phospho(),
        site_metadata=_site_level_site_metadata(),
        input_intensity_scale="linear",
        site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
    )
    assert (
        DatasetBuildRequestValidator().run(safe_default_request) is safe_default_request
    )
    assert DatasetBuildRequestValidator().run(explicit_request) is explicit_request


def test_peptide_evidence_requires_multi_site_policy() -> None:
    request = DatasetBuildRequest(
        site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
        peptide_evidence=_peptide_evidence_frame(),
        peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
        input_intensity_scale="linear",
        organism=Organism.HUMAN,
    )
    with pytest.raises(
        PhosPyInputError,
        match="peptide_evidence input requires multi_site_policy",
    ):
        DatasetBuildRequestValidator().run(request)


def test_peptide_evidence_resolver_preserves_accession_identity_metadata() -> None:
    resolved = PeptideEvidenceDatasetResolver().run(
        evidence=PeptideEvidenceTable(
            frame=_single_site_peptide_evidence_frame(),
            sample_intensity_columns=("sample_a", "sample_b"),
        ),
        multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
        input_intensity_scale="linear",
    )

    assert "protein_accession" in resolved.site_metadata.columns
    site_row = resolved.site_metadata.loc["AKT1;S473;", :]
    assert site_row["protein_accession"] == "P31749"
    assert "protein_id" not in resolved.site_metadata.columns


def test_peptide_evidence_rejects_conflicting_accessions_per_resolved_site() -> None:
    evidence = pd.concat(
        [
            _single_site_peptide_evidence_frame(),
            _single_site_peptide_evidence_frame(),
        ],
        ignore_index=True,
    )
    evidence.loc[1, "peptide_row_id"] = "pep_single_conflict"
    evidence.loc[1, "unique_feature_id"] = "feat_single_conflict"
    evidence.loc[1, "protein_accession"] = "Q9Y243"
    evidence.loc[1, "peptide_sequence"] = "CCCCC"
    evidence.loc[1, "modified_peptide_sequence"] = "CC[+80]CCC"

    with pytest.raises(PhosPyInputError) as exc_info:
        PeptideEvidenceDatasetResolver().run(
            evidence=PeptideEvidenceTable(
                frame=evidence,
                sample_intensity_columns=("sample_a", "sample_b"),
            ),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
            input_intensity_scale="linear",
        )
    message = str(exc_info.value)
    assert "site_id='AKT1;S473;'" in message
    assert "protein_accession" in message
    assert "P31749" in message
    assert "Q9Y243" in message
    assert "disambiguate peptide-site mapping or split rows before building" in message


def test_accession_conflict_helper_returns_none_for_all_empty_values() -> None:
    assert (
        _single_non_empty_string_or_error(
            pd.Series(["", "  ", None, pd.NA]),
            field_name="protein_accession",
            site_id="AKT1;S473;",
        )
        is None
    )


def test_builder_uses_peptide_evidence_accession_for_site_key_identity() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=_single_site_peptide_evidence_frame(),
            peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
            input_intensity_scale="linear",
            organism=Organism.HUMAN,
        )
    )

    assert built.phospho.index.name == "site_key"
    assert built.site_metadata.index.name == "site_key"
    assert "protein_namespace" in built.site_metadata.columns
    assert "protein_identifier" in built.site_metadata.columns
    site_key = _site_key_for_display_id(built, "AKT1;S473;")
    site_row = built.site_metadata.loc[site_key, :]
    assert site_row["protein_namespace"] == "protein_accession"
    assert site_row["protein_identifier"] == "P31749"
    assert "protein_accession" in built.site_metadata.columns
    assert site_row["protein_accession"] == "P31749"

    decoded = decode_site_key(
        site_key,
        field_name="test.builder.peptide_evidence.site_key",
        error_type=ValueError,
    )
    assert decoded.organism == "human"
    assert decoded.protein_namespace == "protein_accession"
    assert decoded.protein_identifier == "P31749"
    assert decoded.residue == "S"
    assert decoded.position == 473


def test_public_builder_rejects_log2_peptide_evidence_default_equal_split() -> None:
    dataset: AnalysisReadyPhosphoDataset | None = None

    with pytest.raises(PhosPyInputError) as exc_info:
        dataset = AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
                peptide_evidence=_peptide_evidence_frame(include_single_site=False),
                peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
                multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
                input_intensity_scale="log2",
                organism=Organism.HUMAN,
            )
        )

    assert dataset is None
    _assert_log2_fractional_allocation_message(str(exc_info.value))


def test_public_builder_keeps_linear_peptide_evidence_default_equal_split_behavior() -> (
    None
):
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=_peptide_evidence_frame(include_single_site=False),
            peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
            input_intensity_scale="linear",
            organism=Organism.HUMAN,
        )
    )

    mapk1_s10 = _site_key_for_display_id(built, "MAPK1;S10;")
    mapk1_t12 = _site_key_for_display_id(built, "MAPK1;T12;")
    assert float(built.phospho.loc[mapk1_s10, "sample_a"]) == pytest.approx(5.0)
    assert float(built.phospho.loc[mapk1_t12, "sample_a"]) == pytest.approx(5.0)
    assert float(built.phospho.loc[mapk1_s10, "sample_b"]) == pytest.approx(6.0)
    assert float(built.phospho.loc[mapk1_t12, "sample_b"]) == pytest.approx(6.0)


def test_public_builder_allows_log2_peptide_evidence_with_unit_mapping() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=_single_site_peptide_evidence_frame(),
            peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
            input_intensity_scale="log2",
            organism=Organism.HUMAN,
        )
    )

    akt1_s473 = _site_key_for_display_id(built, "AKT1;S473;")
    assert float(built.phospho.loc[akt1_s473, "sample_a"]) == pytest.approx(7.0)
    assert built.intensity_scale_state.label == "log2"


@pytest.mark.parametrize("accession_case", ("missing_column", "blank_value"))
def test_peptide_evidence_requires_accession_before_resolution(
    accession_case: str,
) -> None:
    evidence = _single_site_peptide_evidence_frame()
    if accession_case == "missing_column":
        evidence = evidence.drop(columns=["protein_accession"])
    else:
        evidence.loc[:, "protein_accession"] = ""

    with pytest.raises(PhosPyInputError, match="protein_accession"):
        PeptideEvidenceTable(
            frame=evidence,
            sample_intensity_columns=("sample_a", "sample_b"),
        )


def test_peptide_evidence_rejects_site_sequence_center_residue_mismatch() -> None:
    with pytest.raises(PhosPyInputError) as exc_info:
        PeptideEvidenceDatasetResolver().run(
            evidence=PeptideEvidenceTable(
                frame=_peptide_evidence_center_residue_mismatch_frame(),
                sample_intensity_columns=("sample_a", "sample_b"),
            ),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
            input_intensity_scale="linear",
        )

    _assert_site_sequence_center_mismatch_message(str(exc_info.value))


def test_builder_rejects_peptide_evidence_sequence_center_residue_mismatch() -> None:
    dataset: AnalysisReadyPhosphoDataset | None = None

    with pytest.raises(PhosPyInputError) as exc_info:
        dataset = AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
                peptide_evidence=_peptide_evidence_center_residue_mismatch_frame(),
                peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
                multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
                input_intensity_scale="linear",
                organism=Organism.HUMAN,
            )
        )

    assert dataset is None
    _assert_site_sequence_center_mismatch_message(str(exc_info.value))


def test_peptide_evidence_rejects_conflicting_valid_site_sequences() -> None:
    evidence = _same_resolved_site_evidence_frame(("AAASAAA", "CCCSCCC"))

    with pytest.raises(PhosPyInputError) as exc_info:
        PeptideEvidenceDatasetResolver().run(
            evidence=PeptideEvidenceTable(
                frame=evidence,
                sample_intensity_columns=("sample_a", "sample_b"),
            ),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
            input_intensity_scale="linear",
        )

    _assert_site_sequence_conflict_message(str(exc_info.value))


def test_peptide_evidence_conflicting_site_sequence_details_are_order_invariant() -> (
    None
):
    messages: list[str] = []
    for site_sequences in (("AAASAAA", "CCCSCCC"), ("CCCSCCC", "AAASAAA")):
        with pytest.raises(PhosPyInputError) as exc_info:
            PeptideEvidenceDatasetResolver().run(
                evidence=PeptideEvidenceTable(
                    frame=_same_resolved_site_evidence_frame(site_sequences),
                    sample_intensity_columns=("sample_a", "sample_b"),
                ),
                multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
                input_intensity_scale="linear",
            )
        messages.append(str(exc_info.value))

    assert messages[0] == messages[1]
    _assert_site_sequence_conflict_message(messages[0])


def test_peptide_evidence_equivalent_site_sequences_normalize_to_one_value() -> None:
    resolved = PeptideEvidenceDatasetResolver().run(
        evidence=PeptideEvidenceTable(
            frame=_same_resolved_site_evidence_frame((" aaAsaaa ", "AAASAAA")),
            sample_intensity_columns=("sample_a", "sample_b"),
        ),
        multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
        input_intensity_scale="linear",
    )

    assert resolved.site_metadata.loc["AKT1;S473;", "site_sequence"] == "AAASAAA"
    payload = resolved.summary.to_payload()
    assert int(payload["provided_site_sequence_count"]) == 2
    assert int(payload["accepted_site_sequence_count"]) == 1
    assert int(payload["rejected_site_sequence_count"]) == 0
    assert int(payload["provided_site_sequence_used_count"]) == 1


def test_peptide_evidence_rejects_mixed_valid_and_invalid_site_sequences() -> None:
    evidence = _same_resolved_site_evidence_frame(("AAASAAA", "AAATAAA"))

    with pytest.raises(PhosPyInputError) as exc_info:
        PeptideEvidenceDatasetResolver().run(
            evidence=PeptideEvidenceTable(
                frame=evidence,
                sample_intensity_columns=("sample_a", "sample_b"),
            ),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
            input_intensity_scale="linear",
        )

    message = str(exc_info.value)
    assert "Mixed valid and invalid supplied evidence" in message
    assert "site_id='AKT1;S473;'" in message
    assert "valid_normalized_value='AAASAAA'" in message
    assert "AAATAAA" in message
    assert "expected='S'" in message
    assert "observed='T'" in message
    assert "explicit upstream reference-resolution policy" in message


def test_peptide_evidence_ignores_blank_and_null_sequence_values_for_conflicts() -> (
    None
):
    resolved = PeptideEvidenceDatasetResolver().run(
        evidence=PeptideEvidenceTable(
            frame=_same_resolved_site_evidence_frame(("AAASAAA", " ", pd.NA, None)),
            sample_intensity_columns=("sample_a", "sample_b"),
        ),
        multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
        input_intensity_scale="linear",
    )

    assert resolved.site_metadata.loc["AKT1;S473;", "site_sequence"] == "AAASAAA"
    payload = resolved.summary.to_payload()
    assert int(payload["provided_site_sequence_count"]) == 1
    assert int(payload["accepted_site_sequence_count"]) == 1
    assert int(payload["rejected_site_sequence_count"]) == 0
    assert int(payload["provided_site_sequence_used_count"]) == 1


def test_split_policy_derives_site_specific_context_when_shared_window_mismatches() -> (
    None
):
    evidence = _peptide_evidence_frame(include_single_site=False)
    evidence.loc[:, "site_string"] = "S10,T12"
    evidence.loc[:, "peptide_sequence"] = "AAASTTAAAA"
    evidence.loc[:, "modified_peptide_sequence"] = "AAA(ph)ST(ph)TAAAA"
    evidence.loc[:, "site_sequence"] = "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA"

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=evidence,
            peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
            input_intensity_scale="linear",
            organism=Organism.HUMAN,
        )
    )

    mapk1_s10 = _site_key_for_display_id(built, "MAPK1;S10;")
    mapk1_t12 = _site_key_for_display_id(built, "MAPK1;T12;")
    assert built.site_metadata.loc[mapk1_s10, "site_sequence"] == (
        "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA"
    )
    assert built.site_metadata.loc[mapk1_t12, "site_sequence"] == "AASTTAAAA"
    assert built.provenance is not None
    payload = built.provenance.workflow_parameters["peptide_evidence_resolution"]
    assert isinstance(payload, Mapping)
    assert int(payload["provided_site_sequence_count"]) == 1
    assert int(payload["accepted_site_sequence_count"]) == 2
    assert int(payload["rejected_site_sequence_count"]) == 1
    assert int(payload["provided_site_sequence_used_count"]) == 1
    assert int(payload["peptide_context_derived_site_sequence_count"]) == 1
    assert int(payload["missing_site_sequence_count"]) == 0


def test_split_policy_rejects_non_unique_peptide_context_derivation() -> None:
    with pytest.raises(PhosPyInputError) as exc_info:
        PeptideEvidenceDatasetResolver().run(
            evidence=PeptideEvidenceTable(
                frame=_split_multisite_ambiguous_derivation_frame(),
                sample_intensity_columns=("sample_a", "sample_b"),
                multi_site_handling_config=(
                    build_multi_site_handling_config_for_dataset_policy(
                        multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT
                    )
                ),
            ),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
            input_intensity_scale="linear",
        )

    message = str(exc_info.value)
    assert "site_id='MAPK1;T12;'" in message
    assert "invalid_supplied_values" in message
    assert "derived_candidate_count=2" in message
    assert "derived_candidates=['AASTTAAAA', 'CCSATCCCC']" in message
    assert "explicit upstream reference-resolution policy" in message


def test_peptide_evidence_resolution_is_invariant_to_row_permutation() -> None:
    expected_phospho: pd.DataFrame | None = None
    expected_site_metadata: pd.DataFrame | None = None
    expected_payload: dict[str, object] | None = None
    for site_sequences in permutations((" aaAsaaa ", pd.NA, "AAASAAA")):
        resolved = PeptideEvidenceDatasetResolver().run(
            evidence=PeptideEvidenceTable(
                frame=_same_resolved_site_evidence_frame(tuple(site_sequences)),
                sample_intensity_columns=("sample_a", "sample_b"),
            ),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
            input_intensity_scale="linear",
        )
        payload = resolved.summary.to_payload()
        if expected_phospho is None:
            expected_phospho = resolved.phospho
            expected_site_metadata = resolved.site_metadata
            expected_payload = payload
            continue
        pd.testing.assert_frame_equal(resolved.phospho, expected_phospho)
        assert expected_site_metadata is not None
        pd.testing.assert_frame_equal(resolved.site_metadata, expected_site_metadata)
        assert payload == expected_payload

    failure_messages: set[str] = set()
    for site_sequences in permutations(("AAASAAA", "CCCSCCC", " ")):
        with pytest.raises(PhosPyInputError) as exc_info:
            PeptideEvidenceDatasetResolver().run(
                evidence=PeptideEvidenceTable(
                    frame=_same_resolved_site_evidence_frame(tuple(site_sequences)),
                    sample_intensity_columns=("sample_a", "sample_b"),
                ),
                multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
                input_intensity_scale="linear",
            )
        failure_messages.add(str(exc_info.value))

    assert len(failure_messages) == 1
    _assert_site_sequence_conflict_message(next(iter(failure_messages)))


def test_peptide_evidence_preserves_matching_sequence_with_only_text_normalisation() -> (
    None
):
    evidence = _peptide_evidence_frame(include_single_site=True)
    single_site = evidence.loc[
        evidence.loc[:, "peptide_row_id"] == "pep_single", :
    ].copy(deep=True)
    single_site.loc[:, "site_sequence"] = " aaAsaaa "

    resolved = PeptideEvidenceDatasetResolver().run(
        evidence=PeptideEvidenceTable(
            frame=single_site.reset_index(drop=True),
            sample_intensity_columns=("sample_a", "sample_b"),
        ),
        multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
        input_intensity_scale="linear",
    )

    assert resolved.site_metadata.loc["AKT1;S473;", "site_sequence"] == "AAASAAA"
    payload = resolved.summary.to_payload()
    assert payload["site_sequence_column_present"] is True
    assert int(payload["provided_site_sequence_count"]) == 1
    assert int(payload["accepted_site_sequence_count"]) == 1
    assert int(payload["rejected_site_sequence_count"]) == 0
    assert (
        payload["site_sequence_policy"]
        == DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR
    )


def test_peptide_evidence_resolution_records_absent_sequence_context() -> None:
    evidence = _peptide_evidence_frame(include_single_site=True)
    single_site = evidence.loc[
        evidence.loc[:, "peptide_row_id"] == "pep_single", :
    ].copy(deep=True)
    single_site = single_site.drop(columns=["site_sequence"])

    resolved = PeptideEvidenceDatasetResolver().run(
        evidence=PeptideEvidenceTable(
            frame=single_site.reset_index(drop=True),
            sample_intensity_columns=("sample_a", "sample_b"),
        ),
        multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
        input_intensity_scale="linear",
    )

    payload = resolved.summary.to_payload()
    assert payload["site_sequence_column_present"] is False
    assert int(payload["provided_site_sequence_count"]) == 0
    assert int(payload["accepted_site_sequence_count"]) == 0
    assert int(payload["rejected_site_sequence_count"]) == 0
    assert (
        payload["site_sequence_policy"]
        == DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR
    )


def test_reject_policy_fails_on_ambiguous_peptide_evidence() -> None:
    request = DatasetBuildRequest(
        site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
        peptide_evidence=_peptide_evidence_frame(include_single_site=False),
        peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
        multi_site_policy=DATASET_MULTI_SITE_POLICY_REJECT,
        input_intensity_scale="linear",
        organism=Organism.HUMAN,
    )
    with pytest.raises(
        PhosPyInputError,
        match="multi-site observation cannot be resolved under policy='error'",
    ):
        AnalysisReadyDatasetBuilder().run(request)


def test_exclude_policy_records_exclusions_in_report_and_provenance() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=_peptide_evidence_frame(),
            peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
            input_intensity_scale="linear",
            organism=Organism.HUMAN,
        )
    )
    assert built.phospho.index.name == "site_key"
    assert _display_ids(built) == ["AKT1;S473;"]
    assert built.preprocessing_report is not None
    resolution_rows = built.preprocessing_report.operations.loc[
        built.preprocessing_report.operations.loc[:, "stage"]
        == "peptide_evidence_resolution"
    ]
    assert int(resolution_rows.shape[0]) == 1
    parameters = resolution_rows.iloc[0]["parameters"]
    assert isinstance(parameters, dict)
    assert int(parameters["peptide_observations_received"]) == 2
    assert int(parameters["ambiguous_observations"]) == 1
    assert int(parameters["excluded_observations"]) == 1
    assert int(parameters["split_observations"]) == 0
    assert built.provenance is not None
    payload = built.provenance.workflow_parameters["peptide_evidence_resolution"]
    assert isinstance(payload, Mapping)
    assert (
        payload["multi_site_policy"]
        == DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING
    )


def test_keep_joint_policy_rejects_joint_ambiguous_site_representation() -> None:
    with pytest.raises(PhosPyInputError, match="strict residue/position"):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
                peptide_evidence=_peptide_evidence_frame(include_single_site=False),
                peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
                multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
                input_intensity_scale="linear",
                allow_opaque_site_values=True,
                organism=Organism.HUMAN,
            )
        )


def test_keep_joint_policy_without_opaque_opt_in_fails_dataset_validation() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="strict residue/position",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
                peptide_evidence=_peptide_evidence_frame(include_single_site=False),
                peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
                multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
                input_intensity_scale="linear",
                organism=Organism.HUMAN,
            )
        )


def test_split_policy_applies_deterministic_equal_split() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=_peptide_evidence_frame(include_single_site=False),
            peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
            input_intensity_scale="linear",
            organism=Organism.HUMAN,
        )
    )
    assert set(_display_ids(built)) == {"MAPK1;S10;", "MAPK1;T12;"}
    mapk1_s10 = _site_key_for_display_id(built, "MAPK1;S10;")
    mapk1_t12 = _site_key_for_display_id(built, "MAPK1;T12;")
    assert float(built.phospho.loc[mapk1_s10, "sample_a"]) == pytest.approx(5.0)
    assert float(built.phospho.loc[mapk1_t12, "sample_a"]) == pytest.approx(5.0)
    assert built.provenance is not None
    payload = built.provenance.workflow_parameters["peptide_evidence_resolution"]
    assert isinstance(payload, Mapping)
    assert int(payload["split_observations"]) == 1


def test_split_policy_mixed_ambiguous_and_unambiguous_rows_is_deterministic() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=_peptide_evidence_frame(include_single_site=True),
            peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
            input_intensity_scale="linear",
            organism=Organism.HUMAN,
        )
    )
    assert set(_display_ids(built)) == {
        "MAPK1;S10;",
        "MAPK1;T12;",
        "AKT1;S473;",
    }
    mapk1_s10 = _site_key_for_display_id(built, "MAPK1;S10;")
    mapk1_t12 = _site_key_for_display_id(built, "MAPK1;T12;")
    akt1_s473 = _site_key_for_display_id(built, "AKT1;S473;")
    assert float(built.phospho.loc[mapk1_s10, "sample_a"]) == pytest.approx(5.0)
    assert float(built.phospho.loc[mapk1_t12, "sample_a"]) == pytest.approx(5.0)
    assert float(built.phospho.loc[akt1_s473, "sample_a"]) == pytest.approx(7.0)
    assert float(built.phospho.loc[mapk1_s10, "sample_b"]) == pytest.approx(6.0)
    assert float(built.phospho.loc[mapk1_t12, "sample_b"]) == pytest.approx(6.0)
    assert float(built.phospho.loc[akt1_s473, "sample_b"]) == pytest.approx(9.0)


def test_multiple_peptides_mapping_to_one_site_are_mean_aggregated() -> None:
    evidence = pd.DataFrame(
        [
            {
                "peptide_row_id": "pep_1",
                "site_id": "MAPK1;S10;",
                "unique_feature_id": "feat_1",
                "gene_symbol": "MAPK1",
                "protein_accession": "P28482",
                "site_string": "S10",
                "sample_a": 10.0,
                "sample_b": 20.0,
                "peptide_sequence": "AAAAA",
                "modified_peptide_sequence": "AA[+80]AAA",
                "multi_site": False,
                "provenance_source": "maxquant",
                "site_sequence": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                "localisation_confidence": 0.95,
            },
            {
                "peptide_row_id": "pep_2",
                "site_id": "MAPK1;S10;",
                "unique_feature_id": "feat_2",
                "gene_symbol": "MAPK1",
                "protein_accession": "P28482",
                "site_string": "S10",
                "sample_a": 14.0,
                "sample_b": 26.0,
                "peptide_sequence": "BBBBB",
                "modified_peptide_sequence": "BB[+80]BBB",
                "multi_site": False,
                "provenance_source": "maxquant",
                "site_sequence": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                "localisation_confidence": 0.92,
            },
        ]
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=evidence,
            peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
            input_intensity_scale="linear",
            organism=Organism.HUMAN,
        )
    )
    mapk1_s10 = _site_key_for_display_id(built, "MAPK1;S10;")
    assert list(built.phospho.index.astype(str)) == [mapk1_s10]
    assert float(built.phospho.loc[mapk1_s10, "sample_a"]) == pytest.approx(12.0)
    assert float(built.phospho.loc[mapk1_s10, "sample_b"]) == pytest.approx(23.0)


def test_single_peptide_split_allocates_fraction_before_site_mean() -> None:
    resolved = _run_policy_resolution(
        [
            _resolution_policy_evidence_row(
                peptide_row_id="pep_1",
                sample_a=10.0,
                site_string="S10,T12",
                multi_site=True,
            )
        ],
        split=True,
    )

    assert float(resolved.phospho.loc["MAPK1;S10;", "sample_a"]) == pytest.approx(5.0)
    assert float(resolved.phospho.loc["MAPK1;T12;", "sample_a"]) == pytest.approx(5.0)


def test_two_split_rows_use_arithmetic_mean_of_allocated_values() -> None:
    resolved = _run_policy_resolution(
        [
            _resolution_policy_evidence_row(
                peptide_row_id="pep_1",
                sample_a=10.0,
                site_string="S10,T12",
                peptide_sequence="AAASAAA",
                multi_site=True,
            ),
            _resolution_policy_evidence_row(
                peptide_row_id="pep_2",
                sample_a=30.0,
                site_string="S10,T12",
                peptide_sequence="CCCSCCC",
                multi_site=True,
            ),
        ],
        split=True,
    )

    assert float(resolved.phospho.loc["MAPK1;S10;", "sample_a"]) == pytest.approx(10.0)
    assert float(resolved.phospho.loc["MAPK1;T12;", "sample_a"]) == pytest.approx(10.0)


def test_site_summarisation_is_not_conventional_normalized_weighted_mean() -> None:
    mapping = pd.DataFrame(
        {
            "peptide_row_id": ["pep_1", "pep_1", "pep_2", "pep_2"],
            "site_id": ["MAPK1;S10;", "MAPK1;T12;", "MAPK1;S10;", "MAPK1;T12;"],
            "mapping_weight": [0.2, 0.8, 0.8, 0.2],
        }
    )
    resolved = _run_policy_resolution(
        [
            _resolution_policy_evidence_row(
                peptide_row_id="pep_1",
                sample_a=10.0,
                site_string="S10,T12",
                peptide_sequence="AAASAAA",
                multi_site=True,
            ),
            _resolution_policy_evidence_row(
                peptide_row_id="pep_2",
                sample_a=30.0,
                site_string="S10,T12",
                peptide_sequence="CCCSCCC",
                multi_site=True,
            ),
        ],
        split=True,
        site_mapping=mapping,
    )

    current_policy_value = float(resolved.phospho.loc["MAPK1;S10;", "sample_a"])
    conventional_normalized_weighted_mean = ((0.2 * 10.0) + (0.8 * 30.0)) / (0.2 + 0.8)
    assert current_policy_value == pytest.approx(13.0)
    assert current_policy_value != pytest.approx(conventional_normalized_weighted_mean)


def test_unequal_explicit_fractions_use_mean_of_allocated_values() -> None:
    mapping = pd.DataFrame(
        {
            "peptide_row_id": ["pep_1", "pep_1", "pep_2", "pep_2"],
            "site_id": ["MAPK1;S10;", "MAPK1;T12;", "MAPK1;S10;", "MAPK1;T12;"],
            "mapping_weight": [0.25, 0.75, 0.75, 0.25],
        }
    )
    resolved = _run_policy_resolution(
        [
            _resolution_policy_evidence_row(
                peptide_row_id="pep_1",
                sample_a=10.0,
                site_string="S10,T12",
                peptide_sequence="AAASAAA",
                multi_site=True,
            ),
            _resolution_policy_evidence_row(
                peptide_row_id="pep_2",
                sample_a=20.0,
                site_string="S10,T12",
                peptide_sequence="CCCSCCC",
                multi_site=True,
            ),
        ],
        split=True,
        site_mapping=mapping,
    )

    assert float(resolved.phospho.loc["MAPK1;S10;", "sample_a"]) == pytest.approx(
        (2.5 + 15.0) / 2.0
    )
    assert float(resolved.phospho.loc["MAPK1;T12;", "sample_a"]) == pytest.approx(
        (7.5 + 5.0) / 2.0
    )


def test_mapping_weight_sum_tolerance_accepts_near_unit_fraction_total() -> None:
    mapping = pd.DataFrame(
        {
            "peptide_row_id": ["pep_1", "pep_1"],
            "site_id": ["MAPK1;S10;", "MAPK1;T12;"],
            "mapping_weight": [0.7000004, 0.2999999],
        }
    )
    resolved = _run_policy_resolution(
        [
            _resolution_policy_evidence_row(
                peptide_row_id="pep_1",
                sample_a=10.0,
                site_string="S10,T12",
                multi_site=True,
            )
        ],
        split=True,
        site_mapping=mapping,
    )

    assert float(resolved.phospho.loc["MAPK1;S10;", "sample_a"]) == pytest.approx(
        7.000004
    )


def test_allocated_signal_summarisation_is_invariant_to_row_permutation() -> None:
    rows = [
        _resolution_policy_evidence_row(
            peptide_row_id="pep_1",
            sample_a=10.0,
            site_string="S10,T12",
            peptide_sequence="AAASAAA",
            multi_site=True,
        ),
        _resolution_policy_evidence_row(
            peptide_row_id="pep_2",
            sample_a=30.0,
            site_string="S10,T12",
            peptide_sequence="CCCSCCC",
            multi_site=True,
        ),
        _resolution_policy_evidence_row(
            peptide_row_id="pep_3",
            sample_a=20.0,
            peptide_sequence="GGGSGGG",
        ),
    ]
    expected = _run_policy_resolution(rows, split=True).phospho

    for row_order in permutations(rows):
        resolved = _run_policy_resolution(list(row_order), split=True)
        pd.testing.assert_frame_equal(resolved.phospho, expected)


def test_duplicate_evidence_rows_affect_mean_under_retained_duplicate_policy() -> None:
    unduplicated = _run_policy_resolution(
        [
            _resolution_policy_evidence_row(
                peptide_row_id="pep_1",
                sample_a=10.0,
                peptide_sequence="AAASAAA",
            ),
            _resolution_policy_evidence_row(
                peptide_row_id="pep_2",
                sample_a=30.0,
                peptide_sequence="CCCSCCC",
            ),
        ]
    )
    duplicated = _run_policy_resolution(
        [
            _resolution_policy_evidence_row(
                peptide_row_id="pep_1",
                sample_a=10.0,
                peptide_sequence="AAASAAA",
            ),
            _resolution_policy_evidence_row(
                peptide_row_id="pep_1_duplicate",
                sample_a=10.0,
                peptide_sequence="AAASAAA",
            ),
            _resolution_policy_evidence_row(
                peptide_row_id="pep_2",
                sample_a=30.0,
                peptide_sequence="CCCSCCC",
            ),
        ]
    )

    assert float(unduplicated.phospho.loc["MAPK1;S10;", "sample_a"]) == pytest.approx(
        20.0
    )
    assert float(duplicated.phospho.loc["MAPK1;S10;", "sample_a"]) == pytest.approx(
        (10.0 + 10.0 + 30.0) / 3.0
    )
    assert duplicated.summary.duplicate_evidence_policy == (
        DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS
    )
    assert duplicated.summary.duplicate_peptide_rows == 2


def test_mixed_ambiguous_and_unambiguous_rows_share_allocated_signal_mean() -> None:
    resolved = _run_policy_resolution(
        [
            _resolution_policy_evidence_row(
                peptide_row_id="pep_split",
                sample_a=10.0,
                site_string="S10,T12",
                peptide_sequence="AAASAAA",
                multi_site=True,
            ),
            _resolution_policy_evidence_row(
                peptide_row_id="pep_single",
                sample_a=20.0,
                peptide_sequence="CCCSCCC",
            ),
        ],
        split=True,
    )

    assert float(resolved.phospho.loc["MAPK1;S10;", "sample_a"]) == pytest.approx(
        (5.0 + 20.0) / 2.0
    )
    assert float(resolved.phospho.loc["MAPK1;T12;", "sample_a"]) == pytest.approx(5.0)
    assert resolved.summary.mixed_ambiguity_policy == (
        DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS
    )


def test_localisation_aggregation_uses_mean_of_finite_reported_values() -> None:
    resolved = _run_policy_resolution(
        [
            _resolution_policy_evidence_row(
                peptide_row_id="pep_1",
                sample_a=10.0,
                peptide_sequence="AAASAAA",
                localisation_confidence=0.2,
            ),
            _resolution_policy_evidence_row(
                peptide_row_id="pep_2",
                sample_a=20.0,
                peptide_sequence="CCCSCCC",
                localisation_confidence=pd.NA,
            ),
            _resolution_policy_evidence_row(
                peptide_row_id="pep_3",
                sample_a=30.0,
                peptide_sequence="GGGSGGG",
                localisation_confidence=1.0,
            ),
        ]
    )

    assert float(
        resolved.site_metadata.loc["MAPK1;S10;", "localisation_confidence"]
    ) == pytest.approx(0.6)
    assert resolved.summary.localisation_aggregation_policy == (
        DATASET_PEPTIDE_LOCALISATION_AGGREGATION_POLICY_ARITHMETIC_MEAN_OF_FINITE_VALUES
    )


def test_duplicate_peptide_rows_are_retained_as_independent_observations() -> None:
    evidence = _peptide_evidence_frame(include_single_site=False)
    evidence = pd.concat([evidence, evidence.copy(deep=True)], ignore_index=True)
    evidence.loc[0, "peptide_row_id"] = "pep_joint_a"
    evidence.loc[1, "peptide_row_id"] = "pep_joint_b"
    evidence.loc[:, "sample_a"] = [10.0, 30.0]
    evidence.loc[:, "sample_b"] = [12.0, 28.0]

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=evidence,
            peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
            input_intensity_scale="linear",
            organism=Organism.HUMAN,
        )
    )
    mapk1_s10 = _site_key_for_display_id(built, "MAPK1;S10;")
    mapk1_t12 = _site_key_for_display_id(built, "MAPK1;T12;")
    assert float(built.phospho.loc[mapk1_s10, "sample_a"]) == pytest.approx(10.0)
    assert float(built.phospho.loc[mapk1_t12, "sample_a"]) == pytest.approx(10.0)
    assert float(built.phospho.loc[mapk1_s10, "sample_b"]) == pytest.approx(10.0)
    assert float(built.phospho.loc[mapk1_t12, "sample_b"]) == pytest.approx(10.0)
    assert built.provenance is not None
    payload = built.provenance.workflow_parameters["peptide_evidence_resolution"]
    assert isinstance(payload, Mapping)
    assert (
        payload["duplicate_peptide_policy"]
        == DATASET_PEPTIDE_DUPLICATE_POLICY_RETAIN_ALL_ROWS
    )
    assert int(payload["duplicate_peptide_rows"]) == 2


def test_explicit_mapping_weights_are_applied_deterministically() -> None:
    evidence_frame = _peptide_evidence_frame(include_single_site=False)
    mapping = pd.DataFrame(
        {
            "peptide_row_id": ["pep_joint", "pep_joint"],
            "site_id": ["MAPK1;S10;", "MAPK1;T12;"],
            "mapping_weight": [0.7, 0.3],
            "mapping_uncertainty": [True, True],
        }
    )
    resolved = PeptideEvidenceDatasetResolver().run(
        evidence=PeptideEvidenceTable(
            frame=evidence_frame,
            sample_intensity_columns=("sample_a", "sample_b"),
            site_mapping=mapping,
        ),
        multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
        input_intensity_scale="linear",
    )
    assert float(resolved.phospho.loc["MAPK1;S10;", "sample_a"]) == pytest.approx(7.0)
    assert float(resolved.phospho.loc["MAPK1;T12;", "sample_a"]) == pytest.approx(3.0)
    payload = resolved.summary.to_payload()
    assert (
        payload["mapping_weight_source"]
        == DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT
    )


def test_log2_peptide_evidence_rejects_explicit_fractional_mapping_weights() -> None:
    evidence_frame = _peptide_evidence_frame(include_single_site=False)
    mapping = pd.DataFrame(
        {
            "peptide_row_id": ["pep_joint", "pep_joint"],
            "site_id": ["MAPK1;S10;", "MAPK1;T12;"],
            "mapping_weight": [0.5, 0.5],
        }
    )

    with pytest.raises(PhosPyInputError) as exc_info:
        PeptideEvidenceDatasetResolver().run(
            evidence=PeptideEvidenceTable(
                frame=evidence_frame,
                sample_intensity_columns=("sample_a", "sample_b"),
                site_mapping=mapping,
            ),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
            input_intensity_scale="log2",
        )

    _assert_log2_fractional_allocation_message(str(exc_info.value))


def test_log2_peptide_evidence_rejects_former_ten_to_five_corruption_path() -> None:
    evidence_frame = _peptide_evidence_frame(include_single_site=False)
    evidence_frame.loc[:, "sample_a"] = 10.0

    with pytest.raises(PhosPyInputError) as exc_info:
        PeptideEvidenceDatasetResolver().run(
            evidence=PeptideEvidenceTable(
                frame=evidence_frame,
                sample_intensity_columns=("sample_a",),
                multi_site_handling_config=(
                    build_multi_site_handling_config_for_dataset_policy(
                        multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT
                    )
                ),
            ),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
            input_intensity_scale="log2",
        )

    _assert_log2_fractional_allocation_message(str(exc_info.value))


def test_mapping_weights_must_sum_to_one_per_peptide_row() -> None:
    evidence_frame = _peptide_evidence_frame(include_single_site=False)
    mapping = pd.DataFrame(
        {
            "peptide_row_id": ["pep_joint", "pep_joint"],
            "site_id": ["MAPK1;S10;", "MAPK1;T12;"],
            "mapping_weight": [0.7, 0.4],
        }
    )
    with pytest.raises(PhosPyInputError, match="must sum to 1.0 per peptide_row_id"):
        PeptideEvidenceDatasetResolver().run(
            evidence=PeptideEvidenceTable(
                frame=evidence_frame,
                sample_intensity_columns=("sample_a", "sample_b"),
                site_mapping=mapping,
            ),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
            input_intensity_scale="linear",
        )


def test_new_resolution_summary_serializes_explicit_policy_identifiers_in_order() -> (
    None
):
    resolved = _run_policy_resolution(
        [
            _resolution_policy_evidence_row(
                peptide_row_id="pep_1",
                sample_a=10.0,
                site_string="S10,T12",
                multi_site=True,
            )
        ],
        split=True,
    )

    payload = resolved.summary.to_payload()
    ordered_policy_keys = [
        "mapping_weight_source_policy",
        "mapping_weight_normalization_policy",
        "signal_allocation_policy",
        "site_summarisation_policy",
        "duplicate_evidence_policy",
        "mixed_ambiguity_policy",
        "localisation_aggregation_policy",
    ]
    assert [key for key in payload if key in ordered_policy_keys] == (
        ordered_policy_keys
    )
    assert payload["mapping_weight_source_policy"] == (
        DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL
    )
    assert payload["mapping_weight_normalization_policy"] == (
        DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW
    )
    assert payload["signal_allocation_policy"] == (
        DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION
    )
    assert payload["site_summarisation_policy"] == (
        DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS
    )
    assert payload["duplicate_evidence_policy"] == (
        DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS
    )
    assert payload["mixed_ambiguity_policy"] == (
        DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS
    )
    assert payload["localisation_aggregation_policy"] == (
        DATASET_PEPTIDE_LOCALISATION_AGGREGATION_POLICY_ARITHMETIC_MEAN_OF_FINITE_VALUES
    )
    assert (
        payload["aggregation_policy"]
        == DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS
    )
    assert (
        payload["aggregation_policy"]
        != DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN
    )


def test_legacy_mapping_weighted_mean_provenance_reconstructs_policy_tuple() -> None:
    legacy_resolution_payload: dict[str, object] = {
        "input_mode": DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
        "multi_site_policy": DATASET_MULTI_SITE_POLICY_SPLIT,
        "peptide_observations_received": 1,
        "unique_site_ids_produced": 2,
        "ambiguous_observations": 1,
        "excluded_observations": 0,
        "split_observations": 1,
        "aggregation_policy": (
            DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN
        ),
        "aggregation_formula": (
            "site_intensity = mean(per_peptide_intensity * mapping_weight)"
        ),
        "mapping_weight_source": DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT,
        "mapping_weight_normalisation": "sum_to_one_per_peptide_row",
        "duplicate_peptide_policy": (
            "retain_all_peptide_rows_as_independent_observations"
        ),
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

    summary = PeptideEvidenceResolutionSummary.from_payload(legacy_resolution_payload)
    assert summary.signal_allocation_policy == (
        DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION
    )
    assert summary.site_summarisation_policy == (
        DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS
    )
    assert summary.aggregation_policy == (
        DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS
    )

    restored = from_payload(
        {
            "environment": {
                "schema_version": 2,
                "package_name": "phospy",
                "package_version": "0",
                "python_version": "3.13",
                "dependency_versions": {},
                "platform": {},
                "blas_lapack": {},
                "thread_environment": {},
                "timezone": None,
                "locale": {},
                "constraints_fingerprint": {},
            },
            "input_tables": [],
            "preprocessing_stages": [],
            "reference": None,
            "reference_context": None,
            "workflow_name": "dataset_builder",
            "workflow_parameters": {
                "peptide_evidence_resolution": legacy_resolution_payload
            },
            "random_state": None,
            "random_seed_policy": None,
            "output_tables": [],
            "scientific_policies": [],
        }
    )
    restored_resolution = restored.workflow_parameters["peptide_evidence_resolution"]
    assert isinstance(restored_resolution, Mapping)
    assert restored_resolution["mapping_weight_source_policy"] == (
        DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL
    )
    assert restored_resolution["mapping_weight_normalization_policy"] == (
        DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW
    )
    assert restored_resolution["signal_allocation_policy"] == (
        DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION
    )
    assert restored_resolution["site_summarisation_policy"] == (
        DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS
    )
    assert restored_resolution["duplicate_evidence_policy"] == (
        DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS
    )
    assert restored_resolution["mixed_ambiguity_policy"] == (
        DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS
    )
    assert restored_resolution["localisation_aggregation_policy"] == (
        DATASET_PEPTIDE_LOCALISATION_AGGREGATION_POLICY_ARITHMETIC_MEAN_OF_FINITE_VALUES
    )
    assert restored_resolution["aggregation_policy"] == (
        DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS
    )


def test_peptide_evidence_resolution_provenance_records_aggregation_semantics() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=_peptide_evidence_frame(include_single_site=True),
            peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
            input_intensity_scale="linear",
            organism=Organism.HUMAN,
        )
    )
    assert built.provenance is not None
    payload = built.provenance.workflow_parameters["peptide_evidence_resolution"]
    assert isinstance(payload, Mapping)
    assert payload["mapping_weight_source_policy"] == (
        DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL
    )
    assert (
        payload["mapping_weight_normalization_policy"]
        == DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW
    )
    assert payload["signal_allocation_policy"] == (
        DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION
    )
    assert payload["site_summarisation_policy"] == (
        DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS
    )
    assert payload["duplicate_evidence_policy"] == (
        DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS
    )
    assert payload["mixed_ambiguity_policy"] == (
        DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS
    )
    assert payload["localisation_aggregation_policy"] == (
        DATASET_PEPTIDE_LOCALISATION_AGGREGATION_POLICY_ARITHMETIC_MEAN_OF_FINITE_VALUES
    )
    assert (
        payload["aggregation_policy"]
        == DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS
    )
    assert (
        payload["mapping_weight_normalisation"]
        == DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE
    )
    assert payload["site_sequence_column_present"] is True
    assert int(payload["provided_site_sequence_count"]) == 2
    assert int(payload["accepted_site_sequence_count"]) == 3
    assert int(payload["rejected_site_sequence_count"]) == 0
    assert (
        payload["site_sequence_policy"]
        == DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR
    )
    assert built.preprocessing_report is not None
    resolution_row = built.preprocessing_report.operations.loc[
        built.preprocessing_report.operations.loc[:, "stage"]
        == "peptide_evidence_resolution"
    ].iloc[0]
    parameters = resolution_row["parameters"]
    assert isinstance(parameters, dict)
    assert parameters["signal_allocation_policy"] == (
        DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION
    )
    assert parameters["site_summarisation_policy"] == (
        DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS
    )
    assert (
        parameters["aggregation_policy"]
        == DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS
    )
    assert (
        parameters["site_sequence_policy"]
        == DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR
    )
    assert int(parameters["provided_site_sequence_count"]) == 2
    assert int(parameters["accepted_site_sequence_count"]) == 3

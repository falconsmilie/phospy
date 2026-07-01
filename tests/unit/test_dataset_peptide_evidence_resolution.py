from __future__ import annotations

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
from phospy.science.datasets.builders.validator import DatasetBuildRequestValidator
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.evidence import (
    DATASET_PEPTIDE_DUPLICATE_POLICY_RETAIN_ALL_ROWS,
    DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE,
    DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT,
    DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR,
    DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN,
    PeptideEvidenceDatasetResolver,
    PeptideEvidenceTable,
)


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


def _assert_site_sequence_center_mismatch_message(message: str) -> None:
    lower_message = message.lower()
    assert "site_sequence" in message
    assert "centre" in lower_message or "center" in lower_message
    assert "expected='S'" in message
    assert "observed='T'" in message
    assert "AKT1;S473;" in message


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


def test_peptide_evidence_rejects_site_sequence_center_residue_mismatch() -> None:
    with pytest.raises(PhosPyInputError) as exc_info:
        PeptideEvidenceDatasetResolver().run(
            evidence=PeptideEvidenceTable(
                frame=_peptide_evidence_center_residue_mismatch_frame(),
                sample_intensity_columns=("sample_a", "sample_b"),
            ),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
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
    assert isinstance(payload, dict)
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
    assert isinstance(payload, dict)
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
    assert isinstance(payload, dict)
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
    )
    assert float(resolved.phospho.loc["MAPK1;S10;", "sample_a"]) == pytest.approx(7.0)
    assert float(resolved.phospho.loc["MAPK1;T12;", "sample_a"]) == pytest.approx(3.0)
    payload = resolved.summary.to_payload()
    assert (
        payload["mapping_weight_source"]
        == DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT
    )


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
    assert isinstance(payload, dict)
    assert (
        payload["aggregation_policy"]
        == DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN
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
    assert (
        parameters["aggregation_policy"]
        == DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN
    )
    assert (
        parameters["site_sequence_policy"]
        == DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR
    )
    assert int(parameters["provided_site_sequence_count"]) == 2
    assert int(parameters["accepted_site_sequence_count"]) == 3

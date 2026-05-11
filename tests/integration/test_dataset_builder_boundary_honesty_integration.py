from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    Organism,
    PhosPyInputError,
)
from phospy.api.configs import DATASET_SITE_MATRIX_MISSING_DATA_POLICIES
from phospy.errors.references import UnsupportedOrganismError

pytestmark = pytest.mark.integration


class _InterpreterSentinel:
    def __init__(self) -> None:
        self.called = False

    def run(self, request: object) -> object:
        self.called = True
        raise AssertionError(
            "interpreter must not run when invalid site_matrix missing-data policy "
            "is rejected at validator boundary"
        )


@pytest.mark.parametrize(
    "missing_data_policy",
    tuple(sorted(DATASET_SITE_MATRIX_MISSING_DATA_POLICIES)),
)
def test_builder_supports_all_publicly_advertised_site_matrix_missing_data_modes_end_to_end(
    missing_data_policy: str,
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    missing_data_policy=missing_data_policy,
                )
            ),
        )
    )

    assert built.phospho.index.tolist() == ["GSK3B;S9;", "MAPK14;Y182;"]
    assert built.phospho.isna().to_numpy().sum() == 0
    assert built.preprocessing_report is not None
    assert built.preprocessing_report.row_audit is not None
    dropped = built.preprocessing_report.row_audit.loc[
        (built.preprocessing_report.row_audit.loc[:, "stage"] == "site_matrix")
        & (built.preprocessing_report.row_audit.loc[:, "action"] == "dropped")
    ]
    assert dropped.empty


@pytest.mark.parametrize(
    "missing_data_policy",
    ("retain_missing", "require_min_observed_values"),
)
def test_builder_rejects_dead_end_site_matrix_missing_data_modes_before_interpretation(
    missing_data_policy: str,
) -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [float("nan")]},
        index=pd.Index(["row_a"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
        },
        index=phospho.index.copy(),
    )
    interpreter = _InterpreterSentinel()
    builder = AnalysisReadyDatasetBuilder(interpreter=interpreter)

    with pytest.raises(
        PhosPyInputError,
        match=(
            f"missing_data_policy='{missing_data_policy}' is not supported for strict "
            "AnalysisReadyPhosphoDataset construction"
        ),
    ):
        builder.run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(
                        policy="build_from_metadata",
                        missing_data_policy=missing_data_policy,  # type: ignore[arg-type]
                    )
                ),
            )
        )

    assert interpreter.called is False


def test_builder_site_sequence_mixed_support_keeps_resolvable_rows_and_excludes_unresolved_rows() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0, 4.0],
            "sample_b": [1.5, 2.5, 3.5, 4.5],
        },
        index=pd.Index(["row_a", "row_b", "row_c", "row_d"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "FAKE1", "FAKE2"],
            "site": ["Y182", "S9", "S1", "T2"],
            "site_sequence": [pd.NA, "SEQ_MANUAL", pd.NA, "   "],
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            ),
        )
    )

    assert built.phospho.index.tolist() == ["GSK3B;S9;", "MAPK14;Y182;"]
    assert built.site_metadata.index.tolist() == ["GSK3B;S9;", "MAPK14;Y182;"]
    assert built.site_metadata.loc["GSK3B;S9;", "site_sequence"] == "SEQ_MANUAL"
    assert isinstance(built.site_metadata.loc["MAPK14;Y182;", "site_sequence"], str)
    assert built.preprocessing_report is not None
    assert built.preprocessing_report.row_audit is not None
    dropped = built.preprocessing_report.row_audit.loc[
        (built.preprocessing_report.row_audit.loc[:, "stage"] == "site_matrix")
        & (built.preprocessing_report.row_audit.loc[:, "action"] == "dropped")
    ]
    assert set(dropped.loc[:, "source_row_id"].astype(str)) == {"row_c", "row_d"}
    assert "missing or blank" in str(dropped.iloc[0]["reason"])


def test_builder_rejects_missing_sample_label_before_stringification() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0],
            "sample_b": [2.0],
        },
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
        },
        index=phospho.index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {"comparison_group": ["group_1", "group_2"]},
        index=pd.Index(["sample_a", pd.NA], name="sample_id"),
    )

    with pytest.raises(
        PhosPyInputError,
        match="sample_metadata.index must not contain missing labels",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                organism=Organism.RAT,
            )
        )


def test_builder_provenance_exposes_site_identifier_normalisation_records() -> None:
    phospho = pd.DataFrame(
        {
            " sample_a ": [1.0],
            " sample_b ": [2.0],
        },
        index=pd.Index([" mapk14 ; y182 "], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["mapk14"],
            "site": ["y182"],
            "site_sequence": ["SEQ_A"],
        },
        index=pd.Index([" mapk14 ; y182 "], name="site_id"),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )

    assert built.provenance is not None
    payload = built.provenance.workflow_parameters.get("site_identifier_normalisation")
    assert isinstance(payload, dict)
    assert payload["changed_identifier_count"] >= 2


def test_builder_succeeds_with_provided_site_sequence_at_analysis_ready_boundary() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(phospho=phospho, site_metadata=site_metadata)
    )

    assert built.site_metadata.loc["MAPK14;Y182;", "site_sequence"] == "SEQ_A"
    assert built.provenance is not None
    derivation = built.provenance.workflow_parameters.get("site_sequence_derivation")
    assert isinstance(derivation, dict)
    assert derivation["provided_sequence_count"] == 1
    assert derivation["derived_sequence_count"] == 0
    assert derivation["unresolved_sequence_count"] == 0


def test_builder_derives_missing_site_sequence_before_analysis_ready_construction() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [1.5, 2.5]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )

    assert built.site_metadata.loc[:, "site_sequence"].isna().sum() == 0
    assert built.provenance is not None
    derivation = built.provenance.workflow_parameters.get("site_sequence_derivation")
    assert isinstance(derivation, dict)
    assert derivation["provided_sequence_count"] == 0
    assert derivation["derived_sequence_count"] == 2
    assert derivation["unresolved_sequence_count"] == 0
    assert derivation["reference_support"] == "available"
    assert derivation["reference_bundle_id"] == "l6_native"
    assert isinstance(derivation["reference_manifest"], dict)
    assert derivation["reference_manifest"]["bundle_id"] == "l6_native"
    assert str(derivation["reference_source"]).startswith("bundled_reference:rat/")


def test_builder_fails_when_site_sequence_cannot_be_derived_before_dataset_construction() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [1.5]},
        index=pd.Index(["FAKE1;S1;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["FAKE1"],
            "site": ["S1"],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="cannot construct AnalysisReadyPhosphoDataset",
    ) as caught:
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
            )
        )

    message = str(caught.value)
    assert "site_id='FAKE1;S1;'" in message
    assert "gene_symbol='FAKE1'" in message
    assert "site='S1'" in message
    assert "failure_category='missing_reference_support'" in message


def test_builder_fails_clearly_for_unsupported_organism_when_site_sequences_need_derivation() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [1.5]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        UnsupportedOrganismError, match="supported bundled organisms: rat"
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.HUMAN,
            )
        )

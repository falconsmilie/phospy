from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    Organism,
    ReferencePreset,
    SimpleKinaseWorkflow,
    SimpleKinaseWorkflowRequest,
)
from phospy.errors import (
    DatasetValidationError,
    UnsupportedInputFormatError,
    UnsupportedOrganismError,
)
from phospy.references.resolution import ReferenceResolver
from phospy.transformations.models import TransformationState

ROOT = Path(__file__).resolve().parents[2]
RAT_L6_PHOSPHO = (
    ROOT / "tests_legacy" / "fixtures" / "r_reference_l6" / "l6_phospho_matrix.csv"
)
RAT_L6_EXPECTED_PROFILE = (
    ROOT / "tests_legacy" / "fixtures" / "r_reference_l6" / "native_profile_scores.csv"
)
RAT_L6_SITE_SEQUENCES = (
    ROOT
    / "src"
    / "phospy"
    / "data"
    / "reference_bundles"
    / "rat"
    / "l6_native"
    / "site_sequences.csv"
)


@lru_cache(maxsize=1)
def _load_rat_l6_phospho() -> pd.DataFrame:
    return pd.read_csv(RAT_L6_PHOSPHO, index_col=0)


@lru_cache(maxsize=1)
def _load_rat_l6_sequence_table() -> pd.Series:
    sequence_frame = pd.read_csv(RAT_L6_SITE_SEQUENCES)
    return sequence_frame.set_index("site_id").loc[:, "centralized_sequence"]


@lru_cache(maxsize=1)
def _load_expected_profile_scores() -> pd.DataFrame:
    return pd.read_csv(RAT_L6_EXPECTED_PROFILE, index_col=0)


def _site_metadata_for(phospho: pd.DataFrame) -> pd.DataFrame:
    split = phospho.index.to_series().astype(str).str.split(";", expand=True)
    site_sequences = _load_rat_l6_sequence_table().reindex(phospho.index)
    if site_sequences.isna().any():
        missing = int(site_sequences.isna().sum())
        raise AssertionError(
            f"fixture missing site sequences for {missing} phosphosites"
        )
    return pd.DataFrame(
        {
            "gene_symbol": split.loc[:, 0].values,
            "site": split.loc[:, 1].values,
            "site_sequence": site_sequences.values,
        },
        index=phospho.index.copy(),
    )


def _build_dataset(*, n_sites: int | None = 220) -> object:
    phospho = _load_rat_l6_phospho().copy(deep=True)
    if n_sites is not None:
        phospho = phospho.head(n_sites)
    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=_site_metadata_for(phospho),
        organism=Organism.RAT,
    )
    return AnalysisReadyDatasetBuilder().run(request)


def test_builder_happy_path_supports_dataframe_route() -> None:
    phospho = _load_rat_l6_phospho().head(32).copy(deep=True)
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata_for(phospho),
            organism=Organism.RAT,
        )
    )
    assert list(built.phospho.index) == list(phospho.index)
    assert list(built.site_metadata.columns) == ["gene_symbol", "site", "site_sequence"]
    assert built.transformation_state == TransformationState.raw(has_total_matrix=False)


def test_builder_rejects_invalid_site_metadata_structure() -> None:
    phospho = _load_rat_l6_phospho().head(20).copy(deep=True)
    bad_metadata = _site_metadata_for(phospho).iloc[1:, :]
    with pytest.raises(DatasetValidationError, match="must exactly match"):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=bad_metadata,
                organism=Organism.RAT,
            )
        )


def test_builder_rejects_missing_or_blank_site_sequence() -> None:
    phospho = _load_rat_l6_phospho().head(20).copy(deep=True)
    bad_metadata = _site_metadata_for(phospho)
    bad_metadata.loc[bad_metadata.index[0], "site_sequence"] = " "
    with pytest.raises(DatasetValidationError, match="site_sequence"):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=bad_metadata,
                organism=Organism.RAT,
            )
        )


def test_builder_rejects_unsupported_non_dataframe_input_mode() -> None:
    phospho = _load_rat_l6_phospho().head(10).copy(deep=True)
    with pytest.raises(UnsupportedInputFormatError, match="must be a pandas DataFrame"):
        DatasetBuildRequest(
            phospho="tests_legacy/fixtures/r_reference_l6/l6_phospho_matrix.csv",
            site_metadata=_site_metadata_for(phospho),
            organism=Organism.RAT,
        )


def test_builder_preserves_explicit_transformation_state() -> None:
    phospho = _load_rat_l6_phospho().head(20).copy(deep=True)
    total = pd.DataFrame(
        1.0,
        index=["PRKAA1", "MAPK1"],
        columns=phospho.columns.copy(),
    )
    explicit_state = TransformationState.raw(has_total_matrix=True)
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata_for(phospho),
            total=total,
            organism=Organism.RAT,
            transformation_state=explicit_state,
        )
    )
    assert built.transformation_state == explicit_state
    assert built.total is not None
    assert list(built.total.columns) == list(phospho.columns)


def test_reference_resolver_loads_real_bundled_rat_content() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )
    assert bundle.organism is Organism.RAT
    assert "AKT1" in set(bundle.kinase_substrate_map.loc[:, "kinase"])
    assert "MAPK14;Y182;" in set(bundle.site_sequences.index.astype(str))
    assert "KINASE_A" not in set(bundle.kinase_substrate_map.loc[:, "kinase"])


def test_reference_resolver_rejects_unsupported_bundled_human_preset() -> None:
    with pytest.raises(
        UnsupportedOrganismError, match="supported bundled organisms: rat"
    ):
        ReferenceResolver().run(
            ReferencePreset.HUMAN,
            dataset_organism=Organism.HUMAN,
        )


def test_reference_bundle_tables_are_structurally_coherent() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )
    substrate_sites = set(
        bundle.kinase_substrate_map.loc[:, "substrate_site"].astype(str)
    )
    known_sites = set(bundle.site_sequences.index.astype(str))
    assert substrate_sites.issubset(known_sites)


def test_simple_kinase_workflow_runs_end_to_end_with_nested_real_results() -> None:
    dataset = _build_dataset(n_sites=260)
    result = SimpleKinaseWorkflow().run(
        SimpleKinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=1),
            prediction_config=KinasePredictionConfig(top_k=6, ensemble_size=8),
            activity_config=KinaseActivityConfig(enabled=True, threshold=0.6),
        )
    )
    assert result.scoring_result.profile_scores.shape[0] == dataset.phospho.shape[0]
    assert result.scoring_result.profile_scores.shape[1] > 0
    assert result.scoring_result.combined_scores is not None
    assert result.prediction_result.pred_mat.shape[1] <= 8
    assert (result.prediction_result.pred_mat.to_numpy() >= 0.0).all()
    assert result.activity_result is not None
    assert set(result.activity_result.activity_scores.columns) == {
        "activity_score",
        "weighted_signal",
        "n_predicted_sites",
        "is_active",
    }
    assert not hasattr(result, "profile_scores")
    assert not hasattr(result, "combined_scores")
    assert not hasattr(result, "weights")
    assert not hasattr(result, "substrate_list")


def test_simple_kinase_workflow_activity_stage_is_optional() -> None:
    dataset = _build_dataset(n_sites=180)
    result = SimpleKinaseWorkflow().run(
        SimpleKinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=1),
            prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=6),
            activity_config=None,
        )
    )
    assert result.activity_result is None


def test_scoring_outputs_match_selected_legacy_profile_parity_values() -> None:
    dataset = _build_dataset(n_sites=None)
    result = SimpleKinaseWorkflow().run(
        SimpleKinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=1),
            prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=12),
            activity_config=None,
        )
    )
    expected = _load_expected_profile_scores()
    points = [
        ("AAK1;S677;", "AKT1"),
        ("ABCC4;S604;", "MAPK1"),
        ("ABI2;S165;", "PRKAA1"),
    ]
    for site_id, kinase in points:
        assert result.scoring_result.profile_scores.at[
            site_id, kinase
        ] == pytest.approx(
            expected.at[site_id, kinase],
            rel=1e-6,
            abs=1e-8,
        )


def test_prediction_top_sites_align_with_legacy_profile_ranking_subset() -> None:
    dataset = _build_dataset(n_sites=None)
    result = SimpleKinaseWorkflow().run(
        SimpleKinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=1),
            prediction_config=KinasePredictionConfig(top_k=3, ensemble_size=200),
            activity_config=None,
        )
    )
    expected = _load_expected_profile_scores()
    substrate_map = result.references.kinase_substrate_map
    for kinase in ("AKT1", "MAPK1"):
        candidates = [
            site_id
            for site_id in substrate_map.loc[
                substrate_map.loc[:, "kinase"] == kinase, "substrate_site"
            ].astype(str)
            if site_id in expected.index
            and site_id in result.scoring_result.profile_scores.index
        ]
        expected_top = expected.loc[candidates, kinase].astype(float).idxmax()
        observed_top = result.prediction_result.substrate_list.loc[
            (result.prediction_result.substrate_list.loc[:, "kinase"] == kinase)
            & (result.prediction_result.substrate_list.loc[:, "rank"] == 1),
            "substrate_site",
        ].iloc[0]
        assert observed_top == expected_top

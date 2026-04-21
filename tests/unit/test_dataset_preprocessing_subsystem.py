from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api.configs import (
    DatasetComparisonBuildingConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.api.requests import DatasetBuildRequest
from phospy.datasets.builders.contracts import (
    InterpretedDatasetBuildRequest,
    PreprocessedDatasetBuildTables,
)
from phospy.datasets.builders.executor import DatasetBuildExecutor
from phospy.datasets.builders.interpreter import DatasetBuildRequestInterpreter
from phospy.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.datasets.builders.transformation_resolver import ResolvedTransformation
from phospy.datasets.preprocessing.models import PreprocessingPlan, PreprocessingState
from phospy.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.errors.input import PhosPyInputError
from phospy.references.models import Organism
from phospy.transformations.models import TransformationState

ROOT = Path(__file__).resolve().parents[2]


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [2.0, 1.0, 4.0],
            "sample_c": [3.0, 1.5, 5.0],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
    )


def _site_metadata(index: pd.Index | None = None) -> pd.DataFrame:
    data = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
                "RPHFPQFSYSASGTA",
            ],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
    )
    return data if index is None else data.loc[index]


def _sample_metadata(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {"batch": [1, 1, 2]},
        index=columns,
    )


def _comparison_sample_metadata(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "comparison_group": ["group1", "group1", "group4", "group4"],
        },
        index=columns,
    )


def _total(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            str(columns[0]): [10.0, 12.0],
            str(columns[1]): [10.5, 12.5],
            str(columns[2]): [11.0, 13.0],
        },
        index=pd.Index(["MAPK14", "AKT1"], name="protein_id"),
    )


def test_preprocessing_pipeline_applies_plan_order() -> None:
    calls: list[str] = []

    class StageA:
        stage_key = "stage_a"

        def run(self, state: PreprocessingState) -> PreprocessingState:
            calls.append(self.stage_key)
            return state

    class StageB:
        stage_key = "stage_b"

        def run(self, state: PreprocessingState) -> PreprocessingState:
            calls.append(self.stage_key)
            return state

    state = PreprocessingState(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="forbid",
            missing_data_min_observed_values=None,
            total_protein_correction_policy="none",
            site_matrix_policy="as_input",
            comparison_building_policy="none",
            stage_order=("stage_b", "stage_a"),
        ),
    )

    pipeline = PreprocessingPipeline(stage_registry=(StageA(), StageB()))
    observed = pipeline.run(state)

    assert observed is state
    assert calls == ["stage_b", "stage_a"]


def test_preprocessing_pipeline_passes_stage_state_forward() -> None:
    observed_first_value: list[float] = []

    class AddOneStage:
        stage_key = "add_one"

        def run(self, state: PreprocessingState) -> PreprocessingState:
            return replace(state, phospho=state.phospho + 1.0)

    class InspectStage:
        stage_key = "inspect"

        def run(self, state: PreprocessingState) -> PreprocessingState:
            observed_first_value.append(float(state.phospho.iloc[0, 0]))
            return replace(state, site_metadata=state.site_metadata.assign(seen=True))

    phospho = _phospho()
    sample_metadata = _sample_metadata(phospho.columns)
    total = _total(phospho.columns)
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(),
        sample_metadata=sample_metadata,
        total=total,
        plan=PreprocessingPlan(
            missing_data_policy="forbid",
            missing_data_min_observed_values=None,
            total_protein_correction_policy="none",
            site_matrix_policy="as_input",
            comparison_building_policy="none",
            stage_order=("add_one", "inspect"),
        ),
    )
    pipeline = PreprocessingPipeline(stage_registry=(AddOneStage(), InspectStage()))

    observed = pipeline.run(state)

    assert observed_first_value == [2.0]
    assert "seen" in observed.site_metadata.columns
    assert observed.sample_metadata is sample_metadata
    assert observed.total is total


def test_dataset_preprocessor_regression_forbid_policy_is_passthrough() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    sample_metadata = _sample_metadata(phospho.columns)
    total = _total(phospho.columns)

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=total,
        plan=PreprocessingPlan.from_config(DatasetPreprocessingConfig()),
    )

    assert preprocessed.phospho is phospho
    assert preprocessed.site_metadata is site_metadata
    assert preprocessed.sample_metadata is sample_metadata
    assert preprocessed.total is total


def test_dataset_preprocessor_regression_impute_row_median_policy() -> None:
    phospho = _phospho()
    phospho.loc["MAPK14;Y182;", "sample_a"] = float("nan")
    phospho.loc["GSK3B;S9;", :] = float("nan")
    site_metadata = _site_metadata()
    sample_metadata = _sample_metadata(phospho.columns)
    total = _total(phospho.columns)

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=total,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=2,
                )
            )
        ),
    )

    expected_phospho = pd.DataFrame(
        {
            "sample_a": [2.5, 3.0],
            "sample_b": [2.0, 4.0],
            "sample_c": [3.0, 5.0],
        },
        index=["MAPK14;Y182;", "AKT1;T308;"],
    )
    expected_site_metadata = _site_metadata(pd.Index(["MAPK14;Y182;", "AKT1;T308;"]))

    pdt.assert_frame_equal(preprocessed.phospho, expected_phospho)
    pdt.assert_frame_equal(preprocessed.site_metadata, expected_site_metadata)
    assert preprocessed.sample_metadata is sample_metadata
    assert preprocessed.total is total


def test_preprocessing_plan_orders_site_matrix_after_total_correction() -> None:
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy="ratio_to_total"
            ),
            site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata"),
        )
    )
    assert plan.stage_order == (
        "missing_data",
        "total_protein_correction",
        "site_matrix",
    )


def test_preprocessing_plan_orders_comparisons_after_upstream_stages() -> None:
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy="ratio_to_total"
            ),
            site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata"),
            comparisons=DatasetComparisonBuildingConfig(policy="sample_metadata_pairs"),
        )
    )
    assert plan.stage_order == (
        "missing_data",
        "total_protein_correction",
        "site_matrix",
        "comparisons",
    )


def test_dataset_preprocessor_applies_total_protein_correction_ratio_policy() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    total = pd.DataFrame(
        {
            "sample_a": [0.5, 2.0, 1.0],
            "sample_b": [1.0, 1.0, 1.0],
            "sample_c": [1.5, 0.5, 1.0],
        },
        index=pd.Index(["MAPK14", "GSK3B", "AKT1"], name="protein_id"),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=total,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="ratio_to_total"
                )
            )
        ),
    )

    expected = pd.DataFrame(
        {
            "sample_a": [0.5, 0.0, 2.0],
            "sample_b": [1.0, 0.0, 3.0],
            "sample_c": [1.5, 1.0, 4.0],
        },
        index=phospho.index.copy(),
    )
    pdt.assert_frame_equal(preprocessed.phospho, expected)
    assert preprocessed.site_metadata is site_metadata
    assert preprocessed.total is total


def test_dataset_preprocessor_builds_site_matrix_from_metadata_policy() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [2.0, 8.0, 5.0, 1.0],
            "sample_b": [2.5, 8.5, 5.5, float("nan")],
        },
        index=pd.Index(
            ["row_a", "row_b", "row_c", "row_d"],
            name="input_row",
        ),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "Y182", "T308", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_B", "", "SEQ_D"],
            "source_uid": ["UID_A", "UID_B", "UID_C", "UID_D"],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            )
        ),
    )

    expected_phospho = pd.DataFrame(
        {
            "sample_a": [8.0],
            "sample_b": [8.5],
        },
        index=pd.Index(["MAPK14;Y182;"], name="input_row"),
    )
    expected_site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_B"],
            "source_uid": ["UID_B"],
        },
        index=expected_phospho.index.copy(),
    )

    pdt.assert_frame_equal(preprocessed.phospho, expected_phospho)
    pdt.assert_frame_equal(preprocessed.site_metadata, expected_site_metadata)


def test_dataset_preprocessor_site_matrix_retain_missing_policy_keeps_partial_rows() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [2.0, 8.0],
            "sample_b": [float("nan"), 8.5],
        },
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    missing_data_policy="retain_missing",
                )
            )
        ),
    )

    assert preprocessed.phospho.index.tolist() == ["AKT1;T308;", "MAPK14;Y182;"]
    assert pd.isna(preprocessed.phospho.loc["MAPK14;Y182;", "sample_b"])
    stats = preprocessed.phospho.attrs.get("site_matrix_row_drop_stats")
    assert stats is not None
    assert stats["missing_data_policy"] == "retain_missing"
    assert stats["required_observed_count"] == 0


def test_dataset_preprocessor_site_matrix_supports_min_observed_and_duplicate_aggregate_mean() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0, float("nan")],
            "sample_b": [2.0, 4.0, float("nan")],
            "sample_c": [float("nan"), float("nan"), 9.0],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1"],
            "site": ["Y182", "Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C"],
            "uid": ["A", "B", "C"],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    missing_data_policy="require_min_observed_values",
                    minimum_observed_values=1,
                    duplicate_site_strategy="aggregate_mean",
                )
            )
        ),
    )

    assert preprocessed.phospho.index.tolist() == ["AKT1;T308;", "MAPK14;Y182;"]
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_a"] == pytest.approx(2.0)
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_b"] == pytest.approx(3.0)
    assert pd.isna(preprocessed.phospho.loc["MAPK14;Y182;", "sample_c"])
    assert preprocessed.site_metadata.loc["MAPK14;Y182;", "uid"] == "A"


def test_dataset_preprocessor_site_matrix_duplicate_first_strategy_keeps_first_row() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
            "uid": ["A", "B"],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    duplicate_site_strategy="first",
                    missing_data_policy="retain_missing",
                )
            )
        ),
    )

    assert preprocessed.phospho.index.tolist() == ["MAPK14;Y182;"]
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_a"] == pytest.approx(1.0)
    assert preprocessed.site_metadata.loc["MAPK14;Y182;", "uid"] == "A"


def test_dataset_preprocessor_site_matrix_duplicate_aggregate_median_strategy() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 30.0, 3.0],
            "sample_b": [2.0, 40.0, 4.0],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1"],
            "site": ["Y182", "Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C"],
            "uid": ["A", "B", "C"],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    duplicate_site_strategy="aggregate_median",
                    missing_data_policy="retain_missing",
                )
            )
        ),
    )

    assert preprocessed.phospho.index.tolist() == ["AKT1;T308;", "MAPK14;Y182;"]
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_a"] == pytest.approx(15.5)
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_b"] == pytest.approx(21.0)
    assert preprocessed.site_metadata.loc["MAPK14;Y182;", "uid"] == "A"


def test_dataset_preprocessor_rejects_site_matrix_min_observed_above_sample_count() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
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

    with pytest.raises(
        PhosPyInputError,
        match="minimum_observed_values cannot exceed phospho sample count",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(
                        policy="build_from_metadata",
                        missing_data_policy="require_min_observed_values",
                        minimum_observed_values=3,
                    )
                )
            ),
        )


def test_dataset_preprocessor_rejects_site_matrix_duplicate_rows_in_error_mode() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 3.0], "sample_b": [2.0, 4.0]},
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="duplicate_site_strategy='error'",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(
                        policy="build_from_metadata",
                        duplicate_site_strategy="error",
                    )
                )
            ),
        )


def test_dataset_preprocessor_rejects_site_matrix_build_without_site_sequence() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata().drop(columns=["site_sequence"])

    with pytest.raises(
        PhosPyInputError,
        match="site-matrix construction requires site_metadata columns: site_sequence",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
                )
            ),
        )


def test_dataset_preprocessor_rejects_correction_when_total_columns_mismatch() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    total = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [1.0, 2.0, 3.0],
        },
        index=pd.Index(["MAPK14", "GSK3B", "AKT1"], name="protein_id"),
    )

    with pytest.raises(
        PhosPyInputError,
        match="requires total columns to exactly match phospho columns",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=total,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    total_protein_correction=DatasetTotalProteinCorrectionConfig(
                        policy="ratio_to_total"
                    )
                )
            ),
        )


def test_dataset_preprocessor_rejects_correction_when_proteins_are_unmatched() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    total = _total(phospho.columns)

    with pytest.raises(
        PhosPyInputError,
        match="requires complete phospho/total matching but would drop",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=total,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    total_protein_correction=DatasetTotalProteinCorrectionConfig(
                        policy="ratio_to_total"
                    )
                )
            ),
        )


def test_dataset_preprocessor_total_protein_correction_matches_legacy_donor_fixture() -> (
    None
):
    phospho_fixture = pd.read_csv(
        ROOT / "tests_legacy" / "fixtures" / "r_reference" / "df_phospho_filtered.csv"
    )
    total_fixture = pd.read_csv(
        ROOT / "tests_legacy" / "fixtures" / "r_reference" / "df_total_filtered.csv"
    )
    corrected_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "protein_correction"
        / "legacy_r_reference_corrected_matrix.csv"
    )

    phospho_columns = tuple(f"p_group{group}" for group in range(1, 7))
    total_columns = tuple(f"group{group}" for group in range(1, 7))

    site_tokens = (
        phospho_fixture.loc[:, "gene_p_site"]
        .astype(str)
        .str.split("_", n=1, expand=True)
    )
    site_index = pd.Index(
        site_tokens.loc[:, 0].astype(str)
        + ";"
        + site_tokens.loc[:, 1].astype(str)
        + ";",
        name="site_id",
    )

    phospho = phospho_fixture.loc[:, list(phospho_columns)].copy()
    phospho.columns = pd.Index(phospho_columns)
    phospho.index = site_index
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": phospho_fixture.loc[:, "gene_names"].astype(str).tolist(),
            "site": site_tokens.loc[:, 1].astype(str).tolist(),
        },
        index=site_index.copy(),
    )

    total = total_fixture.loc[:, list(total_columns)].copy()
    total.columns = pd.Index(phospho_columns)
    total.index = pd.Index(total_fixture.loc[:, "genes"].astype(str), name="protein_id")

    expected = (
        corrected_fixture.set_index("site_id")
        .loc[:, list(phospho_columns)]
        .astype(float)
        .copy()
    )
    expected.index = site_index.copy()

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=total,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="ratio_to_total"
                )
            )
        ),
    )

    pdt.assert_frame_equal(preprocessed.phospho, expected)


def test_dataset_preprocessor_site_matrix_build_matches_legacy_donor_fixture() -> None:
    corrected_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "site_matrix"
        / "legacy_r_reference_phospho_corrected.csv"
    )
    expected_matrix_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "site_matrix"
        / "legacy_r_reference_expected_matrix.csv"
    )
    expected_input_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "site_matrix"
        / "legacy_r_reference_expected_phosr_input.csv"
    )

    corrected_cols = tuple(f"phospho_corrected_{position}" for position in range(1, 7))
    phospho = corrected_fixture.loc[:, list(corrected_cols)].astype(float).copy()
    phospho.index = pd.Index(
        corrected_fixture.loc[:, "uid"].astype(str),
        name="source_uid",
    )

    site_tokens = (
        corrected_fixture.loc[:, "gene_p_site"]
        .astype(str)
        .str.split("_", n=1, expand=True)
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": corrected_fixture.loc[:, "gene_names"].astype(str).tolist(),
            "site": site_tokens.loc[:, 1].astype(str).tolist(),
            "site_sequence": corrected_fixture.loc[:, "centralized_sequence"]
            .astype(str)
            .tolist(),
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            )
        ),
    )

    expected_phospho = (
        expected_matrix_fixture.set_index(expected_matrix_fixture.columns[0])
        .loc[:, list(corrected_cols)]
        .astype(float)
    )
    expected_phospho.index = pd.Index(
        expected_phospho.index.astype(str), name="source_uid"
    )
    expected_site_metadata = pd.DataFrame(
        {
            "gene_symbol": expected_input_fixture.loc[:, "gene_names"]
            .astype(str)
            .tolist(),
            "site": expected_input_fixture.loc[:, "p_site"].astype(str).tolist(),
            "site_sequence": expected_input_fixture.loc[:, "centralized_sequence"]
            .astype(str)
            .tolist(),
        },
        index=pd.Index(
            expected_input_fixture.loc[:, "site_id"].astype(str),
            name="source_uid",
        ),
    )

    pdt.assert_frame_equal(preprocessed.phospho, expected_phospho)
    pdt.assert_frame_equal(
        preprocessed.site_metadata.loc[:, ["gene_symbol", "site", "site_sequence"]],
        expected_site_metadata,
    )


def test_dataset_preprocessor_builds_inferred_comparisons_from_sample_groups() -> None:
    phospho = pd.DataFrame(
        {
            "sample_1": [8.0, 2.0],
            "sample_2": [8.0, 4.0],
            "sample_3": [5.0, 1.0],
            "sample_4": [5.0, 1.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = _site_metadata(phospho.index)
    sample_metadata = _comparison_sample_metadata(phospho.columns)

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                comparisons=DatasetComparisonBuildingConfig(
                    policy="sample_metadata_pairs"
                )
            )
        ),
    )

    assert preprocessed.comparisons is not None
    expected = pd.DataFrame(
        {"p_group1_group4": [3.0, 2.0]},
        index=phospho.index.copy(),
    )
    pdt.assert_frame_equal(preprocessed.comparisons, expected)


def test_dataset_preprocessor_builds_explicit_comparison_pairs() -> None:
    phospho = pd.DataFrame(
        {
            "sample_1": [8.0],
            "sample_2": [8.0],
            "sample_3": [5.0],
            "sample_4": [5.0],
        },
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = _site_metadata(phospho.index)
    sample_metadata = _comparison_sample_metadata(phospho.columns)

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                comparisons=DatasetComparisonBuildingConfig(
                    policy="sample_metadata_pairs",
                    pairs=(("group4", "group1"),),
                )
            )
        ),
    )

    assert preprocessed.comparisons is not None
    expected = pd.DataFrame(
        {"p_group4_group1": [-3.0]},
        index=phospho.index.copy(),
    )
    pdt.assert_frame_equal(preprocessed.comparisons, expected)


def test_dataset_preprocessor_rejects_comparison_building_without_sample_metadata() -> (
    None
):
    with pytest.raises(
        PhosPyInputError,
        match="policy='sample_metadata_pairs' requires sample_metadata input data",
    ):
        DatasetPreprocessor().run(
            phospho=_phospho().iloc[:, :2].copy(deep=True),
            site_metadata=_site_metadata().iloc[:3].copy(deep=True),
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    comparisons=DatasetComparisonBuildingConfig(
                        policy="sample_metadata_pairs"
                    )
                )
            ),
        )


def test_dataset_preprocessor_rejects_comparison_building_without_group_column() -> (
    None
):
    phospho = _phospho().iloc[:, :2].copy(deep=True)
    site_metadata = _site_metadata()
    sample_metadata = pd.DataFrame({"batch": [1, 2]}, index=phospho.columns)

    with pytest.raises(
        PhosPyInputError,
        match="requires sample_metadata column 'comparison_group'",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    comparisons=DatasetComparisonBuildingConfig(
                        policy="sample_metadata_pairs"
                    )
                )
            ),
        )


def test_dataset_preprocessor_comparison_building_matches_legacy_pairwise_expectation() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [7.0],
            "sample_b": [4.0],
        },
        index=pd.Index(["PRKACA;S339;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["PRKACA"],
            "site": ["S339"],
            "site_sequence": ["AAAAAA"],
        },
        index=phospho.index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {"comparison_group": ["sample_a", "sample_b"]},
        index=phospho.columns.copy(),
    )
    expected_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "comparison_building"
        / "legacy_pairwise_expected.csv"
    ).set_index("site_id")

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                comparisons=DatasetComparisonBuildingConfig(
                    policy="sample_metadata_pairs",
                    pairs=(("sample_a", "sample_b"),),
                )
            )
        ),
    )

    assert preprocessed.comparisons is not None
    expected = expected_fixture.astype(float).copy()
    expected.index = phospho.index.copy()
    pdt.assert_frame_equal(preprocessed.comparisons, expected)


def test_executor_delegates_preprocessing_to_internal_subsystem() -> None:
    calls: list[str] = []
    phospho = _phospho().iloc[:2].copy(deep=True)
    site_metadata = _site_metadata(phospho.index).copy(deep=True)
    sample_metadata = _sample_metadata(phospho.columns)
    total = _total(phospho.columns)

    interpreted = InterpretedDatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=total,
        organism=Organism.RAT,
        preprocessing_plan=PreprocessingPlan.default(),
    )

    preprocessed_tables = PreprocessedDatasetBuildTables(
        phospho=phospho + 10.0,
        site_metadata=site_metadata.assign(processed=True),
        sample_metadata=sample_metadata.assign(processed=True),
        total=total + 5.0,
    )

    class PreprocessorSpy:
        def run(
            self,
            *,
            phospho: pd.DataFrame,
            site_metadata: pd.DataFrame,
            sample_metadata: pd.DataFrame | None,
            total: pd.DataFrame | None,
            plan: PreprocessingPlan,
        ) -> PreprocessedDatasetBuildTables:
            calls.append("preprocessor")
            assert phospho is interpreted.phospho
            assert site_metadata is interpreted.site_metadata
            assert sample_metadata is interpreted.sample_metadata
            assert total is interpreted.total
            assert plan is interpreted.preprocessing_plan
            return preprocessed_tables

    class ResolverSpy:
        def run(
            self,
            *,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None,
        ) -> ResolvedTransformation:
            calls.append("resolver")
            assert phospho is preprocessed_tables.phospho
            assert total is preprocessed_tables.total
            return ResolvedTransformation(
                phospho=phospho,
                total=total,
                transformation_state=TransformationState.established_raw(
                    has_total_matrix=True
                ),
            )

    built = DatasetBuildExecutor(
        preprocessor=PreprocessorSpy(),
        transformation_resolver=ResolverSpy(),
    ).run(interpreted)

    pdt.assert_frame_equal(built.phospho, preprocessed_tables.phospho)
    pdt.assert_frame_equal(built.site_metadata, preprocessed_tables.site_metadata)
    assert built.sample_metadata is preprocessed_tables.sample_metadata
    pdt.assert_frame_equal(built.total, preprocessed_tables.total)
    assert calls == ["preprocessor", "resolver"]


def test_dataset_interpreter_does_not_apply_preprocessing_science() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan")],
            "sample_b": [2.0, float("nan")],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
            ],
        },
        index=phospho.index,
    )
    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        preprocessing_config=DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(
                policy="impute_row_median",
                min_observed_values=2,
            ),
        ),
    )

    interpreted = DatasetBuildRequestInterpreter().run(request)

    assert interpreted.phospho.isna().to_numpy().sum() == 2
    assert interpreted.phospho.index.tolist() == ["MAPK14;Y182;", "GSK3B;S9;"]

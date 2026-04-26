from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api.configs import (
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
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
from phospy.datasets.builders.transformation_resolver import ResolvedIntensityScale
from phospy.datasets.preprocessing.models import PreprocessingPlan, PreprocessingState
from phospy.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.errors.input import PhosPyInputError
from phospy.references.models import Organism
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
)

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


def _internal_site_matrix_config(
    *,
    policy: str = "build_from_metadata",
    duplicate_site_policy: str = "max_mean_signal",
    missing_data_policy: str = "drop_any_missing",
    minimum_observed_values: int | None = None,
) -> DatasetSiteMatrixConfig:
    """Construct a site-matrix config bypassing public constructor validation.

    Subsystem tests use this to exercise internal preprocessing-lane behaviors
    that are intentionally unsupported on the public config surface.
    """

    config = object.__new__(DatasetSiteMatrixConfig)
    object.__setattr__(config, "policy", policy)
    object.__setattr__(config, "duplicate_site_policy", duplicate_site_policy)
    object.__setattr__(config, "missing_data_policy", missing_data_policy)
    object.__setattr__(config, "minimum_observed_values", minimum_observed_values)
    return config


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
    assert preprocessed.row_audit is not None
    assert preprocessed.row_audit.empty


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
    dropped = preprocessed.row_audit.loc[
        (preprocessed.row_audit.loc[:, "stage"] == "missing_data")
        & (preprocessed.row_audit.loc[:, "action"] == "dropped")
    ]
    assert set(dropped.loc[:, "source_row_id"].astype(str)) == {"GSK3B;S9;"}
    imputed = preprocessed.row_audit.loc[
        (preprocessed.row_audit.loc[:, "stage"] == "missing_data")
        & (preprocessed.row_audit.loc[:, "action"] == "imputed")
    ]
    assert set(imputed.loc[:, "source_row_id"].astype(str)) == {"MAPK14;Y182;"}
    assert "imputed_columns" in imputed.iloc[0]["parameter_snapshot"]


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
                site_matrix=_internal_site_matrix_config(
                    policy="build_from_metadata",
                    missing_data_policy="retain_missing",
                )
            )
        ),
    )

    assert preprocessed.phospho.index.tolist() == ["AKT1;T308;", "MAPK14;Y182;"]
    assert pd.isna(preprocessed.phospho.loc["MAPK14;Y182;", "sample_b"])
    dropped = preprocessed.row_audit.loc[
        (preprocessed.row_audit.loc[:, "stage"] == "site_matrix")
        & (preprocessed.row_audit.loc[:, "action"] == "dropped")
    ]
    assert dropped.empty


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
                site_matrix=_internal_site_matrix_config(
                    policy="build_from_metadata",
                    missing_data_policy="require_min_observed_values",
                    minimum_observed_values=1,
                    duplicate_site_policy="aggregate_mean",
                )
            )
        ),
    )

    assert preprocessed.phospho.index.tolist() == ["AKT1;T308;", "MAPK14;Y182;"]
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_a"] == pytest.approx(2.0)
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_b"] == pytest.approx(3.0)
    assert pd.isna(preprocessed.phospho.loc["MAPK14;Y182;", "sample_c"])
    assert preprocessed.site_metadata.loc["MAPK14;Y182;", "uid"] == "A"


def test_dataset_preprocessor_site_matrix_duplicate_first_policy_keeps_first_row() -> (
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
                site_matrix=_internal_site_matrix_config(
                    policy="build_from_metadata",
                    duplicate_site_policy="first",
                    missing_data_policy="retain_missing",
                )
            )
        ),
    )

    assert preprocessed.phospho.index.tolist() == ["MAPK14;Y182;"]
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_a"] == pytest.approx(1.0)
    assert preprocessed.site_metadata.loc["MAPK14;Y182;", "uid"] == "A"


def test_dataset_preprocessor_site_matrix_duplicate_aggregate_median_policy() -> None:
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
                site_matrix=_internal_site_matrix_config(
                    policy="build_from_metadata",
                    duplicate_site_policy="aggregate_median",
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
                    site_matrix=_internal_site_matrix_config(
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
        match="duplicate_site_policy='error'",
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
                        duplicate_site_policy="error",
                    )
                )
            ),
        )


def test_dataset_preprocessor_rejects_site_matrix_build_without_site_sequence() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata().drop(columns=["site_sequence"])

    with pytest.raises(
        PhosPyInputError,
        match=(
            "site-matrix construction produced no retained rows after filtering; "
            "input_rows=3, dropped_missing_sequence=3"
        ),
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


def test_dataset_preprocessor_total_protein_correction_matches_historical_baseline_fixture() -> (
    None
):
    phospho_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "protein_correction"
        / "reference_input_phospho.csv"
    ).set_index("site_id")
    site_metadata_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "protein_correction"
        / "reference_input_site_metadata.csv"
    ).set_index("site_id")
    total_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "protein_correction"
        / "reference_input_total.csv"
    ).set_index("protein_id")
    corrected_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "protein_correction"
        / "reference_corrected_matrix.csv"
    ).set_index("site_id")

    phospho_columns = tuple(phospho_fixture.columns.astype(str))
    phospho = phospho_fixture.astype(float).copy()
    phospho.columns = pd.Index(phospho_columns)
    phospho.index = pd.Index(phospho_fixture.index.astype(str), name="site_id")
    site_metadata = site_metadata_fixture.copy()
    site_metadata.index = pd.Index(
        site_metadata_fixture.index.astype(str),
        name="site_id",
    )
    total = total_fixture.astype(float).copy()
    total.columns = pd.Index(phospho_columns)
    total.index = pd.Index(total_fixture.index.astype(str), name="protein_id")

    expected = corrected_fixture.loc[:, list(phospho_columns)].astype(float).copy()
    expected.index = pd.Index(corrected_fixture.index.astype(str), name="site_id")

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


def test_dataset_preprocessor_site_matrix_build_matches_historical_baseline_fixture() -> (
    None
):
    corrected_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "site_matrix"
        / "reference_phospho_corrected.csv"
    )
    expected_matrix_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "site_matrix"
        / "reference_expected_matrix.csv"
    )
    expected_input_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "site_matrix"
        / "reference_expected_input.csv"
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


def test_dataset_preprocessor_comparison_building_matches_reference_pairwise_expectation() -> (
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
        / "reference_pairwise_expected.csv"
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
    assert preprocessed.comparison_group_stats is not None
    assert preprocessed.comparison_pair_stats is not None
    assert {"n", "mean", "sd", "sem"}.issubset(
        set(preprocessed.comparison_group_stats.columns)
    )
    assert {
        "left_n",
        "right_n",
        "left_mean",
        "right_mean",
        "left_sd",
        "right_sd",
        "left_sem",
        "right_sem",
        "effect_size",
    }.issubset(set(preprocessed.comparison_pair_stats.columns))
    pair_stats = preprocessed.comparison_pair_stats.copy(deep=True)
    assert pair_stats.shape[0] == preprocessed.comparisons.shape[0]
    paired_with_matrix = pair_stats.merge(
        preprocessed.comparisons.reset_index()
        .rename(columns={"index": "site_id"})
        .rename(columns={"p_sample_a_sample_b": "expected_effect_size"}),
        how="inner",
        on="site_id",
    )
    assert paired_with_matrix.shape[0] == preprocessed.comparisons.shape[0]
    assert (
        paired_with_matrix.loc[:, "effect_size"]
        == paired_with_matrix.loc[:, "expected_effect_size"]
    ).all()


def test_dataset_preprocessor_comparison_stats_follow_pandas_single_sample_convention() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_1": [8.0, 2.0],
            "sample_2": [10.0, 4.0],
            "sample_3": [5.0, 1.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = _site_metadata(phospho.index)
    sample_metadata = pd.DataFrame(
        {"comparison_group": ["group1", "group1", "group2"]},
        index=phospho.columns.copy(),
    )

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
        {"p_group1_group2": [4.0, 2.0]},
        index=phospho.index.copy(),
    )
    pdt.assert_frame_equal(preprocessed.comparisons, expected)
    assert preprocessed.comparison_group_stats is not None
    assert preprocessed.comparison_pair_stats is not None

    group_stats = preprocessed.comparison_group_stats
    group2_rows = group_stats.loc[group_stats["group"] == "group2"].sort_values(
        "site_id"
    )
    assert group2_rows["n"].tolist() == [1, 1]
    assert group2_rows["mean"].tolist() == [1.0, 5.0]
    assert group2_rows["sd"].isna().all()
    assert group2_rows["sem"].isna().all()

    pair_stats = preprocessed.comparison_pair_stats.sort_values("site_id")
    assert pair_stats["comparison"].tolist() == ["p_group1_group2", "p_group1_group2"]
    assert pair_stats["left_n"].tolist() == [2, 2]
    assert pair_stats["right_n"].tolist() == [1, 1]
    assert pair_stats["left_mean"].tolist() == [3.0, 9.0]
    assert pair_stats["right_mean"].tolist() == [1.0, 5.0]
    assert pair_stats["right_sd"].isna().all()
    assert pair_stats["right_sem"].isna().all()
    assert pair_stats["effect_size"].tolist() == [2.0, 4.0]


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
            expected_scale_kind: object | None = None,
        ) -> ResolvedIntensityScale:
            calls.append("resolver")
            assert phospho is preprocessed_tables.phospho
            assert total is preprocessed_tables.total
            assert expected_scale_kind is not None
            return ResolvedIntensityScale(
                phospho=phospho,
                total=total,
                intensity_scale_state=supported_linear_intensity_scale_state(
                    has_total_matrix=True
                ),
            )

    built = DatasetBuildExecutor(
        preprocessor=PreprocessorSpy(),
        intensity_scale_resolver=ResolverSpy(),
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


def test_preprocessing_plan_defaults_keep_identity_transform_and_no_normalisation() -> (
    None
):
    plan = PreprocessingPlan.from_config(DatasetPreprocessingConfig())
    assert plan.intensity_transform_policy == "identity"
    assert plan.intensity_transform_pseudocount == pytest.approx(1.0)
    assert plan.normalisation_policy == "none"
    assert "intensity_transform" not in plan.stage_order
    assert "normalisation" not in plan.stage_order


def test_dataset_preprocessor_applies_log2_intensity_transform_policy() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            )
        ),
    )

    expected = pd.DataFrame(
        {
            "sample_a": [1.0, 1.584962500721156, 2.0],
            "sample_b": [1.584962500721156, 1.0, 2.321928094887362],
            "sample_c": [2.0, 1.3219280948873624, 2.584962500721156],
        },
        index=phospho.index.copy(),
    )
    pdt.assert_frame_equal(preprocessed.phospho, expected)
    assert preprocessed.phospho.index.equals(phospho.index)
    assert preprocessed.phospho.columns.equals(phospho.columns)
    assert preprocessed.site_metadata is site_metadata


def test_dataset_preprocessor_rejects_log2_transform_when_values_are_invalid() -> None:
    phospho = _phospho()
    phospho.loc["MAPK14;Y182;", "sample_a"] = -1.0

    with pytest.raises(
        PhosPyInputError,
        match="requires all non-missing values plus pseudocount to be greater than 0",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=_site_metadata(),
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    intensity_transform=DatasetIntensityTransformConfig(
                        policy="log2",
                        pseudocount=1.0,
                    )
                )
            ),
        )


def test_dataset_preprocessor_applies_median_center_normalisation_policy() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [4.0, 2.0, 0.0],
        },
        index=_phospho().index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                normalisation=DatasetNormalisationConfig(policy="median_center")
            )
        ),
    )

    expected = pd.DataFrame(
        {
            "sample_a": [-1.0, 0.0, 1.0],
            "sample_b": [2.0, 0.0, -2.0],
        },
        index=phospho.index.copy(),
    )
    pdt.assert_frame_equal(preprocessed.phospho, expected)
    assert preprocessed.phospho.index.equals(phospho.index)
    assert preprocessed.phospho.columns.equals(phospho.columns)


def test_dataset_preprocessor_quantile_normalisation_equalises_column_distributions() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [5.0, 2.0, 3.0, 4.0],
            "sample_b": [4.0, 1.0, 2.5, 2.0],
            "sample_c": [8.0, 6.0, 7.0, 9.0],
        },
        index=pd.Index(["A;S1;", "B;S1;", "C;S1;", "D;S1;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["A", "B", "C", "D"],
            "site": ["S1", "S1", "S1", "S1"],
            "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C", "SEQ_D"],
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
                normalisation=DatasetNormalisationConfig(policy="quantile")
            )
        ),
    )

    sorted_a = sorted(preprocessed.phospho.loc[:, "sample_a"].tolist())
    sorted_b = sorted(preprocessed.phospho.loc[:, "sample_b"].tolist())
    sorted_c = sorted(preprocessed.phospho.loc[:, "sample_c"].tolist())
    assert sorted_a == pytest.approx(sorted_b)
    assert sorted_b == pytest.approx(sorted_c)
    assert preprocessed.phospho.index.equals(phospho.index)
    assert preprocessed.phospho.columns.equals(phospho.columns)
    assert preprocessed.site_metadata is site_metadata


def test_dataset_preprocessor_quantile_normalisation_is_deterministic_with_ties_and_missing() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 1.0, 3.0, float("nan")],
            "sample_b": [4.0, 2.0, 2.0, 1.0],
            "sample_c": [7.0, 8.0, 9.0, 10.0],
        },
        index=pd.Index(["A;S1;", "B;S1;", "C;S1;", "D;S1;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["A", "B", "C", "D"],
            "site": ["S1", "S1", "S1", "S1"],
            "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C", "SEQ_D"],
        },
        index=phospho.index.copy(),
    )
    config = DatasetPreprocessingConfig(
        normalisation=DatasetNormalisationConfig(policy="quantile")
    )

    first = DatasetPreprocessor().run(
        phospho=phospho.copy(deep=True),
        site_metadata=site_metadata.copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(config),
    )
    second = DatasetPreprocessor().run(
        phospho=phospho.copy(deep=True),
        site_metadata=site_metadata.copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(config),
    )

    pdt.assert_frame_equal(first.phospho, second.phospho)
    assert first.phospho.loc["A;S1;", "sample_a"] == pytest.approx(
        first.phospho.loc["B;S1;", "sample_a"]
    )
    assert pd.isna(first.phospho.loc["D;S1;", "sample_a"])


@pytest.mark.parametrize(
    ("intensity_transform", "normalisation", "expected"),
    [
        (
            DatasetIntensityTransformConfig(policy="log2", pseudocount=1.0),
            DatasetNormalisationConfig(policy="none"),
            "intensity_transform.policy='log2' requires numeric phospho columns",
        ),
        (
            DatasetIntensityTransformConfig(policy="identity", pseudocount=1.0),
            DatasetNormalisationConfig(policy="median_center"),
            "normalisation.policy='median_center' requires numeric phospho columns",
        ),
    ],
)
def test_dataset_preprocessor_rejects_non_numeric_phospho_columns_for_new_methods(
    intensity_transform: DatasetIntensityTransformConfig,
    normalisation: DatasetNormalisationConfig,
    expected: str,
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "metadata_note": ["a", "b"],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(PhosPyInputError, match=expected):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    intensity_transform=intensity_transform,
                    normalisation=normalisation,
                )
            ),
        )

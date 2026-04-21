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
    row_drop_stats = built.phospho.attrs.get("site_matrix_row_drop_stats")
    assert row_drop_stats is not None
    assert row_drop_stats["missing_data_policy"] == missing_data_policy
    assert row_drop_stats["retained_rows"] == 2


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
    row_drop_stats = built.phospho.attrs.get("site_matrix_row_drop_stats")
    assert row_drop_stats is not None
    assert row_drop_stats["input_rows"] == 4
    assert row_drop_stats["dropped_missing_sequence"] == 2
    assert row_drop_stats["retained_rows"] == 2

from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api.configs import (
    DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER,
    DATASET_LOCALISATION_MODE_REQUIRE_THRESHOLD,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
)
from phospy.api.requests import DatasetBuildRequest
from phospy.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.datasets.preprocessing.models import PreprocessingPlan
from phospy.errors.input import PhosPyInputError


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0, 2.0, 3.0]},
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"], name="site_id"),
    )


def _site_metadata(
    *,
    include_localisation_column: bool = True,
    localisation_values: list[object] | None = None,
) -> pd.DataFrame:
    phospho = _phospho()
    payload: dict[str, list[object]] = {
        "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
        "site": ["Y182", "T308", "S9"],
        "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C"],
    }
    if include_localisation_column:
        payload["localisation_confidence"] = (
            [0.95, 0.9, 0.92] if localisation_values is None else localisation_values
        )
    return pd.DataFrame(payload, index=phospho.index.copy())


def _run_with_localisation(
    *,
    mode: str,
    min_confidence: float = 0.75,
    waiver_reason: str | None = None,
    site_metadata: pd.DataFrame | None = None,
):
    phospho = _phospho()
    resolved_site_metadata = (
        _site_metadata() if site_metadata is None else site_metadata.copy(deep=True)
    )
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            localisation=DatasetLocalisationConfig(
                mode=mode,  # type: ignore[arg-type]
                min_confidence=min_confidence,
                waiver_reason=waiver_reason,
            )
        )
    )
    return DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=resolved_site_metadata,
        sample_metadata=None,
        total=None,
        plan=plan,
    )


def test_localisation_missing_column_fails_when_threshold_is_required() -> None:
    with pytest.raises(
        PhosPyInputError,
        match=("mode=require_threshold requires site_metadata.localisation_confidence"),
    ):
        _run_with_localisation(
            mode=DATASET_LOCALISATION_MODE_REQUIRE_THRESHOLD,
            site_metadata=_site_metadata(include_localisation_column=False),
        )


def test_localisation_below_threshold_fails_when_threshold_is_required() -> None:
    with pytest.raises(PhosPyInputError, match="localisation_confidence >= 0.750"):
        _run_with_localisation(
            mode=DATASET_LOCALISATION_MODE_REQUIRE_THRESHOLD,
            min_confidence=0.75,
            site_metadata=_site_metadata(localisation_values=[0.75, 0.7, 0.92]),
        )


def test_localisation_values_at_or_above_threshold_pass() -> None:
    preprocessed = _run_with_localisation(
        mode=DATASET_LOCALISATION_MODE_REQUIRE_THRESHOLD,
        min_confidence=0.75,
        site_metadata=_site_metadata(localisation_values=[0.75, 0.9, 1.0]),
    )
    assert preprocessed.phospho.index.tolist() == [
        "MAPK14;Y182;",
        "AKT1;T308;",
        "GSK3B;S9;",
    ]
    assert preprocessed.preprocessing_trace is not None
    assert preprocessed.preprocessing_trace[0].stage == "localisation_confidence"


@pytest.mark.parametrize("invalid_confidence", [-0.01, 1.01])
def test_localisation_invalid_confidence_values_fail(invalid_confidence: float) -> None:
    with pytest.raises(
        PhosPyInputError,
        match="found invalid values in site_metadata.localisation_confidence",
    ):
        _run_with_localisation(
            mode=DATASET_LOCALISATION_MODE_REQUIRE_THRESHOLD,
            site_metadata=_site_metadata(
                localisation_values=[0.95, invalid_confidence, 0.9]
            ),
        )


def test_localisation_waiver_mode_passes_and_records_diagnostics() -> None:
    preprocessed = _run_with_localisation(
        mode=DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER,
        min_confidence=0.75,
        waiver_reason="legacy import with unresolved probabilities",
        site_metadata=_site_metadata(localisation_values=[0.95, pd.NA, 0.4]),
    )
    assert preprocessed.row_audit is not None
    waived = preprocessed.row_audit.loc[
        preprocessed.row_audit.loc[:, "stage"] == "localisation_confidence"
    ]
    assert waived.shape[0] == 2
    assert set(waived.loc[:, "action"].astype(str)) == {"retained"}

    assert preprocessed.preprocessing_trace is not None
    diagnostics = preprocessed.preprocessing_trace[0].diagnostics
    assert diagnostics.get("mode") == "allow_missing_with_waiver"
    assert diagnostics.get("waiver_reason") == (
        "legacy import with unresolved probabilities"
    )
    assert diagnostics.get("missing_count") == 1
    assert diagnostics.get("below_threshold_count") == 1

    assert preprocessed.preprocessing_operations is not None
    localisation_ops = preprocessed.preprocessing_operations.loc[
        preprocessed.preprocessing_operations.loc[:, "stage"]
        == "localisation_confidence",
        :,
    ]
    assert localisation_ops.shape[0] == 1
    op_row = localisation_ops.iloc[0]
    assert op_row["operation"] == "allow_missing_with_waiver"
    assert (
        op_row["parameters"]["localisation_waiver_reason"]
        == "legacy import with unresolved probabilities"
    )


def test_localisation_waiver_policy_is_recorded_in_dataset_provenance() -> None:
    phospho = _phospho()
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(localisation_values=[0.95, pd.NA, 0.4]),
            preprocessing_config=DatasetPreprocessingConfig(
                localisation=DatasetLocalisationConfig(
                    mode=DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER,
                    min_confidence=0.75,
                    waiver_reason="legacy import with unresolved probabilities",
                )
            ),
        )
    )
    assert built.provenance is not None
    plan_payload = built.provenance.workflow_parameters["preprocessing_plan"]
    assert plan_payload["localisation_mode"] == "allow_missing_with_waiver"
    assert plan_payload["localisation_min_confidence"] == 0.75
    assert plan_payload["localisation_waiver_reason"] == (
        "legacy import with unresolved probabilities"
    )

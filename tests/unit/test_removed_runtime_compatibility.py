from __future__ import annotations

import warnings
from types import SimpleNamespace

import pandas as pd
import pytest

import phospy.differential.planning.replicates as replicate_models
import phospy.io.readers.enrichment_sets as enrichment_readers
from phospy.activities.semantics import ActivityInputMatrix
from phospy.errors import WorkflowBoundaryError
from phospy.kinase.activity.inputs import (
    KinaseActivityInputs,
    PredMatOverlapSummary,
)
from phospy.kinase.activity.models import ActivityMethodSummary
from phospy.kinase.activity.results import KinaseActivityResult
from phospy.kinase.contracts import ValidatedKinaseWorkflowRequest
from phospy.signalome.contracts import ValidatedSignalomeWorkflowRequest


def _activity_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["rat|P12345|S10"], name="site_key"),
    )


def test_removed_activity_input_constructor_forms_fail_without_warning() -> None:
    frame = _activity_frame()
    pred_mat = pd.DataFrame({"K1": [0.9]}, index=frame.index.copy())
    common = {
        "pred_mat": pred_mat,
        "phospho_matrix": frame,
        "threshold": 0.5,
        "min_substrates": 1,
        "top_n_substrates": 1,
        "overlap_summary": PredMatOverlapSummary(1, 1, 1),
    }

    with warnings.catch_warnings(record=True) as caught:
        with pytest.raises(TypeError, match="activity_input"):
            KinaseActivityInputs(**common)  # type: ignore[arg-type]
        with pytest.raises(TypeError, match="effect_matrix"):
            KinaseActivityInputs(
                **common,
                activity_input=ActivityInputMatrix.sample_level_abundance(frame),
                effect_matrix=frame,
            )  # type: ignore[call-arg]
    assert caught == []


def test_removed_activity_result_fallback_and_accessors_fail_without_warning() -> None:
    frame = _activity_frame().rename_axis(index="kinase")
    activity_input = ActivityInputMatrix.sample_level_abundance(frame)

    with warnings.catch_warnings(record=True) as caught:
        with pytest.raises(WorkflowBoundaryError, match="input_semantics"):
            KinaseActivityResult(activity_matrix=frame)
        result = KinaseActivityResult(
            activity_matrix=frame,
            input_semantics=activity_input.semantics,
            profile_metadata=activity_input.profile_metadata,
        )
        assert not hasattr(result, "activity_scores")
        assert not hasattr(result, "weighted_activity")
        assert not hasattr(result, "legacy_condition_statistics_table_dataframe")
    assert caught == []


def test_generated_activity_summary_condition_aliases_are_absent() -> None:
    summary = ActivityMethodSummary(
        kinases_evaluated=1,
        kinase_profile_pairs_evaluated=1,
        kinase_profile_pairs_computed=1,
        kinase_profile_pairs_insufficient_substrates=0,
        kinase_profile_pairs_invalid_background_variance=0,
        kinase_profile_pairs_no_finite_background_values=0,
        kinase_profile_pairs_no_finite_substrate_values=0,
    )

    for suffix in (
        "evaluated",
        "computed",
        "insufficient_substrates",
        "invalid_background_variance",
        "no_finite_background_values",
        "no_finite_substrate_values",
    ):
        assert not hasattr(summary, f"kinase_condition_pairs_{suffix}")


def test_removed_class_and_reader_aliases_are_absent() -> None:
    assert not hasattr(replicate_models, "TechnicalReplicateResolver")
    for name in (
        "load_enrichment_sets_gmt",
        "load_enrichment_sets_table",
        "load_enrichment_sets_csv",
        "load_enrichment_sets_tsv",
    ):
        assert not hasattr(enrichment_readers, name)


@pytest.mark.parametrize(
    "validated_type",
    (ValidatedKinaseWorkflowRequest, ValidatedSignalomeWorkflowRequest),
)
def test_validated_request_wrappers_do_not_delegate_unknown_attributes(
    validated_type: type,
) -> None:
    request = SimpleNamespace(legacy_internal_field="removed")
    validated = validated_type(request=request, dataset_view=object())

    with warnings.catch_warnings(record=True) as caught:
        with pytest.raises(AttributeError, match="legacy_internal_field"):
            _ = validated.legacy_internal_field
    assert caught == []

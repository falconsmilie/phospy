from __future__ import annotations

from collections import Counter
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import pytest

import phospy.science.datasets.internal_frame_store as internal_frame_store_module
import phospy.validation.workflows.method_quantitative as method_quantitative_module
import phospy.workflows.kinase.interpreter as interpreter_module
import phospy.workflows.kinase.validator as validator_module
from phospy.advanced import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api import (
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.api.results import KinasePredictionResult, KinaseScoringResult
from phospy.contracts.results import KinaseWorkflowResult
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseWorkflowRequest,
    ValidatedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    protein_site_key_index,
    site_key_context_columns,
)

_DISPLAY_IDS = ("MAPK14;Y182;", "GSK3B;S9;")
_GENE_SYMBOLS = ("MAPK14", "GSK3B")
_SITES = ("Y182", "S9")


@dataclass
class _SnapshotEvents:
    views: list[DatasetInternalView] = field(default_factory=list)
    dataframe_snapshots: Counter[str] = field(default_factory=Counter)
    optional_dataframe_snapshots: Counter[str] = field(default_factory=Counter)


@dataclass(slots=True)
class _FullMatrixCopyCounts:
    full_matrix_deep: int = 0


@contextmanager
def _count_full_matrix_deep_copies(
    *,
    shape: tuple[int, int],
    columns: tuple[object, ...],
) -> Iterator[_FullMatrixCopyCounts]:
    counts = _FullMatrixCopyCounts()
    original_copy = pd.DataFrame.copy

    def wrapped_copy(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        deep = kwargs.get("deep", args[0] if args else True)
        if bool(deep) and self.shape == shape and tuple(self.columns) == columns:
            counts.full_matrix_deep += 1
        return original_copy(self, *args, **kwargs)

    pd.DataFrame.copy = wrapped_copy
    try:
        yield counts
    finally:
        pd.DataFrame.copy = original_copy


@pytest.fixture
def instrument_dataset_view(
    monkeypatch: pytest.MonkeyPatch,
) -> _SnapshotEvents:
    events = _SnapshotEvents()
    original_view_class = DatasetInternalView
    original_dataframe_snapshot = (
        internal_frame_store_module.immutable_dataframe_snapshot
    )
    original_optional_dataframe_snapshot = (
        internal_frame_store_module.immutable_optional_dataframe_snapshot
    )

    class _CountingDatasetInternalView(original_view_class):
        def __init__(self, dataset: AnalysisReadyPhosphoDataset) -> None:
            events.views.append(self)
            super().__init__(dataset)

    def _counting_dataframe_snapshot(
        value: pd.DataFrame,
        *,
        field_name: str,
        error_type: type[Exception] = TypeError,
    ):
        events.dataframe_snapshots.update((field_name,))
        return original_dataframe_snapshot(
            value,
            field_name=field_name,
            error_type=error_type,
        )

    def _counting_optional_dataframe_snapshot(
        value: pd.DataFrame | None,
        *,
        field_name: str,
        error_type: type[Exception] = TypeError,
    ):
        if isinstance(value, pd.DataFrame):
            events.optional_dataframe_snapshots.update((field_name,))
        return original_optional_dataframe_snapshot(
            value,
            field_name=field_name,
            error_type=error_type,
        )

    monkeypatch.setattr(
        validator_module,
        "DatasetInternalView",
        _CountingDatasetInternalView,
    )
    monkeypatch.setattr(
        method_quantitative_module,
        "DatasetInternalView",
        _CountingDatasetInternalView,
    )
    monkeypatch.setattr(
        interpreter_module,
        "DatasetInternalView",
        _CountingDatasetInternalView,
    )
    monkeypatch.setattr(
        internal_frame_store_module,
        "immutable_dataframe_snapshot",
        _counting_dataframe_snapshot,
    )
    monkeypatch.setattr(
        internal_frame_store_module,
        "immutable_optional_dataframe_snapshot",
        _counting_optional_dataframe_snapshot,
    )
    return events


def _site_index() -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=_GENE_SYMBOLS,
        sites=_SITES,
    )


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0, 0.8], "sample_b": [2.0, 1.2]},
        index=_site_index(),
    )


def _site_metadata(*, object_payload: object | None = None) -> pd.DataFrame:
    site_index = _site_index()
    metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": list(_DISPLAY_IDS),
            **site_key_context_columns(site_index),
            "gene_symbol": list(_GENE_SYMBOLS),
            "protein_id": list(_GENE_SYMBOLS),
            "site": list(_SITES),
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            ],
            "localisation_confidence": [0.95, 0.90],
        },
        index=site_index,
    )
    if object_payload is not None:
        metadata.loc[:, "object_payload"] = pd.Series(
            [object_payload, _mutable_payload()],
            index=metadata.index,
            dtype=object,
        )
    return metadata


def _dataset(*, object_payload: object | None = None) -> AnalysisReadyPhosphoDataset:
    phospho = _phospho()
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=_site_metadata(object_payload=object_payload),
        sample_metadata=pd.DataFrame(
            {"batch": ["batch_1", "batch_1"]},
            index=phospho.columns.copy(),
        ),
        comparisons=pd.DataFrame(
            {"treated_vs_control": [0.25, -0.15]},
            index=phospho.index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": list(_DISPLAY_IDS),
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                    "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                ]
            },
            index=pd.Index(_DISPLAY_IDS, name="site_id"),
        ),
    )


def _request(
    dataset: AnalysisReadyPhosphoDataset | None = None,
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=_dataset() if dataset is None else dataset,
        references=_references(),
        scoring_config=KinaseScoringConfig(
            reliability_profile="custom",
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
    )


def _minimal_result(request: ResolvedKinaseWorkflowRequest) -> KinaseWorkflowResult:
    score_matrix = pd.DataFrame(
        {"MAP2K6": [1.0 for _ in request.scoring_site_index]},
        index=request.scoring_site_index.copy(),
    )
    return KinaseWorkflowResult(
        dataset=request.dataset,
        references=request.references,
        scoring_result=KinaseScoringResult._from_owned(
            profile_scores=score_matrix.copy(deep=True),
        ),
        prediction_result=KinasePredictionResult._from_owned(
            pred_mat=score_matrix.copy(deep=True),
        ),
        activity_result=None,
    )


def _mutable_payload() -> dict[str, object]:
    return {
        "list": ["start"],
        "dict": {"nested": ["start"]},
        "array": np.asarray([1.0, 2.0]),
    }


def _mutate_payload(payload: object) -> None:
    assert isinstance(payload, dict)
    list_value = payload["list"]
    dict_value = payload["dict"]
    array_value = payload["array"]
    assert isinstance(list_value, list)
    assert isinstance(dict_value, dict)
    assert isinstance(array_value, np.ndarray)
    nested_value = dict_value["nested"]
    assert isinstance(nested_value, list)
    list_value.append("mutated")
    nested_value.append("mutated")
    array_value[0] = 99.0


def _immutable_payload_state(payload: object) -> tuple[tuple[object, ...], ...]:
    assert isinstance(payload, Mapping)
    list_value = payload["list"]
    dict_value = payload["dict"]
    array_value = payload["array"]
    assert isinstance(list_value, tuple)
    assert isinstance(dict_value, Mapping)
    assert isinstance(array_value, np.ndarray)
    nested_value = dict_value["nested"]
    assert isinstance(nested_value, tuple)
    return (
        ("list", *list_value),
        ("dict", *nested_value),
        ("array", *tuple(float(value) for value in array_value.tolist())),
    )


def test_full_kinase_run_creates_one_dataset_view_and_one_required_snapshot(
    instrument_dataset_view: _SnapshotEvents,
) -> None:
    result = KinaseWorkflow().run(_request())

    assert not result.scoring_result.profile_scores.empty
    assert len(instrument_dataset_view.views) == 1
    assert (
        instrument_dataset_view.dataframe_snapshots["dataset.phospho internal snapshot"]
        == 1
    )
    assert (
        instrument_dataset_view.dataframe_snapshots[
            "dataset.site_metadata internal snapshot"
        ]
        == 1
    )
    assert instrument_dataset_view.optional_dataframe_snapshots == Counter(
        {
            "dataset.sample_metadata internal snapshot": 1,
            "dataset.comparisons internal snapshot": 1,
        }
    )


def test_validator_and_interpreter_share_one_view_per_run_but_not_across_runs() -> None:
    validator_views: list[DatasetInternalView] = []
    interpreter_views: list[DatasetInternalView] = []

    class _RecordingValidator(KinaseWorkflowValidator):
        def run(self, request: object) -> ValidatedKinaseWorkflowRequest:
            validated = super().run(request)
            validator_views.append(validated.dataset_view)
            return validated

    class _RecordingInterpreter(KinaseWorkflowInterpreter):
        def run(
            self, request: ValidatedKinaseWorkflowRequest
        ) -> ResolvedKinaseWorkflowRequest:
            interpreter_views.append(request.dataset_view)
            return super().run(request)

    request = _request()
    workflow = KinaseWorkflow._with_components(
        validator=_RecordingValidator(),
        interpreter=_RecordingInterpreter(),
    )

    workflow.run(request)
    workflow.run(request)

    assert validator_views[0] is interpreter_views[0]
    assert validator_views[1] is interpreter_views[1]
    assert validator_views[0] is not validator_views[1]


def test_repeated_view_access_within_run_reuses_required_snapshots(
    instrument_dataset_view: _SnapshotEvents,
) -> None:
    class _RepeatedAccessInterpreter(KinaseWorkflowInterpreter):
        def run(
            self, request: ValidatedKinaseWorkflowRequest
        ) -> ResolvedKinaseWorkflowRequest:
            _ = (
                request.dataset_view.phospho,
                request.dataset_view.phospho,
                request.dataset_view.site_metadata,
                request.dataset_view.site_metadata,
            )
            return super().run(request)

    KinaseWorkflow._with_components(interpreter=_RepeatedAccessInterpreter()).run(
        _request()
    )

    assert len(instrument_dataset_view.views) == 1
    assert (
        instrument_dataset_view.dataframe_snapshots["dataset.phospho internal snapshot"]
        == 1
    )
    assert (
        instrument_dataset_view.dataframe_snapshots[
            "dataset.site_metadata internal snapshot"
        ]
        == 1
    )


def test_repeated_kinase_runs_reuse_dataset_owned_snapshots(
    instrument_dataset_view: _SnapshotEvents,
) -> None:
    request = _request()
    workflow = KinaseWorkflow()

    workflow.run(request)
    workflow.run(request)

    assert len(instrument_dataset_view.views) == 2
    assert (
        instrument_dataset_view.dataframe_snapshots["dataset.phospho internal snapshot"]
        == 1
    )
    assert (
        instrument_dataset_view.dataframe_snapshots[
            "dataset.site_metadata internal snapshot"
        ]
        == 1
    )
    assert instrument_dataset_view.optional_dataframe_snapshots == Counter(
        {
            "dataset.sample_metadata internal snapshot": 1,
            "dataset.comparisons internal snapshot": 1,
        }
    )


def test_repeated_kinase_runs_do_not_repeat_full_dataset_matrix_deep_copies() -> None:
    request = _request()

    with _count_full_matrix_deep_copies(
        shape=request.dataset._phospho.shape,
        columns=tuple(request.dataset._phospho.columns),
    ) as counts:
        workflow = KinaseWorkflow()
        workflow.run(request)
        workflow.run(request)

    assert counts.full_matrix_deep == 1


def test_public_dataset_export_mutation_does_not_reach_workflow_inputs() -> None:
    dataset = _dataset()
    request = _request(dataset)
    original_phospho = dataset.phospho
    captured_activity_inputs: list[pd.DataFrame] = []
    captured_display_ids: list[list[str]] = []

    phospho_export = dataset.phospho
    phospho_export.iloc[0, 0] = 999.0
    site_metadata_export = dataset.site_metadata
    site_metadata_export.loc[site_metadata_export.index[0], "display_id"] = "BROKEN;S1;"

    class _CaptureInterpreter(KinaseWorkflowInterpreter):
        def run(
            self, request: ValidatedKinaseWorkflowRequest
        ) -> ResolvedKinaseWorkflowRequest:
            resolved = super().run(request)
            captured_activity_inputs.append(
                resolved.activity_phospho_matrix.copy(deep=True)
            )
            assert resolved.site_identity_map is not None
            captured_display_ids.append(
                resolved.site_identity_map.loc[:, "display_id"].astype(str).tolist()
            )
            return resolved

    class _MinimalExecutor:
        def run(self, request: ResolvedKinaseWorkflowRequest) -> KinaseWorkflowResult:
            return _minimal_result(request)

    KinaseWorkflow._with_components(
        interpreter=_CaptureInterpreter(),
        executor=_MinimalExecutor(),
    ).run(request)

    pd.testing.assert_frame_equal(captured_activity_inputs[0], original_phospho)
    assert captured_display_ids == [list(_DISPLAY_IDS)]
    pd.testing.assert_frame_equal(dataset.phospho, original_phospho)


def test_mutating_result_exports_does_not_alter_dataset_or_later_run() -> None:
    request = _request()
    result = KinaseWorkflow().run(request)
    original_dataset_phospho = request.dataset.phospho
    original_scores = result.scoring_result.profile_scores
    original_predictions = result.prediction_result.pred_mat

    exported_scores = result.scoring_result.profile_scores
    exported_predictions = result.prediction_result.pred_mat
    exported_scores.iloc[0, 0] = 999.0
    exported_predictions.iloc[0, 0] = 999.0

    pd.testing.assert_frame_equal(result.scoring_result.profile_scores, original_scores)
    pd.testing.assert_frame_equal(
        result.prediction_result.pred_mat,
        original_predictions,
    )
    pd.testing.assert_frame_equal(request.dataset.phospho, original_dataset_phospho)

    later_result = KinaseWorkflow().run(request)
    pd.testing.assert_frame_equal(
        later_result.scoring_result.profile_scores,
        original_scores,
    )
    pd.testing.assert_frame_equal(
        later_result.prediction_result.pred_mat,
        original_predictions,
    )


def test_object_bearing_site_metadata_is_isolated_by_reused_view() -> None:
    original_payload = _mutable_payload()
    dataset = _dataset(object_payload=original_payload)
    public_metadata_export = dataset.site_metadata
    exported_payload = public_metadata_export.loc[
        public_metadata_export.index[0],
        "object_payload",
    ]
    observed_states: list[tuple[tuple[object, ...], ...]] = []

    _mutate_payload(original_payload)
    _mutate_payload(exported_payload)

    class _PayloadCaptureInterpreter(KinaseWorkflowInterpreter):
        def run(
            self, request: ValidatedKinaseWorkflowRequest
        ) -> ResolvedKinaseWorkflowRequest:
            metadata = request.dataset_view.site_metadata
            payload = metadata.loc[metadata.index[0], "object_payload"]
            observed_states.append(_immutable_payload_state(payload))
            return super().run(request)

    class _MinimalExecutor:
        def run(self, request: ResolvedKinaseWorkflowRequest) -> KinaseWorkflowResult:
            return _minimal_result(request)

    KinaseWorkflow._with_components(
        interpreter=_PayloadCaptureInterpreter(),
        executor=_MinimalExecutor(),
    ).run(_request(dataset))

    assert observed_states == [
        (
            ("list", "start"),
            ("dict", "start"),
            ("array", 1.0, 2.0),
        )
    ]


def test_injected_validator_and_interpreter_use_updated_validated_protocol() -> None:
    calls: list[str] = []
    captured_view: DatasetInternalView | None = None

    class _Validator:
        def run(self, request: object) -> ValidatedKinaseWorkflowRequest:
            calls.append("validator")
            assert isinstance(request, KinaseWorkflowRequest)
            return ValidatedKinaseWorkflowRequest(
                request=request,
                dataset_view=DatasetInternalView(request.dataset),
            )

    class _Interpreter(KinaseWorkflowInterpreter):
        def run(
            self, request: ValidatedKinaseWorkflowRequest
        ) -> ResolvedKinaseWorkflowRequest:
            nonlocal captured_view
            calls.append("interpreter")
            captured_view = request.dataset_view
            return super().run(request)

    class _Executor:
        def run(self, request: ResolvedKinaseWorkflowRequest) -> KinaseWorkflowResult:
            calls.append("executor")
            return _minimal_result(request)

    result = KinaseWorkflow._with_components(
        validator=_Validator(),
        interpreter=_Interpreter(),
        executor=_Executor(),
    ).run(_request())

    assert isinstance(result, KinaseWorkflowResult)
    assert captured_view is not None
    assert calls == ["validator", "interpreter", "executor"]

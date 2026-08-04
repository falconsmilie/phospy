from __future__ import annotations

import hashlib
import json
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import phospy.io.bundles._signalome.reconstruction as reconstruction_module
import phospy.science.datasets.internal_view as internal_view_module
import phospy.science.signalomes._result_validation as result_validation_module
import phospy.workflows.signalome.result_assembly as result_assembly_module
import phospy.workflows.signalome.validator as validator_module
from phospy import KinaseWorkflow
from phospy.advanced import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api import (
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    SignalomeWorkflowRequest,
)
from phospy.api.results import KinasePredictionResult, KinaseScoringResult
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles.signalome import (
    SignalomeWorkflowConfigSnapshot,
    load_signalome_workflow_bundle,
    save_signalome_workflow_bundle,
)
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.signalomes.constants import DISPLAY_ID_COLUMN
from phospy.science.sites.site_keys import (
    ProteinScopedPhosphositeKey,
    encode_site_key,
)
from phospy.workflows.signalome.context_tables import SignalomeContextTableBuilder
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeWorkflowRequest,
    ValidatedSignalomeWorkflowRequest,
)
from phospy.workflows.signalome.executor import SignalomeWorkflowExecutor
from phospy.workflows.signalome.interpreter import SignalomeWorkflowInterpreter
from phospy.workflows.signalome.module_tables import SignalomeModuleTableBuilder
from phospy.workflows.signalome.protein_resolution import SignalomeProteinResolver
from phospy.workflows.signalome.public import SignalomeWorkflow
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import site_key_context_columns


@dataclass
class _SnapshotEvents:
    views: list[DatasetInternalView] = field(default_factory=list)
    dataframe_snapshots: Counter[str] = field(default_factory=Counter)
    optional_dataframe_snapshots: Counter[str] = field(default_factory=Counter)


@pytest.fixture
def instrument_signalome_dataset_view(
    monkeypatch: pytest.MonkeyPatch,
) -> _SnapshotEvents:
    events = _SnapshotEvents()
    original_view_class = DatasetInternalView
    original_dataframe_snapshot = internal_view_module.immutable_dataframe_snapshot
    original_optional_dataframe_snapshot = (
        internal_view_module.immutable_optional_dataframe_snapshot
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
        result_validation_module,
        "DatasetInternalView",
        _CountingDatasetInternalView,
    )
    monkeypatch.setattr(
        internal_view_module,
        "immutable_dataframe_snapshot",
        _counting_dataframe_snapshot,
    )
    monkeypatch.setattr(
        internal_view_module,
        "immutable_optional_dataframe_snapshot",
        _counting_optional_dataframe_snapshot,
    )
    return events


def _site_key(*, protein_id: str, site: str) -> str:
    return encode_site_key(
        ProteinScopedPhosphositeKey(
            organism=Organism.RAT.value,
            protein_namespace="protein_id",
            protein_identifier=protein_id,
            residue=site[0],
            position=int(site[1:]),
        )
    )


def _site_index() -> pd.Index:
    protein_ids = ("P1", "P1", "P2", "P3")
    sites = ("S1", "S2", "S3", "S4")
    return pd.Index(
        [
            _site_key(protein_id=protein_id, site=site)
            for protein_id, site in zip(protein_ids, sites, strict=True)
        ],
        name="site_key",
    )


def _display_ids() -> tuple[str, ...]:
    return ("P1;S1;", "P1;S2;", "P2;S3;", "P3;S4;")


def _protein_ids() -> tuple[str, ...]:
    return ("P1", "P1", "P2", "P3")


def _sites() -> tuple[str, ...]:
    return ("S1", "S2", "S3", "S4")


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0, 4.0],
            "sample_b": [1.2, 2.2, 3.2, 4.2],
        },
        index=_site_index(),
    )


def _site_metadata(*, object_payload: object | None = None) -> pd.DataFrame:
    site_index = _site_index()
    sites = _sites()
    metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": list(_display_ids()),
            **site_key_context_columns(site_index),
            "gene_symbol": list(_protein_ids()),
            "site": list(sites),
            "site_sequence": [
                ("A" * 15) + site[0].upper() + ("A" * 15) for site in sites
            ],
            "protein_id": list(_protein_ids()),
            "localisation_confidence": [0.95, 0.95, 0.95, 0.95],
        },
        index=site_index.copy(),
    )
    if object_payload is not None:
        metadata.loc[:, "object_payload"] = pd.Series(
            [
                object_payload,
                _mutable_payload(),
                _mutable_payload(),
                _mutable_payload(),
            ],
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
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references() -> ReferenceBundle:
    display_ids = list(_display_ids())
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K1", "K2", "K2"],
                "substrate_site": [
                    display_ids[0],
                    display_ids[2],
                    display_ids[1],
                    display_ids[3],
                ],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": [("A" * 15) + "S" + ("A" * 15) for _ in display_ids]},
            index=pd.Index(display_ids, name="site_id"),
        ),
    )


def _matrix(*, values: list[list[float]]) -> pd.DataFrame:
    return pd.DataFrame(
        values,
        index=_site_index().copy(),
        columns=pd.Index(["K1", "K2"], name="kinase"),
        dtype=float,
    )


def _kinase_result(dataset: AnalysisReadyPhosphoDataset) -> KinaseWorkflowResult:
    prediction_matrix = _matrix(
        values=[
            [0.9, 0.1],
            [0.1, 0.9],
            [0.8, 0.2],
            [0.7, 0.3],
        ]
    )
    score_matrix = _matrix(
        values=[
            [1.0, 0.0],
            [0.0, 1.0],
            [0.9, 0.1],
            [0.2, 0.8],
        ]
    )
    return KinaseWorkflowResult(
        dataset=dataset,
        references=_references(),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            rank_weighted_fusion_scores=score_matrix,
        ),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


def _request(
    dataset: AnalysisReadyPhosphoDataset | None = None,
) -> SignalomeWorkflowRequest:
    resolved_dataset = _dataset() if dataset is None else dataset
    return SignalomeWorkflowRequest(
        kinase_result=_kinase_result(resolved_dataset),
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            network_correlation_threshold=0.3,
            network_min_paired_finite_observations=3,
            module_count=2,
        ),
    )


def _provenanced_request() -> SignalomeWorkflowRequest:
    dataset = _dataset()
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
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
    )
    return SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            network_correlation_threshold=0.3,
            network_min_paired_finite_observations=3,
            module_count=2,
        ),
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


def _sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _refresh_signalome_table_entry(
    *,
    bundle_root: Path,
    manifest: dict[str, object],
    table_key: str,
    table: pd.DataFrame,
) -> None:
    signalome_outputs = manifest["signalome_outputs"]
    assert isinstance(signalome_outputs, dict)
    tables = signalome_outputs["tables"]
    assert isinstance(tables, dict)
    entry = tables[table_key]
    assert isinstance(entry, dict)
    relative_path = entry["path"]
    assert isinstance(relative_path, str)
    path = bundle_root / relative_path
    entry["byte_size"] = path.stat().st_size
    entry["sha256"] = _sha256_path(path)
    entry["shape"] = {"rows": int(table.shape[0]), "columns": int(table.shape[1])}


def test_full_signalome_run_creates_one_workflow_view_and_required_snapshots(
    instrument_signalome_dataset_view: _SnapshotEvents,
) -> None:
    result = SignalomeWorkflow().run(_request())

    assert not result.module_assignments.table.empty
    assert len(instrument_signalome_dataset_view.views) == 1
    assert (
        instrument_signalome_dataset_view.dataframe_snapshots[
            "dataset.phospho internal snapshot"
        ]
        == 1
    )
    assert (
        instrument_signalome_dataset_view.dataframe_snapshots[
            "dataset.site_metadata internal snapshot"
        ]
        == 1
    )
    assert instrument_signalome_dataset_view.optional_dataframe_snapshots == Counter()


def test_signalome_view_identity_is_threaded_through_interpretation_and_execution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    observed: dict[str, object] = {}

    class _RecordingValidator(SignalomeWorkflowValidator):
        def run(
            self, request: SignalomeWorkflowRequest
        ) -> ValidatedSignalomeWorkflowRequest:
            validated = super().run(request)
            observed["validator_view"] = validated.dataset_view
            return validated

    class _RecordingProteinResolver(SignalomeProteinResolver):
        def run(
            self,
            *,
            site_metadata: pd.DataFrame,
            site_index: pd.Index,
            removed_by_score_preconditioning_count: int,
        ) -> pd.Series:
            observed["protein_resolver_site_metadata_columns"] = tuple(
                site_metadata.columns.astype(str).tolist()
            )
            return super().run(
                site_metadata=site_metadata,
                site_index=site_index,
                removed_by_score_preconditioning_count=(
                    removed_by_score_preconditioning_count
                ),
            )

    class _RecordingInterpreter(SignalomeWorkflowInterpreter):
        def run(
            self, request: ValidatedSignalomeWorkflowRequest
        ) -> ResolvedSignalomeWorkflowRequest:
            observed["interpreter_view"] = request.dataset_view
            resolved = super().run(request)
            observed["resolved_view"] = resolved.dataset_view
            return resolved

    class _RecordingModuleTableBuilder(SignalomeModuleTableBuilder):
        def run(self, **kwargs: object):
            request = kwargs["request"]
            assert isinstance(request, ResolvedSignalomeWorkflowRequest)
            observed["module_builder_view"] = request.dataset_view
            return super().run(**kwargs)

    class _RecordingContextTableBuilder(SignalomeContextTableBuilder):
        def run(self, **kwargs: object):
            request = kwargs["request"]
            assert isinstance(request, ResolvedSignalomeWorkflowRequest)
            observed["context_builder_view"] = request.dataset_view
            return super().run(**kwargs)

    original_identity_validator = (
        result_assembly_module.validate_signalome_result_site_level_identity
    )

    def _recording_identity_validator(**kwargs: object) -> None:
        site_metadata = kwargs["site_metadata"]
        assert isinstance(site_metadata, pd.DataFrame)
        observed["result_identity_site_metadata_columns"] = tuple(
            site_metadata.columns.astype(str).tolist()
        )
        original_identity_validator(**kwargs)

    monkeypatch.setattr(
        result_assembly_module,
        "validate_signalome_result_site_level_identity",
        _recording_identity_validator,
    )

    workflow = SignalomeWorkflow._with_components(
        validator=_RecordingValidator(),
        interpreter=_RecordingInterpreter(protein_resolver=_RecordingProteinResolver()),
        executor=SignalomeWorkflowExecutor(
            module_table_builder=_RecordingModuleTableBuilder(),
            context_table_builder=_RecordingContextTableBuilder(),
        ),
    )

    result = workflow.run(_request())

    assert not result.signalome_modules.table.empty
    view = observed["validator_view"]
    assert observed["interpreter_view"] is view
    assert observed["resolved_view"] is view
    assert observed["module_builder_view"] is view
    assert observed["context_builder_view"] is view
    assert "site_key" in observed["protein_resolver_site_metadata_columns"]
    assert "site_key" in observed["result_identity_site_metadata_columns"]


def test_signalome_independent_runs_do_not_share_dataset_views() -> None:
    validator_views: list[DatasetInternalView] = []
    interpreter_views: list[DatasetInternalView] = []

    class _RecordingValidator(SignalomeWorkflowValidator):
        def run(
            self, request: SignalomeWorkflowRequest
        ) -> ValidatedSignalomeWorkflowRequest:
            validated = super().run(request)
            validator_views.append(validated.dataset_view)
            return validated

    class _RecordingInterpreter(SignalomeWorkflowInterpreter):
        def run(
            self, request: ValidatedSignalomeWorkflowRequest
        ) -> ResolvedSignalomeWorkflowRequest:
            interpreter_views.append(request.dataset_view)
            return super().run(request)

    request = _request()
    workflow = SignalomeWorkflow._with_components(
        validator=_RecordingValidator(),
        interpreter=_RecordingInterpreter(),
    )

    workflow.run(request)
    workflow.run(request)

    assert validator_views[0] is interpreter_views[0]
    assert validator_views[1] is interpreter_views[1]
    assert validator_views[0] is not validator_views[1]


def test_signalome_public_exports_are_mutation_isolated_across_runs() -> None:
    dataset = _dataset()
    request = _request(dataset)
    original_phospho = dataset.phospho
    original_site_metadata = dataset.site_metadata

    phospho_export = dataset.phospho
    phospho_export.iloc[0, 0] = 999.0
    metadata_export = dataset.site_metadata
    metadata_export.loc[metadata_export.index[0], DISPLAY_ID_COLUMN] = "BROKEN;S1;"

    result = SignalomeWorkflow().run(request)
    assignments = result.module_assignments.table
    assert "BROKEN;S1;" not in set(assignments.loc[:, DISPLAY_ID_COLUMN].astype(str))
    pd.testing.assert_frame_equal(dataset.phospho, original_phospho)
    pd.testing.assert_frame_equal(dataset.site_metadata, original_site_metadata)

    exported_assignments = result.module_assignments.table
    exported_assignments.loc[exported_assignments.index[0], DISPLAY_ID_COLUMN] = (
        "BROKEN;S1;"
    )

    assert "BROKEN;S1;" not in set(
        result.module_assignments.table.loc[:, DISPLAY_ID_COLUMN].astype(str)
    )
    later_result = SignalomeWorkflow().run(request)
    assert "BROKEN;S1;" not in set(
        later_result.module_assignments.table.loc[:, DISPLAY_ID_COLUMN].astype(str)
    )


def test_signalome_object_bearing_site_metadata_is_isolated_by_reused_view() -> None:
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

    class _PayloadCaptureInterpreter(SignalomeWorkflowInterpreter):
        def run(
            self, request: ValidatedSignalomeWorkflowRequest
        ) -> ResolvedSignalomeWorkflowRequest:
            metadata = request.dataset_view.site_metadata
            payload = metadata.loc[metadata.index[0], "object_payload"]
            observed_states.append(_immutable_payload_state(payload))
            return super().run(request)

    SignalomeWorkflow._with_components(interpreter=_PayloadCaptureInterpreter()).run(
        _request(dataset)
    )

    assert observed_states == [
        (
            ("list", "start"),
            ("dict", "start"),
            ("array", 1.0, 2.0),
        )
    ]


def test_standalone_bundle_reconstruction_validates_identity_with_isolated_snapshot(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    request = _provenanced_request()
    result = SignalomeWorkflow().run(request)
    bundle_root = tmp_path / "signalome_bundle"
    save_signalome_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=SignalomeWorkflowConfigSnapshot.from_request(request),
        output_format="csv",
    )

    manifest_path = bundle_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    signalome_outputs = manifest["signalome_outputs"]
    assert isinstance(signalome_outputs, dict)
    tables = signalome_outputs["tables"]
    assert isinstance(tables, dict)
    site_membership_entry = tables["site_membership"]
    assert isinstance(site_membership_entry, dict)
    site_membership_path = bundle_root / str(site_membership_entry["path"])
    site_membership = pd.read_csv(site_membership_path, index_col=0)
    site_membership.loc[site_membership.index[0], DISPLAY_ID_COLUMN] = "BROKEN;S1;"
    site_membership.to_csv(site_membership_path)
    _refresh_signalome_table_entry(
        bundle_root=bundle_root,
        manifest=manifest,
        table_key="site_membership",
        table=site_membership,
    )
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    views: list[DatasetInternalView] = []
    original_view_class = DatasetInternalView

    class _CountingDatasetInternalView(original_view_class):
        def __init__(self, dataset: AnalysisReadyPhosphoDataset) -> None:
            views.append(self)
            super().__init__(dataset)

    monkeypatch.setattr(
        reconstruction_module,
        "DatasetInternalView",
        _CountingDatasetInternalView,
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "bundle signalome result identity validation failed: "
            "signalome_result.site_membership.display_id values must match"
        ),
    ):
        load_signalome_workflow_bundle(bundle_root)

    assert len(views) == 1

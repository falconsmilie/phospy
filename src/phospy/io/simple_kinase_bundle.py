"""External output-bundle services for simple kinase workflow results."""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from phospy.activities.models import KinaseActivityResult
from phospy.api.configs import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
)
from phospy.api.results import SimpleKinaseWorkflowResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.input import PhosPyInputError
from phospy.io.tables import read_table, table_suffix_for_format, write_table
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.references.models import Organism, ReferenceBundle
from phospy.transformations.models import (
    MatrixTransformationState,
    TransformationKind,
    TransformationState,
)

if TYPE_CHECKING:
    from phospy.api.requests import SimpleKinaseWorkflowRequest

SIMPLE_KINASE_BUNDLE_MANIFEST_VERSION = 1

_SIMPLE_KINASE_BUNDLE_KIND = "simple_kinase_workflow_result"
_MANIFEST_FILENAME = "manifest.json"
_CONFIG_SNAPSHOT_RELATIVE_PATH = "config/snapshot.json"


@dataclass(frozen=True, slots=True)
class SimpleKinaseWorkflowConfigSnapshot:
    """Serializable snapshot of the simple kinase workflow configuration."""

    scoring_config: KinaseScoringConfig
    prediction_config: KinasePredictionConfig
    activity_config: KinaseActivityConfig | None

    @classmethod
    def from_request(
        cls, request: SimpleKinaseWorkflowRequest
    ) -> SimpleKinaseWorkflowConfigSnapshot:
        """Create a config snapshot from a workflow request."""

        from phospy.api.requests import SimpleKinaseWorkflowRequest

        if not isinstance(request, SimpleKinaseWorkflowRequest):
            raise TypeError("request must be a SimpleKinaseWorkflowRequest")
        return cls(
            scoring_config=request.scoring_config,
            prediction_config=request.prediction_config,
            activity_config=request.activity_config,
        )

    def to_payload(self) -> dict[str, object]:
        """Return a manifest-safe JSON payload for this config snapshot."""

        activity_payload: dict[str, object] | None
        if self.activity_config is None:
            activity_payload = None
        else:
            activity_payload = {
                "enabled": self.activity_config.enabled,
                "threshold": float(self.activity_config.threshold),
            }
        return {
            "scoring_config": {
                "min_substrates": self.scoring_config.min_substrates,
            },
            "prediction_config": {
                "top_k": self.prediction_config.top_k,
                "ensemble_size": self.prediction_config.ensemble_size,
            },
            "activity_config": activity_payload,
        }

    @classmethod
    def from_payload(
        cls, payload: Mapping[str, object]
    ) -> SimpleKinaseWorkflowConfigSnapshot:
        """Create a config snapshot from a decoded JSON payload."""

        scope = "config snapshot"
        scoring_payload = _require_mapping(
            payload.get("scoring_config"),
            field_name=f"{scope}.scoring_config",
        )
        prediction_payload = _require_mapping(
            payload.get("prediction_config"),
            field_name=f"{scope}.prediction_config",
        )
        activity_raw = payload.get("activity_config")
        if activity_raw is None:
            activity_config = None
        else:
            activity_payload = _require_mapping(
                activity_raw,
                field_name=f"{scope}.activity_config",
            )
            activity_config = KinaseActivityConfig(
                enabled=_require_bool(
                    activity_payload.get("enabled"),
                    field_name=f"{scope}.activity_config.enabled",
                ),
                threshold=_require_float(
                    activity_payload.get("threshold"),
                    field_name=f"{scope}.activity_config.threshold",
                ),
            )
        return cls(
            scoring_config=KinaseScoringConfig(
                min_substrates=_require_int(
                    scoring_payload.get("min_substrates"),
                    field_name=f"{scope}.scoring_config.min_substrates",
                )
            ),
            prediction_config=KinasePredictionConfig(
                top_k=_require_int(
                    prediction_payload.get("top_k"),
                    field_name=f"{scope}.prediction_config.top_k",
                ),
                ensemble_size=_require_int(
                    prediction_payload.get("ensemble_size"),
                    field_name=f"{scope}.prediction_config.ensemble_size",
                ),
            ),
            activity_config=activity_config,
        )


@dataclass(frozen=True, slots=True)
class LoadedSimpleKinaseWorkflowBundle:
    """Loaded simple kinase output bundle contents."""

    result: SimpleKinaseWorkflowResult
    config_snapshot: SimpleKinaseWorkflowConfigSnapshot
    manifest_version: int


def save_simple_kinase_workflow_bundle(
    result: SimpleKinaseWorkflowResult,
    output_root: Path,
    *,
    config_snapshot: SimpleKinaseWorkflowConfigSnapshot,
    output_format: str = "csv",
) -> dict[str, Path]:
    """Write a reproducible simple kinase output bundle and return written paths."""

    bundle_root = Path(output_root)
    suffix = table_suffix_for_format(output_format)
    normalized_format = output_format.strip().lower()
    written: dict[str, Path] = {}

    dataset_tables = {
        "phospho": _write_bundle_table(
            table=result.dataset.phospho,
            bundle_root=bundle_root,
            relative_path=Path("dataset") / f"phospho{suffix}",
            written=written,
            written_key="dataset.phospho",
        ),
        "site_metadata": _write_bundle_table(
            table=result.dataset.site_metadata,
            bundle_root=bundle_root,
            relative_path=Path("dataset") / f"site_metadata{suffix}",
            written=written,
            written_key="dataset.site_metadata",
        ),
        "sample_metadata": _write_optional_bundle_table(
            table=result.dataset.sample_metadata,
            bundle_root=bundle_root,
            relative_path=Path("dataset") / f"sample_metadata{suffix}",
            written=written,
            written_key="dataset.sample_metadata",
        ),
        "total": _write_optional_bundle_table(
            table=result.dataset.total,
            bundle_root=bundle_root,
            relative_path=Path("dataset") / f"total{suffix}",
            written=written,
            written_key="dataset.total",
        ),
    }

    reference_tables = {
        "kinase_substrate_map": _write_bundle_table(
            table=result.references.kinase_substrate_map,
            bundle_root=bundle_root,
            relative_path=Path("references") / f"kinase_substrate_map{suffix}",
            written=written,
            written_key="references.kinase_substrate_map",
        ),
        "site_sequences": _write_bundle_table(
            table=result.references.site_sequences,
            bundle_root=bundle_root,
            relative_path=Path("references") / f"site_sequences{suffix}",
            written=written,
            written_key="references.site_sequences",
        ),
    }

    scoring_tables = {
        "profile_scores": _write_bundle_table(
            table=result.scoring_result.profile_scores,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"profile_scores{suffix}",
            written=written,
            written_key="scoring.profile_scores",
        ),
        "motif_scores": _write_optional_bundle_table(
            table=result.scoring_result.motif_scores,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"motif_scores{suffix}",
            written=written,
            written_key="scoring.motif_scores",
        ),
        "combined_scores": _write_optional_bundle_table(
            table=result.scoring_result.combined_scores,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"combined_scores{suffix}",
            written=written,
            written_key="scoring.combined_scores",
        ),
        "weights": _write_optional_bundle_table(
            table=result.scoring_result.weights,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"weights{suffix}",
            written=written,
            written_key="scoring.weights",
        ),
    }

    prediction_tables = {
        "pred_mat": _write_bundle_table(
            table=result.prediction_result.pred_mat,
            bundle_root=bundle_root,
            relative_path=Path("prediction") / f"pred_mat{suffix}",
            written=written,
            written_key="prediction.pred_mat",
        ),
        "substrate_list": _write_optional_bundle_table(
            table=result.prediction_result.substrate_list,
            bundle_root=bundle_root,
            relative_path=Path("prediction") / f"substrate_list{suffix}",
            written=written,
            written_key="prediction.substrate_list",
        ),
    }

    activity_tables = {
        "activity_scores": _write_optional_bundle_table(
            table=(
                None
                if result.activity_result is None
                else result.activity_result.activity_scores
            ),
            bundle_root=bundle_root,
            relative_path=Path("activity") / f"activity_scores{suffix}",
            written=written,
            written_key="activity.activity_scores",
        )
    }

    config_path = bundle_root / Path(_CONFIG_SNAPSHOT_RELATIVE_PATH)
    _write_json(config_path, config_snapshot.to_payload(), label="config snapshot")
    written["config_snapshot"] = config_path

    manifest = {
        "bundle_type": _SIMPLE_KINASE_BUNDLE_KIND,
        "manifest_version": SIMPLE_KINASE_BUNDLE_MANIFEST_VERSION,
        "table_format": normalized_format,
        "dataset": {
            "metadata": {
                "organism": (
                    None
                    if result.dataset.organism is None
                    else result.dataset.organism.value
                ),
                "transformation_state": _transformation_state_to_payload(
                    result.dataset.transformation_state
                ),
            },
            "tables": dataset_tables,
        },
        "resolved_references": {
            "metadata": {
                "organism": result.references.organism.value,
            },
            "tables": reference_tables,
        },
        "outputs": {
            "scoring": {
                "tables": scoring_tables,
            },
            "prediction": {
                "tables": prediction_tables,
            },
            "activity": {
                "enabled": result.activity_result is not None,
                "tables": activity_tables,
            },
        },
        "config_snapshot": _CONFIG_SNAPSHOT_RELATIVE_PATH,
    }
    manifest_path = bundle_root / _MANIFEST_FILENAME
    _write_json(manifest_path, manifest, label="bundle manifest")
    written["manifest"] = manifest_path
    return written


def load_simple_kinase_workflow_bundle(
    bundle_root: Path,
) -> LoadedSimpleKinaseWorkflowBundle:
    """Load a simple kinase output bundle from disk."""

    root = Path(bundle_root)
    manifest = _read_json(root / _MANIFEST_FILENAME, label="bundle manifest")
    manifest_payload = _require_mapping(manifest, field_name="bundle manifest")

    bundle_type = _require_str(
        manifest_payload.get("bundle_type"),
        field_name="bundle manifest.bundle_type",
    )
    if bundle_type != _SIMPLE_KINASE_BUNDLE_KIND:
        raise PhosPyInputError(
            "unsupported bundle manifest bundle_type "
            f"'{bundle_type}'; expected '{_SIMPLE_KINASE_BUNDLE_KIND}'"
        )

    manifest_version = _require_int(
        manifest_payload.get("manifest_version"),
        field_name="bundle manifest.manifest_version",
    )
    if manifest_version != SIMPLE_KINASE_BUNDLE_MANIFEST_VERSION:
        raise PhosPyInputError(
            "unsupported bundle manifest version "
            f"'{manifest_version}'; expected {SIMPLE_KINASE_BUNDLE_MANIFEST_VERSION}"
        )

    dataset_payload = _require_mapping(
        manifest_payload.get("dataset"),
        field_name="bundle manifest.dataset",
    )
    dataset_metadata = _require_mapping(
        dataset_payload.get("metadata"),
        field_name="bundle manifest.dataset.metadata",
    )
    dataset_tables = _require_mapping(
        dataset_payload.get("tables"),
        field_name="bundle manifest.dataset.tables",
    )

    references_payload = _require_mapping(
        manifest_payload.get("resolved_references"),
        field_name="bundle manifest.resolved_references",
    )
    references_metadata = _require_mapping(
        references_payload.get("metadata"),
        field_name="bundle manifest.resolved_references.metadata",
    )
    reference_tables = _require_mapping(
        references_payload.get("tables"),
        field_name="bundle manifest.resolved_references.tables",
    )

    outputs_payload = _require_mapping(
        manifest_payload.get("outputs"),
        field_name="bundle manifest.outputs",
    )
    scoring_payload = _require_mapping(
        outputs_payload.get("scoring"),
        field_name="bundle manifest.outputs.scoring",
    )
    scoring_tables = _require_mapping(
        scoring_payload.get("tables"),
        field_name="bundle manifest.outputs.scoring.tables",
    )
    prediction_payload = _require_mapping(
        outputs_payload.get("prediction"),
        field_name="bundle manifest.outputs.prediction",
    )
    prediction_tables = _require_mapping(
        prediction_payload.get("tables"),
        field_name="bundle manifest.outputs.prediction.tables",
    )
    activity_payload = _require_mapping(
        outputs_payload.get("activity"),
        field_name="bundle manifest.outputs.activity",
    )
    activity_tables = _require_mapping(
        activity_payload.get("tables"),
        field_name="bundle manifest.outputs.activity.tables",
    )

    dataset = AnalysisReadyPhosphoDataset(
        phospho=_read_required_table(
            bundle_root=root,
            tables=dataset_tables,
            table_key="phospho",
            field_name="bundle manifest.dataset.tables.phospho",
        ),
        site_metadata=_read_required_table(
            bundle_root=root,
            tables=dataset_tables,
            table_key="site_metadata",
            field_name="bundle manifest.dataset.tables.site_metadata",
        ),
        sample_metadata=_read_optional_table(
            bundle_root=root,
            tables=dataset_tables,
            table_key="sample_metadata",
            field_name="bundle manifest.dataset.tables.sample_metadata",
        ),
        total=_read_optional_table(
            bundle_root=root,
            tables=dataset_tables,
            table_key="total",
            field_name="bundle manifest.dataset.tables.total",
        ),
        organism=_parse_optional_organism(
            dataset_metadata.get("organism"),
            field_name="bundle manifest.dataset.metadata.organism",
        ),
        transformation_state=_transformation_state_from_payload(
            _require_mapping(
                dataset_metadata.get("transformation_state"),
                field_name="bundle manifest.dataset.metadata.transformation_state",
            )
        ),
    )

    references = ReferenceBundle(
        organism=_parse_required_organism(
            references_metadata.get("organism"),
            field_name="bundle manifest.resolved_references.metadata.organism",
        ),
        kinase_substrate_map=_read_required_table(
            bundle_root=root,
            tables=reference_tables,
            table_key="kinase_substrate_map",
            field_name="bundle manifest.resolved_references.tables.kinase_substrate_map",
        ),
        site_sequences=_read_required_table(
            bundle_root=root,
            tables=reference_tables,
            table_key="site_sequences",
            field_name="bundle manifest.resolved_references.tables.site_sequences",
        ),
    )

    scoring_result = KinaseScoringResult(
        profile_scores=_read_required_table(
            bundle_root=root,
            tables=scoring_tables,
            table_key="profile_scores",
            field_name="bundle manifest.outputs.scoring.tables.profile_scores",
        ),
        motif_scores=_read_optional_table(
            bundle_root=root,
            tables=scoring_tables,
            table_key="motif_scores",
            field_name="bundle manifest.outputs.scoring.tables.motif_scores",
        ),
        combined_scores=_read_optional_table(
            bundle_root=root,
            tables=scoring_tables,
            table_key="combined_scores",
            field_name="bundle manifest.outputs.scoring.tables.combined_scores",
        ),
        weights=_read_optional_table(
            bundle_root=root,
            tables=scoring_tables,
            table_key="weights",
            field_name="bundle manifest.outputs.scoring.tables.weights",
        ),
    )

    prediction_result = KinasePredictionResult(
        pred_mat=_read_required_table(
            bundle_root=root,
            tables=prediction_tables,
            table_key="pred_mat",
            field_name="bundle manifest.outputs.prediction.tables.pred_mat",
        ),
        substrate_list=_read_optional_table(
            bundle_root=root,
            tables=prediction_tables,
            table_key="substrate_list",
            field_name="bundle manifest.outputs.prediction.tables.substrate_list",
        ),
    )

    activity_table = _read_optional_table(
        bundle_root=root,
        tables=activity_tables,
        table_key="activity_scores",
        field_name="bundle manifest.outputs.activity.tables.activity_scores",
    )
    activity_result = (
        None
        if activity_table is None
        else KinaseActivityResult(activity_scores=activity_table)
    )

    config_snapshot_path = _resolve_bundle_relative_path(
        root,
        _require_str(
            manifest_payload.get("config_snapshot"),
            field_name="bundle manifest.config_snapshot",
        ),
        field_name="bundle manifest.config_snapshot",
    )
    config_snapshot_payload = _require_mapping(
        _read_json(config_snapshot_path, label="config snapshot"),
        field_name="config snapshot",
    )
    config_snapshot = SimpleKinaseWorkflowConfigSnapshot.from_payload(
        config_snapshot_payload
    )

    return LoadedSimpleKinaseWorkflowBundle(
        result=SimpleKinaseWorkflowResult(
            dataset=dataset,
            references=references,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            activity_result=activity_result,
        ),
        config_snapshot=config_snapshot,
        manifest_version=manifest_version,
    )


def _write_bundle_table(
    *,
    table,
    bundle_root: Path,
    relative_path: Path,
    written: dict[str, Path],
    written_key: str,
) -> str:
    output_path = bundle_root / relative_path
    write_table(table, output_path)
    written[written_key] = output_path
    return relative_path.as_posix()


def _write_optional_bundle_table(
    *,
    table,
    bundle_root: Path,
    relative_path: Path,
    written: dict[str, Path],
    written_key: str,
) -> str | None:
    if table is None:
        return None
    return _write_bundle_table(
        table=table,
        bundle_root=bundle_root,
        relative_path=relative_path,
        written=written,
        written_key=written_key,
    )


def _transformation_state_to_payload(state: TransformationState) -> dict[str, object]:
    return {
        "phospho": _matrix_state_to_payload(state.phospho),
        "total": (
            None if state.total is None else _matrix_state_to_payload(state.total)
        ),
    }


def _matrix_state_to_payload(state: MatrixTransformationState) -> dict[str, object]:
    return {
        "kind": state.kind.value,
        "transformed": state.transformed,
        "established_by": state.established_by,
    }


def _transformation_state_from_payload(
    payload: Mapping[str, object],
) -> TransformationState:
    phospho_payload = _require_mapping(
        payload.get("phospho"),
        field_name="dataset.metadata.transformation_state.phospho",
    )
    total_raw = payload.get("total")
    if total_raw is None:
        total_state = None
    else:
        total_state = _matrix_state_from_payload(
            _require_mapping(
                total_raw,
                field_name="dataset.metadata.transformation_state.total",
            )
        )
    return TransformationState(
        phospho=_matrix_state_from_payload(phospho_payload),
        total=total_state,
    )


def _matrix_state_from_payload(
    payload: Mapping[str, object],
) -> MatrixTransformationState:
    kind_token = _require_str(
        payload.get("kind"),
        field_name="matrix_transformation_state.kind",
    )
    try:
        kind = TransformationKind(kind_token)
    except ValueError as exc:
        supported = ", ".join(member.value for member in TransformationKind)
        raise PhosPyInputError(
            f"unsupported transformation kind '{kind_token}'; supported: {supported}"
        ) from exc
    return MatrixTransformationState(
        kind=kind,
        transformed=_require_bool(
            payload.get("transformed"),
            field_name="matrix_transformation_state.transformed",
        ),
        established_by=_require_str(
            payload.get("established_by"),
            field_name="matrix_transformation_state.established_by",
        ),
    )


def _read_required_table(
    *,
    bundle_root: Path,
    tables: Mapping[str, object],
    table_key: str,
    field_name: str,
):
    table_path = _resolve_bundle_relative_path(
        bundle_root,
        _require_str(tables.get(table_key), field_name=field_name),
        field_name=field_name,
    )
    return read_table(table_path)


def _read_optional_table(
    *,
    bundle_root: Path,
    tables: Mapping[str, object],
    table_key: str,
    field_name: str,
):
    if table_key not in tables:
        raise PhosPyInputError(f"{field_name} must be declared in the bundle manifest")
    raw_value = tables[table_key]
    if raw_value is None:
        return None
    table_path = _resolve_bundle_relative_path(
        bundle_root,
        _require_str(raw_value, field_name=field_name),
        field_name=field_name,
    )
    return read_table(table_path)


def _parse_optional_organism(value: object, *, field_name: str) -> Organism | None:
    if value is None:
        return None
    return _parse_required_organism(value, field_name=field_name)


def _parse_required_organism(value: object, *, field_name: str) -> Organism:
    token = _require_str(value, field_name=field_name)
    try:
        return Organism(token)
    except ValueError as exc:
        supported = ", ".join(member.value for member in Organism)
        raise PhosPyInputError(
            f"unsupported organism '{token}' in {field_name}; supported: {supported}"
        ) from exc


def _resolve_bundle_relative_path(
    bundle_root: Path,
    relative_path: str,
    *,
    field_name: str,
) -> Path:
    candidate = Path(relative_path)
    if candidate.is_absolute():
        raise PhosPyInputError(f"{field_name} must be a relative bundle path")
    resolved_root = bundle_root.resolve()
    resolved_candidate = (bundle_root / candidate).resolve()
    try:
        resolved_candidate.relative_to(resolved_root)
    except ValueError as exc:
        raise PhosPyInputError(
            f"{field_name} points outside the bundle root: {relative_path}"
        ) from exc
    return resolved_candidate


def _read_json(path: Path, *, label: str) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise PhosPyInputError(f"{label} does not exist: {path}") from exc
    except PermissionError as exc:
        raise PhosPyInputError(
            f"permission denied while reading {label}: {path}"
        ) from exc
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise PhosPyInputError(f"failed to parse {label} '{path}': {exc}") from exc


def _write_json(path: Path, payload: object, *, label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True),
            encoding="utf-8",
        )
    except (OSError, TypeError, ValueError) as exc:
        raise PhosPyInputError(f"failed to write {label} '{path}': {exc}") from exc


def _require_mapping(value: object, *, field_name: str) -> Mapping[str, object]:
    if isinstance(value, Mapping):
        return value
    raise PhosPyInputError(f"{field_name} must be an object")


def _require_str(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return normalized


def _require_bool(value: object, *, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be a bool")
    return value


def _require_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be an int")
    return value


def _require_float(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PhosPyInputError(f"{field_name} must be a float")
    return float(value)


__all__ = [
    "LoadedSimpleKinaseWorkflowBundle",
    "SIMPLE_KINASE_BUNDLE_MANIFEST_VERSION",
    "SimpleKinaseWorkflowConfigSnapshot",
    "load_simple_kinase_workflow_bundle",
    "save_simple_kinase_workflow_bundle",
]

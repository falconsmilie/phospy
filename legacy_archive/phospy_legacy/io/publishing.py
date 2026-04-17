from __future__ import annotations

import json
import shutil
from collections.abc import Callable, Mapping
from dataclasses import asdict
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from json import JSONDecodeError
from pathlib import Path
from typing import TYPE_CHECKING
from uuid import uuid4

from ..internal.constants import RUN_MANIFEST_FILENAME

if TYPE_CHECKING:
    from ..activities.results import KinaseActivityResult
    from ..api import DatasetLoadOptions, KinaseActivityConfig, PredictionRunConfig
    from ..api.contracts import (
        SimpleKinaseWorkflowBundleMetadata,
        SimpleKinaseWorkflowConfigSnapshot,
    )
    from ..api.workflow_results import SimpleKinaseWorkflowResult
    from ..io.readers import SimpleKinaseWorkflowOutputBundle
    from ..preprocessing import CorePreprocessingConfig, CoreProcessingResult


def package_version() -> str:
    try:
        return version("phospy")
    except PackageNotFoundError:
        return "unknown"


class RunManifestWriter:
    """Serialize pipeline execution metadata to a JSON manifest."""

    def __init__(
        self,
        *,
        package_version_resolver: Callable[[], str] | None = None,
        clock: Callable[[], datetime] | None = None,
    ) -> None:
        self._package_version_resolver = package_version_resolver or package_version
        self._clock = clock or self._utc_now

    def write(
        self,
        *,
        outdir: Path,
        core: CoreProcessingResult,
        kinase_activity: KinaseActivityResult | None,
        preprocessing_config: CorePreprocessingConfig,
    ) -> None:
        manifest = {
            "status": "success",
            "generated_at_utc": self._clock().isoformat(),
            "package_version": self._package_version_resolver(),
            "has_kinase_activity": kinase_activity is not None,
            "core_rows": {
                "total_unique": int(core.total_unique.shape[0]),
                "total_filtered": int(core.total_filtered.shape[0]),
                "phospho_filtered": int(core.phospho_filtered.shape[0]),
                "phospho_corrected": int(core.phospho_corrected.shape[0]),
                "site_matrix": int(core.site_matrix.matrix.shape[0]),
            },
            "preprocessing_config": asdict(preprocessing_config),
        }
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / RUN_MANIFEST_FILENAME).write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @staticmethod
    def _utc_now() -> datetime:
        return datetime.now(timezone.utc)


class OutputPublisher:
    """Publish staged output directories via crash-recoverable replacement.

    Replacement is not atomic when the target directory already exists. The
    publisher renames the current target to a backup directory, then renames the
    staging directory into place. If the process is interrupted between those
    steps, the target may be temporarily absent. To make that state recoverable,
    the publisher writes a small recovery marker and restores or cleans up the
    previous state on the next publish attempt when possible.
    """

    def publish(self, *, staging_dir: Path, target_dir: Path) -> None:
        target_dir = self._resolve_path(target_dir)
        staging_dir = self._resolve_path(staging_dir)
        self._recover_if_needed(target_dir)
        if not target_dir.exists():
            self._replace_directory(staging_dir, target_dir)
            return

        backup_dir = self._backup_dir_for(target_dir)
        recovery_marker = self._recovery_marker_for(target_dir)
        self._write_recovery_marker(
            marker_path=recovery_marker,
            target_dir=target_dir,
            backup_dir=backup_dir,
        )
        self._replace_directory(target_dir, backup_dir)
        try:
            self._replace_directory(staging_dir, target_dir)
        except OSError:
            self._replace_directory(backup_dir, target_dir)
            self._remove_recovery_marker(recovery_marker)
            raise
        else:
            self._remove_directory(backup_dir)
            self._remove_recovery_marker(recovery_marker)

    def _recover_if_needed(self, target_dir: Path) -> None:
        target_dir = self._resolve_path(target_dir)
        marker_path = self._recovery_marker_for(target_dir)
        if not marker_path.exists():
            return

        try:
            recovery_state = self._read_recovery_marker(marker_path)
        except (JSONDecodeError, OSError, TypeError, ValueError):
            self._quarantine_recovery_marker(marker_path)
            return

        backup_dir = self._validated_backup_dir(
            marker_path=marker_path,
            target_dir=target_dir,
            recovery_state=recovery_state,
        )

        if not target_dir.exists() and backup_dir.exists():
            self._replace_directory(backup_dir, target_dir)
            self._remove_recovery_marker(marker_path)
            return

        if target_dir.exists() and backup_dir.exists():
            self._remove_directory(backup_dir)
            self._remove_recovery_marker(marker_path)
            return

        if target_dir.exists() and not backup_dir.exists():
            self._remove_recovery_marker(marker_path)
            return

        if not target_dir.exists() and not backup_dir.exists():
            self._remove_recovery_marker(marker_path)
            return

    @staticmethod
    def _backup_dir_for(target_dir: Path) -> Path:
        return target_dir.with_name(f".{target_dir.name}.backup-{uuid4().hex}")

    @staticmethod
    def _recovery_marker_for(target_dir: Path) -> Path:
        return target_dir.with_name(f".{target_dir.name}.publish-state.json")

    @staticmethod
    def _write_recovery_marker(
        *,
        marker_path: Path,
        target_dir: Path,
        backup_dir: Path,
    ) -> None:
        marker_path.write_text(
            json.dumps(
                {
                    "target_dir": str(OutputPublisher._resolve_path(target_dir)),
                    "backup_dir": str(OutputPublisher._resolve_path(backup_dir)),
                    "created_at_utc": datetime.now(timezone.utc).isoformat(),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )

    @staticmethod
    def _read_recovery_marker(marker_path: Path) -> dict[str, str]:
        raw_state = json.loads(marker_path.read_text(encoding="utf-8"))
        if not isinstance(raw_state, dict):
            raise ValueError(
                f"Recovery marker {marker_path} must contain a JSON object."
            )

        target_dir = raw_state.get("target_dir")
        backup_dir = raw_state.get("backup_dir")
        if not isinstance(target_dir, str) or not target_dir:
            raise ValueError(
                f"Recovery marker {marker_path} is missing a valid 'target_dir'."
            )
        if not isinstance(backup_dir, str) or not backup_dir:
            raise ValueError(
                f"Recovery marker {marker_path} is missing a valid 'backup_dir'."
            )

        return {
            "target_dir": target_dir,
            "backup_dir": backup_dir,
        }

    @staticmethod
    def _validated_backup_dir(
        *,
        marker_path: Path,
        target_dir: Path,
        recovery_state: dict[str, str],
    ) -> Path:
        resolved_target_dir = OutputPublisher._resolve_path(target_dir)
        expected_target_dir = OutputPublisher._resolve_path(
            Path(recovery_state["target_dir"])
        )
        backup_dir = OutputPublisher._resolve_path(Path(recovery_state["backup_dir"]))

        if expected_target_dir != resolved_target_dir:
            raise OSError(
                f"Recovery marker {marker_path} does not match target directory {resolved_target_dir}."
            )

        if backup_dir.parent != resolved_target_dir.parent:
            raise OSError(
                f"Recovery marker {marker_path} points outside target parent directory {resolved_target_dir.parent}."
            )

        expected_backup_prefix = f".{resolved_target_dir.name}.backup-"
        if not backup_dir.name.startswith(expected_backup_prefix):
            raise OSError(
                f"Recovery marker {marker_path} has unexpected backup directory name {backup_dir.name}."
            )

        return backup_dir

    @staticmethod
    def _quarantine_recovery_marker(marker_path: Path) -> Path:
        quarantined_path = marker_path.with_suffix(marker_path.suffix + ".corrupt")
        counter = 1
        while quarantined_path.exists():
            quarantined_path = marker_path.with_suffix(
                marker_path.suffix + f".corrupt.{counter}"
            )
            counter += 1
        marker_path.replace(quarantined_path)
        return quarantined_path

    @staticmethod
    def _remove_recovery_marker(marker_path: Path) -> None:
        marker_path.unlink(missing_ok=True)

    @staticmethod
    def _replace_directory(source: Path, target: Path) -> None:
        source.replace(target)

    @staticmethod
    def _remove_directory(target: Path) -> None:
        shutil.rmtree(target)

    @staticmethod
    def _resolve_path(path: Path) -> Path:
        return path.resolve(strict=False)


def save_simple_kinase_workflow_output_bundle(
    *,
    result: SimpleKinaseWorkflowResult,
    outdir: str | Path,
    config_snapshot: SimpleKinaseWorkflowConfigSnapshot
    | Mapping[str, object]
    | None = None,
    dataset_options: DatasetLoadOptions | Mapping[str, object] | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    prediction_config: PredictionRunConfig | Mapping[str, object] | None = None,
    activity_config: KinaseActivityConfig | Mapping[str, object] | None = None,
) -> Path:
    from ..api.contracts import SimpleKinaseWorkflowConfigSnapshot
    from .writers import SimpleKinaseWorkflowBundleWriter

    resolved_snapshot = config_snapshot
    if resolved_snapshot is None:
        resolved_snapshot = SimpleKinaseWorkflowConfigSnapshot.from_workflow_inputs(
            dataset_options=dataset_options,
            preprocessing_config=preprocessing_config,
            prediction_config=prediction_config,
            activity_config=activity_config,
        )

    return SimpleKinaseWorkflowBundleWriter().write(
        result=result,
        outdir=outdir,
        config_snapshot=resolved_snapshot,
    )


def load_simple_kinase_workflow_output_bundle_metadata(
    bundle_dir: str | Path,
) -> SimpleKinaseWorkflowBundleMetadata:
    from .readers import load_simple_kinase_workflow_output_bundle_metadata

    return load_simple_kinase_workflow_output_bundle_metadata(bundle_dir)


def load_simple_kinase_workflow_output_bundle(
    bundle_dir: str | Path,
    *,
    table_ids: tuple[str, ...] | list[str] | None = None,
) -> SimpleKinaseWorkflowOutputBundle:
    from .readers import load_simple_kinase_workflow_output_bundle

    return load_simple_kinase_workflow_output_bundle(
        bundle_dir,
        table_ids=table_ids,
    )

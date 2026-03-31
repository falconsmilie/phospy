from __future__ import annotations

import json
import shutil
from collections.abc import Callable
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING
from uuid import uuid4

if TYPE_CHECKING:
    from .analysis import KinaseActivityResult
    from .core_processing import CorePreprocessingConfig, CoreProcessingResult


def package_version() -> str:
    try:
        from importlib.metadata import version

        return version("phospy")
    except Exception:
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
        (outdir / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @staticmethod
    def _utc_now() -> datetime:
        return datetime.now(timezone.utc)


class OutputPublisher:
    """Publish staged output directories to their final location atomically."""

    def publish(self, *, staging_dir: Path, target_dir: Path) -> None:
        if not target_dir.exists():
            self._replace_directory(staging_dir, target_dir)
            return

        backup_dir = self._backup_dir_for(target_dir)
        self._replace_directory(target_dir, backup_dir)
        try:
            self._replace_directory(staging_dir, target_dir)
        except Exception:
            self._replace_directory(backup_dir, target_dir)
            raise
        else:
            self._remove_directory(backup_dir)

    @staticmethod
    def _backup_dir_for(target_dir: Path) -> Path:
        return target_dir.with_name(f".{target_dir.name}.backup-{uuid4().hex}")

    @staticmethod
    def _replace_directory(source: Path, target: Path) -> None:
        source.replace(target)

    @staticmethod
    def _remove_directory(target: Path) -> None:
        shutil.rmtree(target)

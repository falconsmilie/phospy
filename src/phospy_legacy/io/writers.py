from __future__ import annotations

import json
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Literal, Protocol, TypeAlias

import pandas as pd

if TYPE_CHECKING:
    from ..activities.results import KinaseActivityResult
    from ..api.contracts import (
        SimpleKinaseWorkflowConfigSnapshot,
        WorkflowOutputInventoryItem,
    )
    from ..api.workflow_results import SimpleKinaseWorkflowResult

from ..internal.constants import (
    CENTRALIZED_SEQUENCE_COLUMN,
    CORE_PHOSPHO_CORRECTED_BASENAME,
    CORE_PHOSPHO_FILTERED_BASENAME,
    CORE_PHOSR_INPUT_BASENAME,
    CORE_SITE_MATRIX_BASENAME,
    CORE_SITE_SEQUENCES_BASENAME,
    CORE_TOTAL_FILTERED_BASENAME,
    CORE_TOTAL_UNIQUE_BASENAME,
    KINASE_ACTIVITY_MATRIX_FILENAME,
    KINASE_TARGET_COUNTS_FILENAME,
    KINASE_TARGET_TABLE_FILENAME,
    KSEA_COUNTS_FILENAME,
    KSEA_SCORES_FILENAME,
    SIMPLE_KINASE_WORKFLOW_BUNDLE_FORMAT,
    SIMPLE_KINASE_WORKFLOW_RESULT_TYPE,
    WORKFLOW_OUTPUT_BUNDLE_MANIFEST_FILENAME,
)
from ..preprocessing import CoreProcessingResult
from .publishing import package_version

CoreOutputFormat: TypeAlias = Literal["csv", "tsv", "parquet"]


class TabularOutputWriter(Protocol):
    """Format-specific writer for pandas tabular outputs."""

    file_extension: str

    def write_table(
        self,
        table: pd.DataFrame,
        destination: Path,
        *,
        include_index: bool,
    ) -> None: ...


class CoreProcessingResultWriter(Protocol):
    """Writer contract for persisted core preprocessing outputs."""

    def write(
        self,
        result: CoreProcessingResult,
        outdir: str | Path,
        *,
        format: CoreOutputFormat | str = "csv",
    ) -> None: ...


class KinaseActivityResultWriter(Protocol):
    """Writer contract for persisted kinase activity outputs."""

    def write(self, result: KinaseActivityResult, outdir: str | Path) -> None: ...


@dataclass(frozen=True, slots=True)
class DelimitedTabularWriter:
    """Persist tabular outputs as a delimited text file."""

    separator: str
    file_extension: str

    def write_table(
        self,
        table: pd.DataFrame,
        destination: Path,
        *,
        include_index: bool,
    ) -> None:
        table.to_csv(destination, index=include_index, sep=self.separator)


@dataclass(frozen=True, slots=True)
class ParquetTabularWriter:
    """Persist tabular outputs as parquet files."""

    file_extension: str = ".parquet"

    def write_table(
        self,
        table: pd.DataFrame,
        destination: Path,
        *,
        include_index: bool,
    ) -> None:
        try:
            table.to_parquet(destination, index=include_index)
        except (ImportError, ModuleNotFoundError) as exc:
            msg = (
                "Parquet output requires an optional pandas parquet engine such as "
                "'pyarrow' or 'fastparquet'."
            )
            raise RuntimeError(msg) from exc


@dataclass(frozen=True, slots=True)
class CoreOutputArtifact:
    """Logical core preprocessing output and its persistence metadata."""

    basename: str
    table: pd.DataFrame
    include_index: bool = False


def _default_tabular_writers() -> dict[str, TabularOutputWriter]:
    return {
        "csv": DelimitedTabularWriter(separator=",", file_extension=".csv"),
        "tsv": DelimitedTabularWriter(separator="\t", file_extension=".tsv"),
        "parquet": ParquetTabularWriter(),
    }


@dataclass(frozen=True, slots=True)
class CoreOutputWriter:
    """Persist core preprocessing outputs using an explicit tabular format."""

    tabular_writers: Mapping[str, TabularOutputWriter] = field(
        default_factory=_default_tabular_writers
    )

    def write(
        self,
        result: CoreProcessingResult,
        outdir: str | Path,
        *,
        format: CoreOutputFormat | str = "csv",
    ) -> None:
        target_dir = Path(outdir)
        target_dir.mkdir(parents=True, exist_ok=True)
        tabular_writer = self._resolve_tabular_writer(format)

        for artifact in self._core_output_artifacts(result):
            tabular_writer.write_table(
                artifact.table,
                target_dir / f"{artifact.basename}{tabular_writer.file_extension}",
                include_index=artifact.include_index,
            )

    def _resolve_tabular_writer(
        self,
        format: CoreOutputFormat | str,
    ) -> TabularOutputWriter:
        normalized_format = format.strip().lower()
        try:
            return self.tabular_writers[normalized_format]
        except KeyError as exc:
            supported_formats = ", ".join(sorted(self.tabular_writers))
            msg = (
                f"Unsupported core output format {format!r}. "
                f"Supported formats are: {supported_formats}."
            )
            raise ValueError(msg) from exc

    @staticmethod
    def _core_output_artifacts(
        result: CoreProcessingResult,
    ) -> tuple[CoreOutputArtifact, ...]:
        return (
            CoreOutputArtifact(CORE_TOTAL_UNIQUE_BASENAME, result.total_unique),
            CoreOutputArtifact(CORE_TOTAL_FILTERED_BASENAME, result.total_filtered),
            CoreOutputArtifact(CORE_PHOSPHO_FILTERED_BASENAME, result.phospho_filtered),
            CoreOutputArtifact(
                CORE_PHOSPHO_CORRECTED_BASENAME, result.phospho_corrected
            ),
            CoreOutputArtifact(
                CORE_PHOSR_INPUT_BASENAME, result.site_matrix.phosr_input
            ),
            CoreOutputArtifact(
                CORE_SITE_MATRIX_BASENAME,
                result.site_matrix.matrix,
                include_index=True,
            ),
            CoreOutputArtifact(
                CORE_SITE_SEQUENCES_BASENAME,
                result.site_matrix.sequences.rename(
                    CENTRALIZED_SEQUENCE_COLUMN
                ).to_frame(),
                include_index=True,
            ),
        )


@dataclass(frozen=True, slots=True)
class KinaseActivityWriter:
    """Persist downstream kinase activity summaries to disk."""

    def write(self, result: KinaseActivityResult, outdir: str | Path) -> None:
        target_dir = Path(outdir)
        target_dir.mkdir(parents=True, exist_ok=True)
        result.weighted_activity.to_csv(target_dir / KINASE_ACTIVITY_MATRIX_FILENAME)
        result.ksea_scores.to_csv(target_dir / KSEA_SCORES_FILENAME)
        result.ksea_counts.to_csv(target_dir / KSEA_COUNTS_FILENAME)
        result.target_counts.to_csv(target_dir / KINASE_TARGET_COUNTS_FILENAME)
        result.target_table.to_csv(
            target_dir / KINASE_TARGET_TABLE_FILENAME,
            index=False,
        )


BundleTableValueType: TypeAlias = Literal["dataframe", "series"]


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


@dataclass(frozen=True, slots=True)
class BundleTableArtifact:
    """Logical output-bundle table artifact and persistence metadata."""

    table_id: str
    table: pd.DataFrame | pd.Series
    value_type: BundleTableValueType = "dataframe"
    include_index: bool = True


@dataclass(frozen=True, slots=True)
class SimpleKinaseWorkflowBundleWriter:
    """Persist one ``SimpleKinaseWorkflowResult`` as a reproducible bundle."""

    tabular_writer: TabularOutputWriter = field(
        default_factory=lambda: DelimitedTabularWriter(
            separator=",",
            file_extension=".csv",
        )
    )
    package_version_resolver: Callable[[], str] = package_version
    clock: Callable[[], datetime] = _utc_now

    def write(
        self,
        *,
        result: SimpleKinaseWorkflowResult,
        outdir: str | Path,
        config_snapshot: SimpleKinaseWorkflowConfigSnapshot
        | Mapping[str, object]
        | None = None,
    ) -> Path:
        from ..api.contracts import (
            SimpleKinaseWorkflowBundleMetadata,
            WorkflowOutputInventoryItem,
        )

        bundle_dir = Path(outdir)
        tables_dir = bundle_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)

        inventory: list[WorkflowOutputInventoryItem] = []
        for artifact in self._bundle_table_artifacts(result):
            destination = tables_dir / f"{artifact.table_id}.csv"
            frame = (
                artifact.table.to_frame(name=artifact.table.name)
                if isinstance(artifact.table, pd.Series)
                else artifact.table
            )
            self.tabular_writer.write_table(
                frame,
                destination,
                include_index=artifact.include_index,
            )
            inventory.append(
                WorkflowOutputInventoryItem(
                    table_id=artifact.table_id,
                    path=str(destination.relative_to(bundle_dir).as_posix()),
                    value_type=artifact.value_type,
                )
            )

        metadata = SimpleKinaseWorkflowBundleMetadata(
            workflow_type=SIMPLE_KINASE_WORKFLOW_RESULT_TYPE,
            bundle_format=SIMPLE_KINASE_WORKFLOW_BUNDLE_FORMAT,
            generated_at_utc=self.clock().isoformat(),
            package_version=str(self.package_version_resolver()),
            config_snapshot=self._resolve_config_snapshot(config_snapshot),
            reference_identity={
                "species": result.reference_bundle.species,
                "source": result.reference_bundle.source_metadata.source,
                "reference": result.reference_bundle.source_metadata.reference,
                "version": result.reference_bundle.source_metadata.version,
                "provider": result.reference_bundle.provenance.provider,
                "notes": list(result.reference_bundle.provenance.notes),
            },
            output_inventory=tuple(inventory),
        )
        (bundle_dir / WORKFLOW_OUTPUT_BUNDLE_MANIFEST_FILENAME).write_text(
            json.dumps(metadata.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return bundle_dir

    @staticmethod
    def _resolve_config_snapshot(
        config_snapshot: SimpleKinaseWorkflowConfigSnapshot
        | Mapping[str, object]
        | None,
    ) -> dict[str, object]:
        from ..api.contracts import SimpleKinaseWorkflowConfigSnapshot

        if config_snapshot is None:
            return SimpleKinaseWorkflowConfigSnapshot.from_workflow_inputs().to_dict()
        if isinstance(config_snapshot, SimpleKinaseWorkflowConfigSnapshot):
            return config_snapshot.to_dict()
        return dict(config_snapshot)

    @staticmethod
    def _bundle_table_artifacts(
        result: SimpleKinaseWorkflowResult,
    ) -> tuple[BundleTableArtifact, ...]:
        scoring_result = result.to_owned_scoring_result()
        artifacts: list[BundleTableArtifact] = [
            BundleTableArtifact(
                table_id="analysis_ready_phospho_matrix",
                table=result.analysis_ready_dataset.to_owned_phospho_matrix(),
                value_type="dataframe",
            ),
            BundleTableArtifact(
                table_id="analysis_ready_site_metadata",
                table=result.analysis_ready_dataset.to_owned_site_metadata(),
                value_type="dataframe",
            ),
            BundleTableArtifact(
                table_id="analysis_ready_site_sequences",
                table=result.analysis_ready_dataset.to_owned_site_sequences(),
                value_type="series",
            ),
            BundleTableArtifact(
                table_id="analysis_ready_phospho_corrected",
                table=result.analysis_ready_dataset.to_owned_phospho_corrected(),
                value_type="dataframe",
            ),
            BundleTableArtifact(
                table_id="pred_mat",
                table=result.pred_mat_result.to_owned_frame(),
                value_type="dataframe",
            ),
            BundleTableArtifact(
                table_id="profile_scores",
                table=scoring_result.profile_scores,
                value_type="dataframe",
            ),
            BundleTableArtifact(
                table_id="kinase_activity_matrix",
                table=result.kinase_activity_result.weighted_activity,
                value_type="dataframe",
            ),
            BundleTableArtifact(
                table_id="ksea_scores",
                table=result.kinase_activity_result.ksea_scores,
                value_type="dataframe",
            ),
            BundleTableArtifact(
                table_id="ksea_counts",
                table=result.kinase_activity_result.ksea_counts,
                value_type="series",
            ),
            BundleTableArtifact(
                table_id="kinase_target_counts",
                table=result.kinase_activity_result.target_counts,
                value_type="series",
            ),
            BundleTableArtifact(
                table_id="kinase_target_table",
                table=result.kinase_activity_result.target_table,
                value_type="dataframe",
            ),
        ]

        if scoring_result.combined_scores is not None:
            artifacts.append(
                BundleTableArtifact(
                    table_id="combined_scores",
                    table=scoring_result.combined_scores,
                    value_type="dataframe",
                )
            )
        if scoring_result.weights is not None:
            artifacts.append(
                BundleTableArtifact(
                    table_id="scoring_weights",
                    table=scoring_result.weights,
                    value_type="dataframe",
                )
            )
        return tuple(artifacts)

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Literal, Protocol, TypeAlias

import pandas as pd

if TYPE_CHECKING:
    from .activities import KinaseActivityResult

from .constants import (
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
)
from .preprocessing import CoreProcessingResult

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

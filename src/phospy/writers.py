from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Literal, Protocol, TypeAlias

import pandas as pd

if TYPE_CHECKING:
    from .analysis import KinaseActivityResult
    from .core_processing import CoreProcessingResult

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

    @staticmethod
    def write(result: KinaseActivityResult, outdir: str | Path) -> None: ...


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
            CoreOutputArtifact("df_total_unique", result.total_unique),
            CoreOutputArtifact("df_total_filtered", result.total_filtered),
            CoreOutputArtifact("df_phospho_filtered", result.phospho_filtered),
            CoreOutputArtifact("df_phospho_corrected", result.phospho_corrected),
            CoreOutputArtifact("phosr_input", result.site_matrix.phosr_input),
            CoreOutputArtifact(
                "mat_phospho_corrected",
                result.site_matrix.matrix,
                include_index=True,
            ),
            CoreOutputArtifact(
                "site_sequences",
                result.site_matrix.sequences.rename("centralized_sequence").to_frame(),
                include_index=True,
            ),
        )


class KinaseActivityWriter:
    """Persist downstream kinase activity summaries to disk."""

    @staticmethod
    def write(result: KinaseActivityResult, outdir: str | Path) -> None:
        target_dir = Path(outdir)
        target_dir.mkdir(parents=True, exist_ok=True)
        result.weighted_activity.to_csv(target_dir / "kinase_activity_matrix.csv")
        result.ksea_scores.to_csv(target_dir / "ksea_scores.csv")
        result.ksea_counts.to_csv(target_dir / "ksea_counts.csv")
        result.target_counts.to_csv(target_dir / "kinase_target_counts.csv")
        result.target_table.to_csv(
            target_dir / "kinase_target_table.csv",
            index=False,
        )

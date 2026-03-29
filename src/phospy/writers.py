from __future__ import annotations

from pathlib import Path

from .analysis import KinaseActivityResult
from .dataset import CoreProcessingResult


class CoreProcessingWriter:
    """Persist core preprocessing outputs to disk."""

    @staticmethod
    def write(result: CoreProcessingResult, outdir: str | Path) -> None:
        target_dir = Path(outdir)
        target_dir.mkdir(parents=True, exist_ok=True)
        result.total_unique.to_csv(target_dir / "df_total_unique.csv", index=False)
        result.total_filtered.to_csv(target_dir / "df_total_filtered.csv", index=False)
        result.phospho_filtered.to_csv(
            target_dir / "df_phospho_filtered.csv",
            index=False,
        )
        result.phospho_corrected.to_csv(
            target_dir / "df_phospho_corrected.csv",
            index=False,
        )
        result.site_matrix.phosr_input.to_csv(
            target_dir / "phosr_input.csv",
            index=False,
        )
        result.site_matrix.matrix.to_csv(target_dir / "mat_phospho_corrected.csv")
        result.site_matrix.sequences.rename("centralized_sequence").to_csv(
            target_dir / "site_sequences.csv"
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

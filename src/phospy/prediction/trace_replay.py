from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..validation.errors import TableSchemaError
from .models import SamplingTraceOverrideEnsemble


def _read_trace_table_from_directory(
    trace_dir: Path,
    table_name: str,
) -> tuple[pd.DataFrame | None, str | None]:
    csv_path = trace_dir / f"{table_name}.csv"
    if csv_path.exists():
        return pd.read_csv(csv_path), csv_path.name

    parquet_path = trace_dir / f"{table_name}.parquet"
    parquet_parts = sorted(trace_dir.glob(f"{table_name}.part-*.parquet"))
    if parquet_path.exists() or parquet_parts:
        parquet_sources = [parquet_path] if parquet_path.exists() else parquet_parts
        try:
            frame = pd.concat(
                [pd.read_parquet(source) for source in parquet_sources],
                ignore_index=True,
            )
        except Exception as error:  # pragma: no cover - engine availability varies
            msg = (
                "Unable to read parquet trace replay input. Install a supported "
                "parquet engine such as 'pyarrow', or provide CSV trace files."
            )
            raise RuntimeError(msg) from error
        label = (
            parquet_path.name
            if parquet_path.exists()
            else f"{table_name}.part-*.parquet"
        )
        return frame, label

    return None, None


class PredictionSamplingTrace:
    def __init__(
        self, ensembles_by_kinase: dict[str, dict[int, SamplingTraceOverrideEnsemble]]
    ) -> None:
        self.ensembles_by_kinase = ensembles_by_kinase

    @classmethod
    def from_trace_directory(cls, trace_dir: str | Path) -> PredictionSamplingTrace:
        path = Path(trace_dir)
        initial_df, initial_label = _read_trace_table_from_directory(
            path, "trace_initial_negatives"
        )
        samples_df, samples_label = _read_trace_table_from_directory(
            path, "trace_iteration_samples"
        )
        if initial_df is None and samples_df is None:
            msg = (
                "sampling trace directory must contain trace_initial_negatives "
                "and/or trace_iteration_samples in CSV or parquet format"
            )
            raise TableSchemaError(msg)

        ensembles_by_kinase: dict[str, dict[int, SamplingTraceOverrideEnsemble]] = {}

        if initial_df is not None:
            required_initial_cols = {"kinase", "ensemble", "draw", "site"}
            if not required_initial_cols.issubset(initial_df.columns):
                missing = sorted(required_initial_cols.difference(initial_df.columns))
                msg = f"{initial_label} is missing required columns: " + ", ".join(
                    missing
                )
                raise TableSchemaError(msg)
            initial_df = initial_df.sort_values(
                ["kinase", "ensemble", "draw"], kind="mergesort"
            )
            for (kinase, ensemble), group in initial_df.groupby(
                ["kinase", "ensemble"], sort=False
            ):
                ensemble_map = ensembles_by_kinase.setdefault(str(kinase), {})
                ensemble_map[int(ensemble)] = SamplingTraceOverrideEnsemble(
                    initial_negative_sites=group.loc[:, "site"].astype(str).tolist(),
                    iteration_sample_sites={},
                )

        if samples_df is not None:
            required_sample_cols = {
                "kinase",
                "ensemble",
                "iteration",
                "class_label",
                "draw",
                "site",
            }
            if not required_sample_cols.issubset(samples_df.columns):
                missing = sorted(required_sample_cols.difference(samples_df.columns))
                msg = f"{samples_label} is missing required columns: " + ", ".join(
                    missing
                )
                raise TableSchemaError(msg)
            samples_df = samples_df.sort_values(
                ["kinase", "ensemble", "iteration", "class_label", "draw"],
                kind="mergesort",
            )
            for (kinase, ensemble, iteration, class_label), group in samples_df.groupby(
                ["kinase", "ensemble", "iteration", "class_label"], sort=False
            ):
                ensemble_map = ensembles_by_kinase.setdefault(str(kinase), {})
                ensemble_override = ensemble_map.setdefault(
                    int(ensemble),
                    SamplingTraceOverrideEnsemble(
                        initial_negative_sites=None,
                        iteration_sample_sites={},
                    ),
                )
                iteration_map = ensemble_override.iteration_sample_sites.setdefault(
                    int(iteration), {}
                )
                iteration_map[int(class_label)] = (
                    group.loc[:, "site"].astype(str).tolist()
                )

        return cls(ensembles_by_kinase=ensembles_by_kinase)

    def get_ensemble_override(
        self, kinase: str, ensemble_index: int
    ) -> SamplingTraceOverrideEnsemble | None:
        return self.ensembles_by_kinase.get(kinase, {}).get(ensemble_index)

    def subset_kinases(self, kinases: list[str] | set[str]) -> PredictionSamplingTrace:
        kinase_set = {str(kinase) for kinase in kinases}
        return PredictionSamplingTrace(
            ensembles_by_kinase={
                kinase: ensemble_map
                for kinase, ensemble_map in self.ensembles_by_kinase.items()
                if kinase in kinase_set
            }
        )


__all__ = [
    "PredictionSamplingTrace",
]

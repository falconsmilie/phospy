from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import pandas as pd

from ..errors import TableSchemaError
from .results import SamplingTraceOverrideEnsemble


def _iter_trace_table_frames(
    trace_dir: Path,
    table_name: str,
) -> Iterator[tuple[pd.DataFrame, str]]:
    csv_path = trace_dir / f"{table_name}.csv"
    if csv_path.exists():
        yield pd.read_csv(csv_path), csv_path.name
        return

    parquet_path = trace_dir / f"{table_name}.parquet"
    parquet_parts = sorted(trace_dir.glob(f"{table_name}.part-*.parquet"))
    parquet_sources: list[Path] = []
    if parquet_path.exists():
        parquet_sources.append(parquet_path)
    else:
        parquet_sources.extend(parquet_parts)

    if not parquet_sources:
        return

    for source in parquet_sources:
        try:
            yield pd.read_parquet(source), source.name
        except (
            ImportError,
            ModuleNotFoundError,
            ValueError,
        ) as error:  # pragma: no cover - engine availability varies
            msg = (
                "Unable to read parquet trace replay input. Install a supported "
                "parquet engine such as 'pyarrow', or provide CSV trace files."
            )
            raise RuntimeError(msg) from error


def _require_trace_columns(
    frame: pd.DataFrame,
    *,
    required_columns: tuple[str, ...],
    source_label: str,
) -> pd.DataFrame:
    required_set = set(required_columns)
    if not required_set.issubset(frame.columns):
        missing = sorted(required_set.difference(frame.columns))
        msg = f"{source_label} is missing required columns: " + ", ".join(missing)
        raise TableSchemaError(msg)
    return frame.loc[:, list(required_columns)]


class PredictionSamplingTrace:
    """Replay-ready view of previously exported sampling trace tables.

    This advanced seam is intended for parity and deterministic replay
    workflows where sampling draws must be pinned from recorded traces.
    """

    def __init__(
        self, ensembles_by_kinase: dict[str, dict[int, SamplingTraceOverrideEnsemble]]
    ) -> None:
        self.ensembles_by_kinase = ensembles_by_kinase

    @classmethod
    def from_trace_directory(cls, trace_dir: str | Path) -> PredictionSamplingTrace:
        """Load replay overrides from a trace-output directory."""

        path = Path(trace_dir)
        initial_required_cols = ("kinase", "ensemble", "draw", "site")
        sample_required_cols = (
            "kinase",
            "ensemble",
            "iteration",
            "class_label",
            "draw",
            "site",
        )

        initial_records: dict[tuple[str, int], list[tuple[int, str]]] = {}
        sample_records: dict[tuple[str, int, int, int], list[tuple[int, str]]] = {}
        found_initial = False
        found_samples = False

        for frame, label in _iter_trace_table_frames(path, "trace_initial_negatives"):
            found_initial = True
            valid_frame = _require_trace_columns(
                frame,
                required_columns=initial_required_cols,
                source_label=label,
            )
            for kinase, ensemble, draw, site in valid_frame.itertuples(
                index=False, name=None
            ):
                key = (str(kinase), int(ensemble))
                initial_records.setdefault(key, []).append((int(draw), str(site)))

        for frame, label in _iter_trace_table_frames(path, "trace_iteration_samples"):
            found_samples = True
            valid_frame = _require_trace_columns(
                frame,
                required_columns=sample_required_cols,
                source_label=label,
            )
            for (
                kinase,
                ensemble,
                iteration,
                class_label,
                draw,
                site,
            ) in valid_frame.itertuples(index=False, name=None):
                key = (str(kinase), int(ensemble), int(iteration), int(class_label))
                sample_records.setdefault(key, []).append((int(draw), str(site)))

        if not found_initial and not found_samples:
            msg = (
                "sampling trace directory must contain trace_initial_negatives "
                "and/or trace_iteration_samples in CSV or parquet format"
            )
            raise TableSchemaError(msg)

        ensembles_by_kinase: dict[str, dict[int, SamplingTraceOverrideEnsemble]] = {}

        for (kinase, ensemble), draws in sorted(initial_records.items()):
            ordered_sites = [
                site for _, site in sorted(draws, key=lambda item: item[0])
            ]
            ensemble_map = ensembles_by_kinase.setdefault(kinase, {})
            ensemble_map[ensemble] = SamplingTraceOverrideEnsemble(
                initial_negative_sites=ordered_sites,
                iteration_sample_sites={},
            )

        for (kinase, ensemble, iteration, class_label), draws in sorted(
            sample_records.items()
        ):
            ensemble_map = ensembles_by_kinase.setdefault(kinase, {})
            ensemble_override = ensemble_map.setdefault(
                ensemble,
                SamplingTraceOverrideEnsemble(
                    initial_negative_sites=None,
                    iteration_sample_sites={},
                ),
            )
            iteration_map = ensemble_override.iteration_sample_sites.setdefault(
                iteration, {}
            )
            iteration_map[class_label] = [
                site for _, site in sorted(draws, key=lambda item: item[0])
            ]

        return cls(ensembles_by_kinase=ensembles_by_kinase)

    def get_ensemble_override(
        self, kinase: str, ensemble_index: int
    ) -> SamplingTraceOverrideEnsemble | None:
        """Return replay overrides for one kinase-ensemble pair."""

        return self.ensembles_by_kinase.get(kinase, {}).get(ensemble_index)

    def subset_kinases(self, kinases: list[str] | set[str]) -> PredictionSamplingTrace:
        """Return a replay trace containing only the requested kinases."""

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

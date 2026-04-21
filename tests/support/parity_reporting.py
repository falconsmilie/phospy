from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import pytest
from _pytest.config import Config
from _pytest.terminal import TerminalReporter

PARITY_FAMILY_ORDER: tuple[str, ...] = (
    "prediction_science",
    "kinase_workflow",
    "l6_prediction",
    "adaptive_prediction",
    "adaptive_replay",
    "public_predmat",
    "order_invariance",
    "activity_stage",
    "preprocessing_science",
    "signalome_workflow",
)

PARITY_FAMILY_TITLES: Mapping[str, str] = {
    "prediction_science": "Prediction-science parity",
    "kinase_workflow": "Kinase workflow parity",
    "l6_prediction": "L6 core kinase scoring/prediction parity",
    "adaptive_prediction": "Adaptive prediction parity",
    "adaptive_replay": "Core kinase adaptive replay-trace parity",
    "public_predmat": "Public end-to-end predMat parity",
    "order_invariance": "Public predMat order-invariance parity",
    "activity_stage": "Activity-stage parity",
    "preprocessing_science": "Preprocessing-science parity",
    "signalome_workflow": "Signalome workflow parity",
}

PARITY_FILE_TO_FAMILY: Mapping[str, str] = {
    "test_prediction_science_parity.py": "prediction_science",
    "test_kinase_workflow_parity.py": "kinase_workflow",
    "test_l6_prediction_parity.py": "l6_prediction",
    "test_adaptive_prediction_parity.py": "adaptive_prediction",
    "test_adaptive_replay_parity.py": "adaptive_replay",
    "test_public_predmat_parity.py": "public_predmat",
    "test_activity_stage_parity.py": "activity_stage",
    "test_preprocessing_science_parity.py": "preprocessing_science",
    "test_signalome_workflow_parity.py": "signalome_workflow",
}


def infer_parity_family_from_path(
    path: Path, *, test_name: str | None = None
) -> str | None:
    if "parity" not in path.parts:
        return None
    if (
        path.name == "test_public_predmat_parity.py"
        and test_name is not None
        and "order_invariant" in test_name
    ):
        return "order_invariance"
    return PARITY_FILE_TO_FAMILY.get(path.name)


def format_shape(rows: int, columns: int) -> str:
    return f"{rows} x {columns}"


def format_float(value: float | int | None, *, precision: int = 6) -> str:
    if value is None:
        return "n/a"
    numeric = float(value)
    if math.isnan(numeric):
        return "n/a"
    return f"{numeric:.{precision}g}"


def format_percent(value: float | None, *, digits: int = 2) -> str:
    if value is None:
        return "n/a"
    numeric = float(value)
    if math.isnan(numeric):
        return "n/a"
    return f"{numeric * 100:.{digits}f}%"


def format_fraction(
    numerator: int,
    denominator: int,
    *,
    include_percent: bool = False,
    digits: int = 2,
) -> str:
    if denominator <= 0:
        return f"{numerator}/{denominator}"
    if not include_percent:
        return f"{numerator}/{denominator}"
    ratio = numerator / denominator
    return f"{numerator}/{denominator} ({ratio * 100:.{digits}f}%)"


def format_bool(value: bool) -> str:
    return "yes" if value else "no"


@dataclass(slots=True)
class _FamilyMetrics:
    passed: int = 0
    failed: int = 0
    skipped: int = 0
    metric_order: list[str] = field(default_factory=list)
    metrics: dict[str, str] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def record_outcome(self, outcome: str) -> None:
        if outcome == "passed":
            self.passed += 1
        elif outcome == "failed":
            self.failed += 1
        elif outcome == "skipped":
            self.skipped += 1

    def record_metrics(
        self,
        metrics: Sequence[tuple[str, object]],
        notes: Sequence[str] | None = None,
    ) -> None:
        for key, raw_value in metrics:
            if key not in self.metrics:
                self.metric_order.append(key)
            self.metrics[key] = _format_metric_value(raw_value)
        if notes:
            for note in notes:
                text = str(note).strip()
                if text and text not in self.notes:
                    self.notes.append(text)

    @property
    def total_tests(self) -> int:
        return self.passed + self.failed + self.skipped


@dataclass(slots=True)
class ParityReporter:
    families: dict[str, _FamilyMetrics] = field(
        default_factory=lambda: {
            family: _FamilyMetrics() for family in PARITY_FAMILY_ORDER
        }
    )

    def record_outcome(self, family: str, outcome: str) -> None:
        self._family(family).record_outcome(outcome)

    def record_metrics(
        self,
        family: str,
        metrics: Sequence[tuple[str, object]],
        notes: Sequence[str] | None = None,
    ) -> None:
        self._family(family).record_metrics(metrics, notes)

    def has_parity_activity(self) -> bool:
        for family in self.families.values():
            if family.total_tests > 0 or family.metrics:
                return True
        return False

    def iter_rendered_families(self) -> list[str]:
        families = [
            family
            for family in PARITY_FAMILY_ORDER
            if self.families[family].total_tests > 0 or self.families[family].metrics
        ]
        if families:
            return families
        return []

    def _family(self, family: str) -> _FamilyMetrics:
        if family not in self.families:
            self.families[family] = _FamilyMetrics()
        return self.families[family]


def _format_metric_value(value: object) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, bool):
        return format_bool(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return format_float(value)
    if value is None:
        return "n/a"
    return str(value)


PARITY_REPORTER_KEY: pytest.StashKey[ParityReporter] = pytest.StashKey()


def get_or_create_parity_reporter(config: Config) -> ParityReporter:
    reporter = config.stash.get(PARITY_REPORTER_KEY, None)
    if reporter is None:
        reporter = ParityReporter()
        config.stash[PARITY_REPORTER_KEY] = reporter
    return reporter


def record_parity_metrics(
    config: Config,
    *,
    family: str,
    metrics: Mapping[str, object] | Sequence[tuple[str, object]],
    notes: Sequence[str] | None = None,
) -> None:
    items = list(metrics.items()) if isinstance(metrics, Mapping) else list(metrics)
    get_or_create_parity_reporter(config).record_metrics(
        family=family,
        metrics=items,
        notes=notes,
    )


def render_parity_terminal_summary(
    terminal_reporter: TerminalReporter,
    reporter: ParityReporter,
) -> None:
    families = reporter.iter_rendered_families()
    if not families:
        return
    terminal_reporter.section("Rewrite parity summary", sep="=")
    for family in families:
        family_metrics = reporter.families[family]
        title = PARITY_FAMILY_TITLES.get(family, family.replace("_", " ").title())
        terminal_reporter.write_line(f"{title}:")
        terminal_reporter.write_line(
            "  tests: "
            f"passed={family_metrics.passed} "
            f"failed={family_metrics.failed} "
            f"skipped={family_metrics.skipped}"
        )
        if family_metrics.metric_order:
            for key in family_metrics.metric_order:
                terminal_reporter.write_line(f"  {key}: {family_metrics.metrics[key]}")
        else:
            terminal_reporter.write_line(
                "  metrics: not reported (family selected without metric capture)"
            )
        for note in family_metrics.notes:
            terminal_reporter.write_line(f"  note: {note}")

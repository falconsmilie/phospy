from __future__ import annotations

import sys
from pathlib import Path

import pytest

from tests.support.parity_reporting import (
    get_or_create_parity_reporter,
    infer_parity_family_from_path,
    render_parity_terminal_summary,
)

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def pytest_configure(config: pytest.Config) -> None:
    get_or_create_parity_reporter(config)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item: pytest.Item, call: pytest.CallInfo[None]):
    outcome = yield
    report = outcome.get_result()
    if "parity" not in item.keywords:
        return

    family = infer_parity_family_from_path(
        Path(str(item.path)),
        test_name=item.name,
    )
    if family is None:
        return

    status: str | None = None
    if report.when == "call":
        if report.passed:
            status = "passed"
        elif report.failed:
            status = "failed"
        elif report.skipped:
            status = "skipped"
    elif report.when == "setup" and report.skipped:
        status = "skipped"

    if status is not None:
        get_or_create_parity_reporter(item.config).record_outcome(family, status)


def pytest_terminal_summary(
    terminalreporter,
    exitstatus: int,
    config: pytest.Config,
) -> None:
    del exitstatus
    reporter = get_or_create_parity_reporter(config)
    if reporter.has_parity_activity():
        render_parity_terminal_summary(terminalreporter, reporter)

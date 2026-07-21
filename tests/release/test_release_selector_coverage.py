from __future__ import annotations

from pathlib import Path

import pytest

from tools.testing.release_selector_coverage import (
    audit_authoritative_release_coverage,
)

ROOT = Path(__file__).resolve().parents[2]

pytestmark = pytest.mark.release_gate


def test_authoritative_release_targets_select_every_release_blocking_node() -> None:
    audit = audit_authoritative_release_coverage(ROOT)

    assert audit.missing_nodes == {}, audit.format_failure()

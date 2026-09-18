from __future__ import annotations

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
RELEASE_NOTES = ROOT / "docs" / "release-notes.md"
POLICY_SECTION_HEADING = "## Kinase Scientific-Policy Versions"
POLICY_TABLE_HEADER = "| Policy | Implemented version |"

pytestmark = pytest.mark.release_gate


def test_current_release_notes_exclude_preexisting_kinase_policy_inventory() -> None:
    release_notes = _read_release_notes()
    assert POLICY_SECTION_HEADING not in release_notes
    assert POLICY_TABLE_HEADER not in release_notes
    assert "KSEA activity policy" not in release_notes
    assert "Membership-selection policy" not in release_notes
    assert "Inferential policy" not in release_notes
    assert "Membership payload schema" not in release_notes
    assert "Membership-independence policy" not in release_notes


def _read_release_notes() -> str:
    return RELEASE_NOTES.read_text(encoding="utf-8").replace("\r\n", "\n")

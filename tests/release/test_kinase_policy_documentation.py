from __future__ import annotations

import re
from pathlib import Path

import pytest

from phospy.science.activities.membership import (
    ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION,
    ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
    KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_VERSION,
    KSEA_MEMBERSHIP_INFERENTIAL_POLICY_VERSION,
)
from phospy.science.activities.scientific_policies import (
    KSEA_ZSCORE_ACTIVITY_POLICY_VERSION,
)

ROOT = Path(__file__).resolve().parents[2]
RELEASE_NOTES = ROOT / "docs" / "release-notes.md"
POLICY_SECTION_HEADING = "## Kinase Scientific-Policy Versions"
POLICY_TABLE_HEADER = "| Policy | Implemented version |"

pytestmark = pytest.mark.release_gate


def test_current_release_notes_kinase_policy_versions_match_source_constants() -> None:
    release_notes = _read_release_notes()
    assert release_notes.count(POLICY_TABLE_HEADER) == 1
    section = _single_current_policy_section(release_notes)
    table = _single_policy_table(section)

    documented_versions = _parse_policy_table(table)

    assert documented_versions == _source_policy_versions()


def test_current_release_notes_do_not_duplicate_kinase_policy_version_claims() -> None:
    release_notes = _read_release_notes()
    section = _single_current_policy_section(release_notes)
    table = _single_policy_table(section)
    release_notes_without_table = release_notes.replace(table, "", 1)

    duplicate_claims = _policy_version_claims(release_notes_without_table)

    assert duplicate_claims == []


def _read_release_notes() -> str:
    return RELEASE_NOTES.read_text(encoding="utf-8").replace("\r\n", "\n")


def _source_policy_versions() -> dict[str, str]:
    return {
        "KSEA activity policy": KSEA_ZSCORE_ACTIVITY_POLICY_VERSION,
        "Membership-selection policy": ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
        "Inferential policy": KSEA_MEMBERSHIP_INFERENTIAL_POLICY_VERSION,
        "Membership payload schema": (
            ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION
        ),
        "Membership-independence policy": (KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_VERSION),
    }


def _single_current_policy_section(release_notes: str) -> str:
    pattern = re.compile(
        rf"(?ms)^{re.escape(POLICY_SECTION_HEADING)}\n(?P<section>.*?)(?=^## |\Z)"
    )
    matches = tuple(pattern.finditer(release_notes))
    assert len(matches) == 1
    return matches[0].group("section")


def _single_policy_table(section: str) -> str:
    lines = section.splitlines()
    header_indexes = [
        index for index, line in enumerate(lines) if line.strip() == POLICY_TABLE_HEADER
    ]
    assert len(header_indexes) == 1

    start = header_indexes[0]
    end = start
    while end < len(lines) and lines[end].strip().startswith("|"):
        end += 1
    return "\n".join(lines[start:end])


def _parse_policy_table(table: str) -> dict[str, str]:
    documented: dict[str, str] = {}
    rows = table.splitlines()
    for line in rows[2:]:
        cells = [cell.strip().strip("`") for cell in line.strip().strip("|").split("|")]
        assert len(cells) == 2
        policy, version = cells
        assert policy not in documented
        assert re.fullmatch(r"\d+", version) is not None
        documented[policy] = version
    return documented


def _policy_version_claims(text: str) -> list[str]:
    fragments = (
        r"KSEA activity policy",
        r"membership-selection policy",
        r"(?:KSEA membership )?inferential policy",
        r"membership(?:-selection)? (?:payload )?schema",
        r"(?:membership-independence policy|independence-evidence|"
        r"KSEA membership independence-policy evidence)",
    )
    claims: list[str] = []
    for fragment in fragments:
        pattern = re.compile(
            rf"{fragment}\s+(?:is\s+)?version\s+`?\d+`?",
            flags=re.IGNORECASE,
        )
        for match in pattern.finditer(text):
            claims.append(" ".join(match.group(0).split()))
    return claims

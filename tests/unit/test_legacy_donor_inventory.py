from __future__ import annotations

from pathlib import Path

from tests.support.legacy_donor_inventory import (
    CLOSED_SCIENCE_GAP_TICKETS,
    LEGACY_SCIENCE_AREAS,
    LEGACY_SCIENCE_STATUS_VALUES,
    OPEN_LEGACY_SCIENCE_AREAS,
    OPEN_SCIENCE_GAP_TICKETS,
    REQUIRED_LEGACY_SCIENCE_AREAS,
    STATUS_CONTRACT_CHANGED,
    STATUS_OPEN_GAP,
    STATUS_PORTED,
    TRACKED_SCIENCE_GAP_TICKETS,
)

ROOT = Path(__file__).resolve().parents[2]


def _assert_test_reference_exists(
    reference: str, *, allowed_prefixes: tuple[str, ...]
) -> None:
    file_path, _, test_name = reference.partition("::")
    assert test_name, f"missing test function qualifier in reference: {reference}"
    assert file_path.startswith(allowed_prefixes), (
        f"unexpected reference prefix for {reference}: "
        f"expected one of {allowed_prefixes}"
    )
    file_on_disk = ROOT / file_path
    assert file_on_disk.exists(), f"missing referenced test file: {file_path}"
    source = file_on_disk.read_text(encoding="utf-8")
    assert f"def {test_name}(" in source, (
        f"missing referenced test function: {reference}"
    )


def _parse_legacy_science_inventory_table(doc_text: str) -> dict[str, dict[str, str]]:
    lines = doc_text.splitlines()
    heading_index = next(
        (
            index
            for index, line in enumerate(lines)
            if line.startswith("## Legacy Science Coverage Inventory")
        ),
        None,
    )
    assert heading_index is not None, "missing legacy-science inventory heading"

    table_start = next(
        (
            index
            for index in range(heading_index + 1, len(lines))
            if lines[index].lstrip().startswith("|")
        ),
        None,
    )
    assert table_start is not None, "missing legacy-science inventory table"
    assert table_start + 1 < len(lines), "missing inventory table separator"

    header_cells = [
        cell.strip() for cell in lines[table_start].strip().strip("|").split("|")
    ]
    assert "Legacy science area" in header_cells
    assert "Coverage tier" in header_cells

    rows: dict[str, dict[str, str]] = {}
    for line in lines[table_start + 2 :]:
        if not line.lstrip().startswith("|"):
            break
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) != len(header_cells):
            continue
        row = dict(zip(header_cells, cells, strict=True))
        area = row.get("Legacy science area")
        if area:
            rows[area] = row
    return rows


def _inventory_by_area():
    return {entry.area: entry for entry in LEGACY_SCIENCE_AREAS}


def test_required_legacy_science_areas_are_all_in_inventory() -> None:
    inventory_areas = {entry.area for entry in LEGACY_SCIENCE_AREAS}
    assert inventory_areas == set(REQUIRED_LEGACY_SCIENCE_AREAS)


def test_inventory_statuses_use_supported_vocabulary() -> None:
    for entry in LEGACY_SCIENCE_AREAS:
        assert entry.status in LEGACY_SCIENCE_STATUS_VALUES
        assert entry.status_summary


def test_non_open_areas_have_rewrite_owned_evidence() -> None:
    for entry in LEGACY_SCIENCE_AREAS:
        rewrite_tests = (
            entry.rewrite_unit_tests
            + entry.rewrite_parity_tests
            + entry.rewrite_integration_tests
        )
        if entry.status == STATUS_OPEN_GAP:
            assert not rewrite_tests, (
                "open-gap areas should not be represented as rewrite-ported blockers: "
                f"{entry.area}"
            )
            continue
        assert rewrite_tests, f"no rewrite-owned evidence mapped for {entry.area}"
        for reference in rewrite_tests:
            _assert_test_reference_exists(reference, allowed_prefixes=("tests/",))


def test_tracked_science_gap_tickets_map_to_rewrite_tests_in_inventory() -> None:
    ticket_to_entry = {
        entry.science_gap_ticket: entry
        for entry in LEGACY_SCIENCE_AREAS
        if entry.science_gap_ticket is not None
    }
    assert set(TRACKED_SCIENCE_GAP_TICKETS) == set(ticket_to_entry)
    for ticket, entry in ticket_to_entry.items():
        rewrite_tests = (
            entry.rewrite_unit_tests
            + entry.rewrite_parity_tests
            + entry.rewrite_integration_tests
        )
        assert rewrite_tests, f"{ticket} must map to at least one rewrite-owned test"
        assert entry.status in {STATUS_PORTED, STATUS_CONTRACT_CHANGED}


def test_open_gap_truth_source_distinguishes_ticket_status_from_coverage_status() -> (
    None
):
    assert OPEN_SCIENCE_GAP_TICKETS == ()
    assert OPEN_LEGACY_SCIENCE_AREAS == ()
    assert set(OPEN_SCIENCE_GAP_TICKETS).issubset(TRACKED_SCIENCE_GAP_TICKETS)
    assert set(CLOSED_SCIENCE_GAP_TICKETS) == set(TRACKED_SCIENCE_GAP_TICKETS)


def test_open_gap_areas_remain_linked_to_archival_donor_evidence() -> None:
    for entry in LEGACY_SCIENCE_AREAS:
        if entry.status != STATUS_OPEN_GAP:
            continue
        assert entry.archival_only_tests, (
            f"missing archival references for {entry.area}"
        )
        for reference in entry.archival_only_tests:
            _assert_test_reference_exists(
                reference, allowed_prefixes=("tests_legacy/",)
            )


def test_ported_areas_pin_promoted_fixtures_with_provenance() -> None:
    for entry in LEGACY_SCIENCE_AREAS:
        if entry.status != STATUS_PORTED:
            continue
        assert entry.promoted_fixture_paths, (
            f"ported area should pin promoted fixtures: {entry.area}"
        )
        for fixture_path in entry.promoted_fixture_paths:
            assert not fixture_path.startswith("tests_legacy/")
            assert (ROOT / fixture_path).exists(), (
                f"missing promoted fixture: {fixture_path}"
            )
        assert entry.provenance_paths
        for provenance_path in entry.provenance_paths:
            assert (ROOT / provenance_path).exists(), (
                f"missing provenance: {provenance_path}"
            )


def test_key_legacy_science_areas_pin_expected_rewrite_parity_blockers() -> None:
    expected_parity_tests = {
        "adaptive sampling / svm_mode": (
            "tests/parity/test_adaptive_prediction_parity.py::"
            "test_adaptive_prediction_cross_policy_divergence_stable_vs_r_parity",
        ),
        "expanded signalome outputs": (
            "tests/parity/test_signalome_workflow_parity.py::"
            "test_signalome_expanded_signalome_matches_l6_full_fixture_table_with_tolerance",
        ),
    }
    inventory_by_area = _inventory_by_area()
    for area, required_references in expected_parity_tests.items():
        entry = inventory_by_area[area]
        assert entry.rewrite_parity_tests, f"{area} must have rewrite parity blockers"
        for reference in required_references:
            assert reference in entry.rewrite_parity_tests, (
                f"{area} missing required rewrite parity blocker: {reference}"
            )
            _assert_test_reference_exists(reference, allowed_prefixes=("tests/",))


def test_governance_docs_keep_legacy_science_inventory_statuses_in_sync() -> None:
    parity_doc = (ROOT / "docs" / "parity.md").read_text(encoding="utf-8")
    audit_doc = (
        ROOT / "docs" / "architecture" / "legacy_science_gap_audit.md"
    ).read_text(encoding="utf-8")
    assert "Legacy Science Coverage Inventory" in parity_doc
    assert "Legacy Science Coverage Inventory" in audit_doc
    for entry in LEGACY_SCIENCE_AREAS:
        row = f"| {entry.area} | {entry.status} |"
        assert row in parity_doc, f"missing parity row: {row}"
        assert row in audit_doc, f"missing audit row: {row}"


def test_governance_docs_keep_legacy_science_coverage_tiers_in_sync() -> None:
    parity_doc = (ROOT / "docs" / "parity.md").read_text(encoding="utf-8")
    audit_doc = (
        ROOT / "docs" / "architecture" / "legacy_science_gap_audit.md"
    ).read_text(encoding="utf-8")
    parity_rows = _parse_legacy_science_inventory_table(parity_doc)
    audit_rows = _parse_legacy_science_inventory_table(audit_doc)
    for entry in LEGACY_SCIENCE_AREAS:
        assert entry.area in parity_rows, f"missing parity inventory area: {entry.area}"
        assert entry.area in audit_rows, f"missing audit inventory area: {entry.area}"
        parity_tier = parity_rows[entry.area]["Coverage tier"]
        audit_tier = audit_rows[entry.area]["Coverage tier"]
        assert parity_tier == audit_tier, (
            f"coverage tier drift for {entry.area}: "
            f"parity={parity_tier!r}, audit={audit_tier!r}"
        )


def test_parity_gated_areas_require_active_rewrite_parity_tests() -> None:
    audit_doc = (
        ROOT / "docs" / "architecture" / "legacy_science_gap_audit.md"
    ).read_text(encoding="utf-8")
    audit_rows = _parse_legacy_science_inventory_table(audit_doc)
    inventory_by_area = _inventory_by_area()
    for area, row in audit_rows.items():
        if row["Coverage tier"] != "PARITY_GATED_ACTIVE_SCIENCE":
            continue
        entry = inventory_by_area[area]
        assert entry.rewrite_parity_tests, (
            f"parity-gated area must map to active rewrite parity tests: {area}"
        )
        assert "tests/parity/" in row.get("Active rewrite test evidence", ""), (
            f"parity-gated area must cite tests/parity evidence in audit table: {area}"
        )
        for reference in entry.rewrite_parity_tests:
            assert reference.startswith("tests/parity/"), (
                "parity-gated area test reference must be in tests/parity: "
                f"{area} -> {reference}"
            )
            _assert_test_reference_exists(
                reference, allowed_prefixes=("tests/parity/",)
            )


def test_release_facing_parity_docs_do_not_claim_blanket_no_open_legacy_science_gaps() -> (
    None
):
    parity_doc = (ROOT / "docs" / "parity.md").read_text(encoding="utf-8")
    roadmap_doc = (ROOT / "docs" / "roadmap.md").read_text(encoding="utf-8")
    readme_doc = (ROOT / "README.md").read_text(encoding="utf-8")
    audit_doc = (
        ROOT / "docs" / "architecture" / "legacy_science_gap_audit.md"
    ).read_text(encoding="utf-8")

    forbidden_phrases = (
        "No open science-parity gap tickets are confirmed in the supported rewrite lane",
        "None confirmed in this scoped pass.",
        "does not apply additional preprocessing transforms to quantitative matrices.",
        "The remaining roadmap is documentation and governance alignment, not core",
    )
    for phrase in forbidden_phrases:
        assert phrase not in parity_doc
        assert phrase not in roadmap_doc
        assert phrase not in readme_doc
        assert phrase not in audit_doc

    assert "| site-matrix construction | PORTED |" in parity_doc
    assert "| site-matrix construction | PORTED |" in audit_doc

from __future__ import annotations

from pathlib import Path

from tests.support.legacy_donor_inventory import (
    LEGACY_DONOR_AREAS,
    OPEN_SCIENCE_GAP_TICKETS,
    REQUIRED_DONOR_AREAS,
)

ROOT = Path(__file__).resolve().parents[2]


def _assert_test_reference_exists(reference: str) -> None:
    file_path, _, test_name = reference.partition("::")
    assert file_path.startswith("tests/")
    file_on_disk = ROOT / file_path
    assert file_on_disk.exists(), f"missing referenced test file: {file_path}"
    source = file_on_disk.read_text(encoding="utf-8")
    assert f"def {test_name}(" in source, (
        f"missing referenced test function: {reference}"
    )


def _inventory_by_area():
    return {entry.area: entry for entry in LEGACY_DONOR_AREAS}


def test_required_legacy_donor_areas_are_all_in_inventory() -> None:
    inventory_areas = {entry.area for entry in LEGACY_DONOR_AREAS}
    assert inventory_areas == set(REQUIRED_DONOR_AREAS)


def test_each_donor_area_has_rewrite_owned_blocking_tests() -> None:
    for entry in LEGACY_DONOR_AREAS:
        rewrite_tests = (
            entry.rewrite_unit_tests
            + entry.rewrite_parity_tests
            + entry.rewrite_integration_tests
        )
        assert rewrite_tests, f"no rewrite-owned blocker mapped for {entry.area}"
        for reference in rewrite_tests:
            _assert_test_reference_exists(reference)


def test_open_science_gap_tickets_map_to_rewrite_tests_in_inventory() -> None:
    ticket_to_entry = {entry.science_gap_ticket: entry for entry in LEGACY_DONOR_AREAS}
    assert set(OPEN_SCIENCE_GAP_TICKETS) == set(ticket_to_entry)
    for ticket, entry in ticket_to_entry.items():
        rewrite_tests = (
            entry.rewrite_unit_tests
            + entry.rewrite_parity_tests
            + entry.rewrite_integration_tests
        )
        assert rewrite_tests, f"{ticket} must map to at least one rewrite-owned test"


def test_inventory_fixtures_are_promoted_under_rewrite_paths_with_provenance() -> None:
    for entry in LEGACY_DONOR_AREAS:
        assert entry.promoted_fixture_paths
        for fixture_path in entry.promoted_fixture_paths:
            assert not fixture_path.startswith("tests_legacy/")
            assert (ROOT / fixture_path).exists(), (
                f"missing promoted fixture: {fixture_path}"
            )
        for provenance_path in entry.provenance_paths:
            assert (ROOT / provenance_path).exists(), (
                f"missing provenance: {provenance_path}"
            )


def test_key_donor_areas_pin_expected_rewrite_parity_blockers() -> None:
    expected_parity_tests = {
        "adaptive sampling / svm_mode": (
            "tests/parity/test_adaptive_prediction_parity.py::"
            "test_adaptive_ensemble_outputs_match_promoted_fixture_tolerances",
        ),
        "expanded signalome outputs": (
            "tests/parity/test_signalome_workflow_parity.py::"
            "test_signalome_expanded_slice_matches_l6_selected_akt1_fixture",
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
            _assert_test_reference_exists(reference)


def test_parity_doc_keeps_donor_inventory_visible() -> None:
    parity_doc = (ROOT / "docs" / "parity.md").read_text(encoding="utf-8")
    assert "Legacy Donor Promotion Inventory" in parity_doc
    assert "adaptive prediction parity" in parity_doc
    assert "`expanded_signalome`" in parity_doc
    for area in REQUIRED_DONOR_AREAS:
        assert area in parity_doc
    for ticket in OPEN_SCIENCE_GAP_TICKETS:
        assert ticket in parity_doc

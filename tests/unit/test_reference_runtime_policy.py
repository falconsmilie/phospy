from __future__ import annotations

import pytest

from phospy.errors import UnsupportedOrganismError
from phospy.science.references.models import Organism, ReferencePreset
from phospy.science.references.resolution import ReferenceResolver
from phospy.science.references.resources import (
    available_bundled_reference_lanes,
    supported_bundled_organisms,
)


def test_bundled_runtime_support_includes_existing_rat_lane() -> None:
    assert Organism.RAT in supported_bundled_organisms()


def test_bundled_runtime_lane_inventory_reports_existing_rat_lane() -> None:
    lanes = available_bundled_reference_lanes()

    rat_lanes = [lane for lane in lanes if lane.organism is Organism.RAT]
    assert len(rat_lanes) == 1
    lane = rat_lanes[0]
    assert lane.bundle_id == "l6_native"
    assert lane.source_version == "PhosR 1.20.0"
    assert lane.retrieved_at.isoformat() == "2026-04-16"
    assert lane.redistribution_status
    assert lane.supports
    assert lane.limitations


def test_rat_preset_resolves_to_packaged_bundled_lane() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )
    assert bundle.organism is Organism.RAT
    assert not bundle.kinase_substrate_map.empty
    assert not bundle.site_sequences.empty
    assert len(bundle.kinase_substrate_map) == 117
    assert len(bundle.site_sequences) == 589


@pytest.mark.parametrize(
    ("preset", "dataset_organism"),
    [
        (ReferencePreset.HUMAN, Organism.HUMAN),
        (ReferencePreset.MOUSE, Organism.MOUSE),
    ],
)
def test_human_mouse_explicit_presets_match_committed_bundle_presence(
    preset: ReferencePreset,
    dataset_organism: Organism,
) -> None:
    if dataset_organism in supported_bundled_organisms():
        bundle = ReferenceResolver().run(
            preset,
            dataset_organism=dataset_organism,
        )
        assert bundle.organism is dataset_organism
        return

    with pytest.raises(UnsupportedOrganismError) as exc_info:
        ReferenceResolver().run(
            preset,
            dataset_organism=dataset_organism,
        )
    message = str(exc_info.value)
    assert "supported bundled organisms:" in message
    assert "ReferenceBundle(organism=..." in message
    assert (
        "reference-bundle docs: https://phospy.com/docs/api/guide/#references"
        in message
    )


@pytest.mark.parametrize("dataset_organism", [Organism.HUMAN, Organism.MOUSE])
def test_human_mouse_auto_resolution_matches_committed_bundle_presence(
    dataset_organism: Organism,
) -> None:
    if dataset_organism in supported_bundled_organisms():
        bundle = ReferenceResolver().run(
            ReferencePreset.AUTO,
            dataset_organism=dataset_organism,
        )
        assert bundle.organism is dataset_organism
        return

    with pytest.raises(UnsupportedOrganismError) as exc_info:
        ReferenceResolver().run(
            ReferencePreset.AUTO,
            dataset_organism=dataset_organism,
        )
    message = str(exc_info.value)
    assert "supported bundled organisms:" in message
    assert "ReferenceBundle(organism=..." in message
    assert (
        "reference-bundle docs: https://phospy.com/docs/api/guide/#references"
        in message
    )

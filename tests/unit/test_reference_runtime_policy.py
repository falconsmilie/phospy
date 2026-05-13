from __future__ import annotations

import pytest

from phospy.errors import UnsupportedOrganismError
from phospy.science.references.models import Organism, ReferencePreset
from phospy.science.references.resolution import ReferenceResolver
from phospy.science.references.resources import supported_bundled_organisms


def test_bundled_runtime_support_is_rat_only() -> None:
    assert supported_bundled_organisms() == (Organism.RAT,)


def test_rat_preset_resolves_to_packaged_bundled_lane() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )
    assert bundle.organism is Organism.RAT
    assert not bundle.kinase_substrate_map.empty
    assert not bundle.site_sequences.empty


@pytest.mark.parametrize(
    ("preset", "dataset_organism"),
    [
        (ReferencePreset.HUMAN, Organism.HUMAN),
        (ReferencePreset.MOUSE, Organism.MOUSE),
    ],
)
def test_non_rat_explicit_presets_fail_for_bundled_resolution(
    preset: ReferencePreset,
    dataset_organism: Organism,
) -> None:
    with pytest.raises(UnsupportedOrganismError) as exc_info:
        ReferenceResolver().run(
            preset,
            dataset_organism=dataset_organism,
        )
    message = str(exc_info.value)
    assert "supported bundled organisms: rat" in message
    assert "caller-supplied ReferenceBundle" in message


@pytest.mark.parametrize("dataset_organism", [Organism.HUMAN, Organism.MOUSE])
def test_non_rat_auto_resolution_fails_for_bundled_resolution(
    dataset_organism: Organism,
) -> None:
    with pytest.raises(UnsupportedOrganismError) as exc_info:
        ReferenceResolver().run(
            ReferencePreset.AUTO,
            dataset_organism=dataset_organism,
        )
    message = str(exc_info.value)
    assert "supported bundled organisms: rat" in message
    assert "caller-supplied ReferenceBundle" in message

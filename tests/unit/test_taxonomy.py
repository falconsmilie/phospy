from __future__ import annotations

import pytest

from phospy.api import Organism as PublicOrganism
from phospy.taxonomy import Organism, normalize_optional_organism, normalize_organism


@pytest.mark.parametrize(
    ("value", "expected"),
    (
        ("human", Organism.HUMAN),
        ("Homo sapiens", Organism.HUMAN),
        ("h_sapiens", Organism.HUMAN),
        ("hsapiens", Organism.HUMAN),
        ("hsa", Organism.HUMAN),
        ("9606", Organism.HUMAN),
        ("mouse", Organism.MOUSE),
        ("Mus-musculus", Organism.MOUSE),
        ("m musculus", Organism.MOUSE),
        ("mmusculus", Organism.MOUSE),
        ("mmu", Organism.MOUSE),
        ("10090", Organism.MOUSE),
        ("rat", Organism.RAT),
        ("Rattus norvegicus", Organism.RAT),
        ("r norvegicus", Organism.RAT),
        ("rnorvegicus", Organism.RAT),
        ("rno", Organism.RAT),
        ("10116", Organism.RAT),
    ),
)
def test_organism_normalisation_preserves_supported_spellings(
    value: str,
    expected: Organism,
) -> None:
    assert (
        normalize_organism(
            value,
            field_name="request.organism",
            error_type=ValueError,
        )
        is expected
    )


def test_organism_public_identity_and_values_are_stable() -> None:
    assert PublicOrganism is Organism
    assert tuple(member.value for member in Organism) == ("human", "mouse", "rat")


def test_organism_normalisation_preserves_supported_input_errors() -> None:
    with pytest.raises(
        ValueError,
        match=(
            "request\\.organism has unsupported organism 'yeast'; "
            "supported organisms: human, mouse, rat"
        ),
    ):
        normalize_organism(
            "yeast",
            field_name="request.organism",
            error_type=ValueError,
        )

    with pytest.raises(
        TypeError,
        match="request\\.organism must be a supported organism string",
    ):
        normalize_organism(
            9606,
            field_name="request.organism",
            error_type=TypeError,
        )

    assert (
        normalize_optional_organism(
            " ",
            field_name="request.organism",
            error_type=ValueError,
        )
        is None
    )

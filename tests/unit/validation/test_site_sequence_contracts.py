from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import WorkflowValidationError
from phospy.science.sites.identity_contracts import (
    PhosphositeIdentityContract as SciencePhosphositeIdentityContract,
)
from phospy.science.sites.identity_contracts import (
    SequenceContextContract as ScienceIdentityRouteSequenceContextContract,
)
from phospy.science.sites.sequence_context import (
    WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT as SCIENCE_WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT,
)
from phospy.science.sites.sequence_context import (
    SequenceContextContract as ScienceSequenceContextContract,
)
from phospy.science.sites.sequence_context import (
    enforce_site_sequence_context_contract as science_enforce_site_sequence_context_contract,
)
from phospy.validation.datasets.site_metadata import (
    SequenceContextContract,
    enforce_site_sequence_context_contract,
)
from phospy.validation.identity_contracts import (
    WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT as VALIDATION_WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT,
)
from phospy.validation.identity_contracts import (
    PhosphositeIdentityContract as ValidationPhosphositeIdentityContract,
)
from phospy.validation.identity_contracts import (
    SequenceContextContract as ValidationIdentityRouteSequenceContextContract,
)


def _contract(
    *,
    allow_lowercase: bool = False,
    allow_terminal_padding: bool = False,
    requires_known_sequence_source: bool = False,
    required_window_length: object = 15,
    center_index: object = 7,
    allowed_residues: object = frozenset("ACDEFGHIKLMNPQRSTVWY"),
    required_center_residues: object = frozenset({"S", "T", "Y"}),
) -> SequenceContextContract:
    return SequenceContextContract(
        requires_site_sequence=True,
        requires_centered_site=True,
        required_window_length=required_window_length,  # type: ignore[arg-type]
        center_index=center_index,  # type: ignore[arg-type]
        allowed_residues=allowed_residues,  # type: ignore[arg-type]
        allow_terminal_padding=allow_terminal_padding,
        allow_lowercase=allow_lowercase,
        allow_modified_residue_symbols=False,
        required_center_residues=required_center_residues,  # type: ignore[arg-type]
        requires_known_sequence_source=requires_known_sequence_source,
        contract_id="test_fixed_15aa_window",
    )


def _site_metadata(sequence: object = "AAAAAAAYAAAAAAA") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site": ["Y182"],
            "site_sequence": [sequence],
        },
        index=pd.Index(["site-key-1"], name="site_key"),
    )


def _enforce(
    site_metadata: pd.DataFrame,
    *,
    contract: SequenceContextContract | None = None,
    sequence_source_by_site: dict[str, str] | None = None,
) -> None:
    enforce_site_sequence_context_contract(
        site_metadata=site_metadata,
        field_name="test.site_metadata",
        workflow_name="kinase workflow request",
        scoring_mode="kinase_library_contextual_motif",
        contract=_contract() if contract is None else contract,
        error_type=WorkflowValidationError,
        sequence_source_by_site=sequence_source_by_site,
    )


def test_sequence_context_contract_import_routes_are_identity_preserving() -> None:
    assert SequenceContextContract is ScienceSequenceContextContract
    assert (
        ValidationIdentityRouteSequenceContextContract is ScienceSequenceContextContract
    )
    assert ScienceIdentityRouteSequenceContextContract is ScienceSequenceContextContract
    assert (
        enforce_site_sequence_context_contract
        is science_enforce_site_sequence_context_contract
    )
    assert (
        VALIDATION_WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT
        is SCIENCE_WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT
    )


def test_phosphosite_identity_contract_import_routes_are_identity_preserving() -> None:
    assert ValidationPhosphositeIdentityContract is SciencePhosphositeIdentityContract


def test_fixed_centered_phosphosite_window_passes() -> None:
    _enforce(_site_metadata())


@pytest.mark.parametrize(
    ("overrides", "expected"),
    [
        ({"required_window_length": 0}, "required_window_length must be > 0"),
        ({"required_window_length": -1}, "required_window_length must be > 0"),
        (
            {"required_window_length": 15.0},
            "required_window_length must be an int or None",
        ),
        ({"center_index": -1}, "center_index must be >= 0"),
        ({"center_index": 15}, "center_index must be within required_window_length"),
        ({"center_index": None}, "center_index is required"),
        (
            {"allowed_residues": frozenset({"ST"})},
            "allowed_residues must contain one-character tokens",
        ),
        (
            {"required_center_residues": frozenset({1})},
            "required_center_residues must contain string residue tokens",
        ),
    ],
)
def test_sequence_context_contract_rejects_invalid_policy_values(
    overrides: dict[str, object],
    expected: str,
) -> None:
    with pytest.raises(ValueError, match=expected):
        _contract(**overrides)


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [
        (pd.NA, "missing or blank"),
        ("", "missing or blank"),
        ("AAAAAAAYAAAAAA", "expected_length=15"),
        ("AAAAAAAAAAAAAAA", "expected_center_residues=S/T/Y"),
        ("AAAAAAA*AAAAAAA", "unsupported_characters='*'"),
        ("aaaaaaayaaaaaaa", "lowercase characters are not allowed"),
        ("__AAAAAYAAAAAAA", "terminal padding is not allowed"),
    ],
)
def test_fixed_window_contract_rejects_invalid_sequence_context(
    sequence: object,
    expected: str,
) -> None:
    with pytest.raises(WorkflowValidationError) as exc_info:
        _enforce(_site_metadata(sequence))

    message = str(exc_info.value)
    assert "site-key-1" in message
    assert "kinase_library_contextual_motif" in message
    assert "expected_length=15" in message
    assert "expected_center_index=7" in message
    assert expected in message


def test_lowercase_behavior_follows_contract() -> None:
    _enforce(
        _site_metadata("aaaaaaayaaaaaaa"),
        contract=_contract(allow_lowercase=True),
    )


def test_terminal_padding_behavior_follows_contract() -> None:
    _enforce(
        _site_metadata("__AAAAAYAAAAAAA"),
        contract=_contract(allow_terminal_padding=True),
    )


def test_site_token_centre_disagreement_is_rejected() -> None:
    site_metadata = _site_metadata("AAAAAAAYAAAAAAA")
    site_metadata.loc["site-key-1", "site"] = "S182"

    with pytest.raises(WorkflowValidationError) as exc_info:
        _enforce(site_metadata)

    message = str(exc_info.value)
    assert "expected='S'" in message
    assert "observed='Y'" in message


def test_known_sequence_source_is_required_when_contract_requires_it() -> None:
    contract = _contract(requires_known_sequence_source=True)

    with pytest.raises(WorkflowValidationError) as exc_info:
        _enforce(_site_metadata(), contract=contract)

    message = str(exc_info.value)
    assert "sequence source is unknown" in message
    assert "site-key-1" in message

    _enforce(
        _site_metadata(),
        contract=contract,
        sequence_source_by_site={"site-key-1": "reference"},
    )


def test_science_owned_enforcement_uses_sequence_source_arguments() -> None:
    contract = _contract(requires_known_sequence_source=True)

    science_enforce_site_sequence_context_contract(
        site_metadata=_site_metadata(),
        field_name="test.site_metadata",
        workflow_name="kinase workflow request",
        scoring_mode="kinase_library_contextual_motif",
        contract=contract,
        error_type=WorkflowValidationError,
        sequence_source_by_site={"site-key-1": "reference"},
    )

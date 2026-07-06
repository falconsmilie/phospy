from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import WorkflowValidationError
from phospy.validation.datasets.site_metadata import (
    SequenceContextContract,
    enforce_site_sequence_context_contract,
)


def _contract(
    *,
    allow_lowercase: bool = False,
    allow_terminal_padding: bool = False,
    requires_known_sequence_source: bool = False,
) -> SequenceContextContract:
    return SequenceContextContract(
        requires_site_sequence=True,
        requires_centered_site=True,
        required_window_length=15,
        center_index=7,
        allowed_residues=frozenset("ACDEFGHIKLMNPQRSTVWY"),
        allow_terminal_padding=allow_terminal_padding,
        allow_lowercase=allow_lowercase,
        allow_modified_residue_symbols=False,
        required_center_residues=frozenset({"S", "T", "Y"}),
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


def test_fixed_centered_phosphosite_window_passes() -> None:
    _enforce(_site_metadata())


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

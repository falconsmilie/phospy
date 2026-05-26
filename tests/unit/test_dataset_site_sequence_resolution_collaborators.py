from __future__ import annotations

import pandas as pd
import pytest

from phospy.science.datasets.preprocessing.policy_models import (
    SiteSequenceConflictPolicy,
)
from phospy.science.datasets.preprocessing.site_sequence import (
    SiteSequenceConflictResolver,
    SiteSequenceMetadataUpdater,
)


@pytest.mark.parametrize(
    ("policy", "expected_action", "expected_sequence", "expected_error"),
    [
        (
            SiteSequenceConflictPolicy.PRESERVE_EXISTING,
            "preserve_existing",
            "XXXXX",
            False,
        ),
        (
            SiteSequenceConflictPolicy.REPLACE_EXISTING,
            "replace_existing",
            "AASAA",
            False,
        ),
        (SiteSequenceConflictPolicy.ERROR, "error", "XXXXX", True),
    ],
)
def test_conflict_resolver_policies_are_independent(
    policy: SiteSequenceConflictPolicy,
    expected_action: str,
    expected_sequence: str,
    expected_error: bool,
) -> None:
    outcome = SiteSequenceConflictResolver().resolve(
        existing_sequence="XXXXX",
        fasta_sequence="AASAA",
        conflict_policy=policy,
    )

    assert outcome.status == "existing_sequence_conflict"
    assert outcome.action == expected_action
    assert outcome.resolved_site_sequence == expected_sequence
    assert outcome.is_conflict is True
    assert outcome.is_error is expected_error


def test_conflict_resolver_matching_sequences_validate_existing() -> None:
    outcome = SiteSequenceConflictResolver().resolve(
        existing_sequence="AASAA",
        fasta_sequence="AASAA",
        conflict_policy=SiteSequenceConflictPolicy.ERROR,
    )

    assert outcome.status == "resolved"
    assert outcome.action == "validate_existing"
    assert outcome.resolved_site_sequence == "AASAA"
    assert outcome.is_conflict is False
    assert outcome.is_error is False


def test_metadata_updater_isolated_and_does_not_mutate_input() -> None:
    site_metadata = pd.DataFrame(
        {
            "site_sequence": [pd.NA, "  CCTCC  "],
            "protein_accession": ["P1", "P2"],
        },
        index=pd.Index(["A;S5;", "B;T6;"], name="site_id"),
    )
    original = site_metadata.copy(deep=True)
    updater = SiteSequenceMetadataUpdater(
        site_metadata=site_metadata,
        existing_site_sequence=site_metadata.loc[:, "site_sequence"],
    )

    updater.assign(row_id="A;S5;", site_sequence="AASAA")
    updated = updater.build()

    assert pd.isna(site_metadata.loc["A;S5;", "site_sequence"])
    assert site_metadata.equals(original)
    assert updated.loc["A;S5;", "site_sequence"] == "AASAA"
    assert updated.loc["B;T6;", "site_sequence"] == "  CCTCC  "

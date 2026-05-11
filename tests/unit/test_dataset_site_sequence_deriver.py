from __future__ import annotations

import pandas as pd
import pytest

from phospy.datasets.builders.sequence_derivation import SiteSequenceDeriver
from phospy.errors.input import UnsupportedInputFormatError
from phospy.errors.references import UnsupportedOrganismError
from phospy.references.models import Organism


def _site_metadata_without_sequences() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )


def test_site_sequence_deriver_preserves_provided_sequences() -> None:
    site_metadata = _site_metadata_without_sequences()
    site_metadata.loc[:, "site_sequence"] = ["SEQ_A", "SEQ_B"]
    deriver = SiteSequenceDeriver()

    resolved = deriver.run(
        site_metadata.copy(deep=True),
        organism=None,
        allow_partial=False,
        derive_missing_from_reference=False,
    )

    assert resolved.loc[:, "site_sequence"].tolist() == ["SEQ_A", "SEQ_B"]
    report = deriver.last_report
    assert report is not None
    assert report.provided_sequence_count == 2
    assert report.derived_sequence_count == 0
    assert report.unresolved_sequence_count == 0


def test_site_sequence_deriver_derives_missing_sequences_from_supported_reference() -> (
    None
):
    deriver = SiteSequenceDeriver()
    resolved = deriver.run(
        _site_metadata_without_sequences().copy(deep=True),
        organism=Organism.RAT,
        allow_partial=False,
        derive_missing_from_reference=True,
    )

    assert "site_sequence" in resolved.columns
    assert resolved.loc[:, "site_sequence"].isna().sum() == 0
    report = deriver.last_report
    assert report is not None
    assert report.provided_sequence_count == 0
    assert report.derived_sequence_count == 2
    assert report.unresolved_sequence_count == 0
    assert report.reference_support == "available"
    assert report.reference_source is not None
    assert report.reference_bundle_id == "l6_native"
    assert report.reference_manifest is not None
    assert report.reference_manifest["bundle_id"] == "l6_native"


def test_site_sequence_deriver_fails_clearly_for_blank_provided_sequence() -> None:
    site_metadata = _site_metadata_without_sequences()
    site_metadata.loc[:, "site_sequence"] = ["", "SEQ_B"]
    deriver = SiteSequenceDeriver()

    with pytest.raises(
        UnsupportedInputFormatError,
        match="invalid sites: MAPK14;Y182;",
    ):
        deriver.run(
            site_metadata.copy(deep=True),
            organism=Organism.RAT,
            allow_partial=False,
            derive_missing_from_reference=True,
        )


def test_site_sequence_deriver_fails_for_unsupported_organism_when_derivation_is_required() -> (
    None
):
    deriver = SiteSequenceDeriver()
    with pytest.raises(
        UnsupportedOrganismError, match="supported bundled organisms: rat"
    ):
        deriver.run(
            _site_metadata_without_sequences().copy(deep=True),
            organism=Organism.HUMAN,
            allow_partial=False,
            derive_missing_from_reference=True,
        )

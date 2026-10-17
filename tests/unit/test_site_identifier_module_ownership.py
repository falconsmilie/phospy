from __future__ import annotations

import importlib

import pytest

from phospy.phosphosites.identifiers import (
    ParsedSiteToken,
    SiteIdentifierNormalisationRecord,
    SiteIdentifierNormalisationReport,
    canonicalize_site_identifier,
    parse_canonical_site_identifier,
    try_parse_site_token,
)


def test_site_identifier_helpers_are_owned_by_sites_identifiers_module() -> None:
    assert canonicalize_site_identifier.__module__ == "phospy.phosphosites.identifiers"
    assert (
        parse_canonical_site_identifier.__module__ == "phospy.phosphosites.identifiers"
    )
    assert try_parse_site_token.__module__ == "phospy.phosphosites.identifiers"


def test_site_identifier_dataclasses_are_owned_by_sites_identifiers_module() -> None:
    assert ParsedSiteToken.__module__ == "phospy.phosphosites.identifiers"
    assert (
        SiteIdentifierNormalisationRecord.__module__
        == "phospy.phosphosites.identifiers"
    )
    assert (
        SiteIdentifierNormalisationReport.__module__
        == "phospy.phosphosites.identifiers"
    )


def test_root_site_ids_module_is_removed() -> None:
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("phospy.site_ids")

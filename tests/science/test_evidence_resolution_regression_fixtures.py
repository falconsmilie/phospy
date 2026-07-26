from __future__ import annotations

import hashlib
import json
from itertools import permutations
from pathlib import Path
from typing import Any

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api import PhosPyInputError
from phospy.api.requests import DATASET_MULTI_SITE_POLICY_SPLIT
from phospy.science.evidence import PeptideEvidenceDatasetResolver, PeptideEvidenceTable

pytestmark = pytest.mark.release_gate

FIXTURE_DIR = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "release_validation_regression"
    / "evidence_resolution"
)


def _manifest() -> dict[str, Any]:
    return json.loads((FIXTURE_DIR / "MANIFEST.json").read_text(encoding="utf-8"))


def _read_evidence(name: str = "peptide_evidence.csv") -> pd.DataFrame:
    frame = pd.read_csv(FIXTURE_DIR / name)
    frame.loc[:, "multi_site"] = frame.loc[:, "multi_site"].map(
        lambda value: str(value).strip().lower() == "true"
    )
    return frame


def _read_mapping() -> pd.DataFrame:
    frame = pd.read_csv(FIXTURE_DIR / "site_mapping.csv")
    frame.loc[:, "mapping_uncertainty"] = frame.loc[:, "mapping_uncertainty"].map(
        lambda value: str(value).strip().lower() == "true"
    )
    return frame


def _run_resolution(
    evidence: pd.DataFrame,
    *,
    site_mapping: pd.DataFrame | None = None,
) -> object:
    return PeptideEvidenceDatasetResolver().run(
        evidence=PeptideEvidenceTable(
            frame=evidence,
            sample_intensity_columns=("sample_a", "sample_b"),
            site_mapping=site_mapping,
        ),
        multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
    )


def test_evidence_resolution_fixture_manifest_hashes_match_files() -> None:
    manifest = _manifest()

    assert manifest["classification"] == "regression"
    assert manifest["fixture_family"] == "evidence_resolution"
    assert "not external parity" in manifest["source_policy"]
    assert manifest["seed"] == 20260724

    for file_entry in manifest["files"]:
        path = FIXTURE_DIR / str(file_entry["relative_path"])
        assert path.is_file()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == file_entry["sha256"]


def test_checked_in_evidence_fixture_resolves_equal_unequal_and_many_to_one_cases() -> (
    None
):
    resolved = _run_resolution(_read_evidence(), site_mapping=_read_mapping())
    observed = resolved.phospho.sort_index()
    expected = pd.read_csv(FIXTURE_DIR / "expected_split_phospho.csv").set_index(
        "site_id"
    )
    expected = expected.reindex(observed.index)

    pdt.assert_frame_equal(
        observed.loc[:, ["sample_a", "sample_b"]],
        expected,
        check_exact=False,
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    assert resolved.site_metadata.loc["MAPK1;S10;", "site_sequence"] == "AAASAAA"
    assert pd.isna(resolved.site_metadata.loc["MAPK1;T12;", "site_sequence"])


def test_evidence_resolution_fixture_is_row_order_invariant() -> None:
    evidence = _read_evidence()
    expected = _run_resolution(
        evidence, site_mapping=_read_mapping()
    ).phospho.sort_index()

    for row_order in permutations(range(evidence.shape[0])):
        permuted = evidence.iloc[list(row_order), :].reset_index(drop=True)
        observed = _run_resolution(
            permuted, site_mapping=_read_mapping()
        ).phospho.sort_index()
        pdt.assert_frame_equal(observed, expected)


def test_sequence_conflict_fixture_fails_with_order_stable_message() -> None:
    messages: list[str] = []
    evidence = _read_evidence("sequence_conflict.csv")

    for row_order in ((0, 1), (1, 0)):
        with pytest.raises(PhosPyInputError) as exc_info:
            _run_resolution(evidence.iloc[list(row_order), :].reset_index(drop=True))
        messages.append(str(exc_info.value))

    assert messages[0] == messages[1]
    assert "site_sequence values conflict" in messages[0]
    assert "site_id='AKT1;S473;'" in messages[0]
    assert "row order" in messages[0]


def test_mixed_valid_invalid_sequence_context_fixture_fails_clearly() -> None:
    with pytest.raises(PhosPyInputError) as exc_info:
        _run_resolution(_read_evidence("mixed_sequence_contexts.csv"))

    message = str(exc_info.value)
    assert "Mixed valid and invalid supplied evidence" in message
    assert "site_id='AKT1;S473;'" in message
    assert "valid_normalized_value='AAASAAA'" in message
    assert "expected='S'" in message
    assert "observed='T'" in message

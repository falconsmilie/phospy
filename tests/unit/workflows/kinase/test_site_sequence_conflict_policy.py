from __future__ import annotations

from typing import cast

import pandas as pd
import pytest

from phospy.contracts.configs import (
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
    KinaseSiteSequenceConflictPolicy,
)
from phospy.workflows.kinase.site_sequence_support import (
    KinaseSiteSequenceConflictError,
    KinaseSiteSequenceSupportBuilder,
)
from tests.support.site_keys import (
    protein_site_key,
    site_key_from_display_id,
    site_key_index_from_display_ids,
)

_DATASET_SEQUENCE = "AAAAAAASAAAAAAA"
_REFERENCE_SEQUENCE = "AAAAAAATAAAAAAA"
_MATCHING_SEQUENCE = "AAAAAAAYAAAAAAA"


def _dataset_frame(site_index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0 + index for index, _ in enumerate(site_index)],
            "sample_b": [2.0 + index for index, _ in enumerate(site_index)],
        },
        index=site_index.copy(),
    )


def _site_metadata(
    *,
    site_index: pd.Index,
    display_ids: list[str],
    sequences: list[object],
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            "site_sequence": sequences,
        },
        index=site_index.copy(),
    )


def _reference_sequences(sequences_by_display_id: dict[str, object]) -> pd.DataFrame:
    return pd.DataFrame(
        {"site_sequence": list(sequences_by_display_id.values())},
        index=pd.Index(list(sequences_by_display_id), name="site_id"),
    )


def _single_site_inputs(
    *,
    dataset_sequence: object = _DATASET_SEQUENCE,
    reference_sequence: object = _REFERENCE_SEQUENCE,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, str]:
    display_id = "GENE1;S10;"
    site_index = site_key_index_from_display_ids([display_id])
    return (
        _dataset_frame(site_index),
        _site_metadata(
            site_index=site_index,
            display_ids=[display_id],
            sequences=[dataset_sequence],
        ),
        _reference_sequences({display_id: reference_sequence}),
        str(site_index[0]),
    )


@pytest.mark.parametrize("policy", tuple(KinaseSiteSequenceConflictPolicy))
def test_matching_dataset_reference_sequence_passes_under_all_policies(
    policy: KinaseSiteSequenceConflictPolicy,
) -> None:
    dataset, site_metadata, references, site_key = _single_site_inputs(
        dataset_sequence=_MATCHING_SEQUENCE,
        reference_sequence=_MATCHING_SEQUENCE,
    )

    result = KinaseSiteSequenceSupportBuilder().run(
        dataset=dataset,
        site_metadata=site_metadata,
        reference_site_sequences=references,
        conflict_policy=policy,
    )

    assert result.dataset_reference_conflict_count == 0
    assert result.conflicts == ()
    assert result.site_sequences.at[site_key, "site_sequence"] == _MATCHING_SEQUENCE
    payload = result.diagnostics_payload()
    source_rows = cast(list[dict[str, object]], payload["selected_sequence_sources"])
    assert source_rows == [
        {
            "site_key": site_key,
            "display_id": "GENE1;S10;",
            "selected_sequence": _MATCHING_SEQUENCE,
            "selected_sequence_source": "reference",
            "policy": policy,
            "dataset_sequence": _MATCHING_SEQUENCE,
            "reference_sequence": _MATCHING_SEQUENCE,
            "diagnostic": (
                "reference sequence selected; dataset sequence matches "
                "reference sequence"
            ),
            "interpreter_version": result.interpreter_version,
        }
    ]


def test_conflicting_sequence_raises_under_error_policy() -> None:
    dataset, site_metadata, references, site_key = _single_site_inputs()

    with pytest.raises(KinaseSiteSequenceConflictError) as exc_info:
        KinaseSiteSequenceSupportBuilder().run(
            dataset=dataset,
            site_metadata=site_metadata,
            reference_site_sequences=references,
            conflict_policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
        )

    message = str(exc_info.value)
    assert site_key in message
    assert _DATASET_SEQUENCE in message
    assert _REFERENCE_SEQUENCE in message
    assert "conflict_policy=error" in message
    diagnostics = cast(
        list[dict[str, object]], exc_info.value.details["conflict_diagnostics"]
    )
    assert diagnostics[0]["site_key"] == site_key
    assert diagnostics[0]["dataset_sequence"] == _DATASET_SEQUENCE
    assert diagnostics[0]["reference_sequence"] == _REFERENCE_SEQUENCE
    assert diagnostics[0]["selected_sequence"] is None
    assert diagnostics[0]["selected_sequence_source"] == "unresolved"
    assert diagnostics[0]["policy"] == KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR


def test_conflicting_sequence_uses_dataset_under_prefer_dataset() -> None:
    dataset, site_metadata, references, site_key = _single_site_inputs()

    result = KinaseSiteSequenceSupportBuilder().run(
        dataset=dataset,
        site_metadata=site_metadata,
        reference_site_sequences=references,
        conflict_policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET,
    )

    assert result.site_sequences.at[site_key, "site_sequence"] == _DATASET_SEQUENCE
    conflict = result.conflicts[0]
    assert conflict.selected_sequence == _DATASET_SEQUENCE
    assert conflict.selected_sequence_source == "dataset"
    payload = result.diagnostics_payload()
    conflict_rows = cast(list[dict[str, object]], payload["conflict_diagnostics"])
    assert conflict_rows[0]["dataset_sequence"] == _DATASET_SEQUENCE
    assert conflict_rows[0]["reference_sequence"] == _REFERENCE_SEQUENCE
    assert conflict_rows[0]["selected_sequence_source"] == "dataset"
    source_rows = cast(list[dict[str, object]], payload["selected_sequence_sources"])
    assert source_rows[0]["site_key"] == site_key
    assert source_rows[0]["selected_sequence_source"] == "dataset"


def test_conflicting_sequence_uses_reference_under_prefer_reference() -> None:
    dataset, site_metadata, references, site_key = _single_site_inputs()

    result = KinaseSiteSequenceSupportBuilder().run(
        dataset=dataset,
        site_metadata=site_metadata,
        reference_site_sequences=references,
        conflict_policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
    )

    assert result.site_sequences.at[site_key, "site_sequence"] == _REFERENCE_SEQUENCE
    conflict = result.conflicts[0]
    assert conflict.selected_sequence == _REFERENCE_SEQUENCE
    assert conflict.selected_sequence_source == "reference"
    payload = result.diagnostics_payload()
    conflict_rows = cast(list[dict[str, object]], payload["conflict_diagnostics"])
    assert conflict_rows[0]["dataset_sequence"] == _DATASET_SEQUENCE
    assert conflict_rows[0]["reference_sequence"] == _REFERENCE_SEQUENCE
    assert conflict_rows[0]["selected_sequence_source"] == "reference"
    source_rows = cast(list[dict[str, object]], payload["selected_sequence_sources"])
    assert source_rows[0]["site_key"] == site_key
    assert source_rows[0]["selected_sequence_source"] == "reference"


def test_duplicate_display_ids_with_distinct_site_keys_do_not_collide() -> None:
    display_id = "GENE1;S10;"
    site_index = pd.Index(
        [
            protein_site_key(
                protein_identifier="GENE1",
                site="S10",
                protein_namespace="gene_symbol",
            ),
            protein_site_key(
                protein_identifier="GENE1",
                site="S10",
                protein_namespace="uniprot_accession",
            ),
        ],
        name="site_key",
    )
    dataset = _dataset_frame(site_index)
    site_metadata = _site_metadata(
        site_index=site_index,
        display_ids=[display_id, display_id],
        sequences=[_REFERENCE_SEQUENCE, _DATASET_SEQUENCE],
    )
    references = _reference_sequences({display_id: _REFERENCE_SEQUENCE})

    result = KinaseSiteSequenceSupportBuilder().run(
        dataset=dataset,
        site_metadata=site_metadata,
        reference_site_sequences=references,
        conflict_policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET,
    )

    assert result.dataset_reference_conflict_count == 1
    assert result.conflicts[0].site_key == str(site_index[1])
    assert result.site_sequences.at[str(site_index[0]), "site_sequence"] == (
        _REFERENCE_SEQUENCE
    )
    assert result.site_sequences.at[str(site_index[1]), "site_sequence"] == (
        _DATASET_SEQUENCE
    )
    assert result.display_reference_multi_matches == (
        {"display_id": display_id, "site_keys": tuple(site_index.astype(str))},
    )


def test_missing_reference_sequence_is_not_treated_as_conflict() -> None:
    display_id = "GENE1;S10;"
    site_key = site_key_from_display_id(display_id)
    site_index = pd.Index([site_key], name="site_key")
    dataset = _dataset_frame(site_index)
    site_metadata = _site_metadata(
        site_index=site_index,
        display_ids=[display_id],
        sequences=[_DATASET_SEQUENCE],
    )
    references = _reference_sequences({display_id: pd.NA})

    result = KinaseSiteSequenceSupportBuilder().run(
        dataset=dataset,
        site_metadata=site_metadata,
        reference_site_sequences=references,
        conflict_policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    )

    assert result.dataset_reference_conflict_count == 0
    assert result.conflicts == ()
    assert result.dataset_sequences_added == 1
    assert result.site_sequences.at[site_key, "site_sequence"] == _DATASET_SEQUENCE
    payload = result.diagnostics_payload()
    source_rows = cast(list[dict[str, object]], payload["selected_sequence_sources"])
    assert source_rows[0]["reference_sequence"] is None
    assert source_rows[0]["selected_sequence_source"] == "dataset"

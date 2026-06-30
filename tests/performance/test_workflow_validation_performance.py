from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import WorkflowValidationError
from phospy.validation.datasets.site_metadata import (
    SequenceContextContract,
    enforce_site_sequence_context_contract,
)
from tests.support.performance_contracts import (
    DATASET_VALIDATION_LARGE_N_SITES,
    WORKFLOW_VALIDATION_LARGE_SEQUENCE_RUNTIME_SECONDS_MAX,
    deterministic_analysis_ready_site_keys,
    deterministic_analysis_ready_site_metadata,
    measure_runtime_and_peak_mib,
)

pytestmark = [pytest.mark.performance, pytest.mark.release_gate]


@pytest.fixture(scope="module")
def large_kinase_sequence_metadata() -> pd.DataFrame:
    site_keys = deterministic_analysis_ready_site_keys(
        DATASET_VALIDATION_LARGE_N_SITES,
        start=50_000,
        gene_prefix="SEQPERF",
    )
    return deterministic_analysis_ready_site_metadata(
        site_keys,
        start=50_000,
        gene_prefix="SEQPERF",
        sequence_width=15,
    )


def test_large_kinase_sequence_validation_completes_under_generous_threshold(
    large_kinase_sequence_metadata: pd.DataFrame,
) -> None:
    source_by_site = dict.fromkeys(
        large_kinase_sequence_metadata.index.astype(str).tolist(),
        "reference",
    )

    _result, runtime_seconds, _peak_mib = measure_runtime_and_peak_mib(
        lambda: enforce_site_sequence_context_contract(
            site_metadata=large_kinase_sequence_metadata,
            field_name="performance.selected_site_sequences",
            workflow_name="kinase workflow request",
            scoring_mode="kinase_library_motif",
            contract=_kinase_library_sequence_contract(),
            error_type=WorkflowValidationError,
            sequence_source_by_site=source_by_site,
        ),
        warmup=False,
    )

    assert large_kinase_sequence_metadata.shape[0] == DATASET_VALIDATION_LARGE_N_SITES
    assert runtime_seconds < WORKFLOW_VALIDATION_LARGE_SEQUENCE_RUNTIME_SECONDS_MAX


def test_kinase_sequence_validation_remains_strict_for_invalid_sequences() -> None:
    site_keys = deterministic_analysis_ready_site_keys(
        2,
        start=60_000,
        gene_prefix="SEQSTRICT",
    )
    site_metadata = deterministic_analysis_ready_site_metadata(
        site_keys,
        start=60_000,
        gene_prefix="SEQSTRICT",
        sequence_width=15,
    )
    site_metadata.loc[site_keys[0], "site_sequence"] = "AAAAAAA*AAAAAAA"

    with pytest.raises(WorkflowValidationError) as exc_info:
        enforce_site_sequence_context_contract(
            site_metadata=site_metadata,
            field_name="performance.selected_site_sequences",
            workflow_name="kinase workflow request",
            scoring_mode="kinase_library_motif",
            contract=_kinase_library_sequence_contract(),
            error_type=WorkflowValidationError,
            sequence_source_by_site={
                str(site_key): "reference" for site_key in site_keys.astype(str)
            },
        )

    message = str(exc_info.value)
    assert "workflow-specific sequence context contract failed" in message
    assert "unsupported_characters='*'" in message
    assert str(site_keys[0]) in message


def _kinase_library_sequence_contract() -> SequenceContextContract:
    return SequenceContextContract(
        requires_site_sequence=True,
        requires_centered_site=True,
        required_window_length=15,
        center_index=7,
        allowed_residues=frozenset("ACDEFGHIKLMNPQRSTVWY"),
        allow_terminal_padding=False,
        allow_lowercase=False,
        allow_modified_residue_symbols=False,
        required_center_residues=frozenset({"S", "T", "Y"}),
        requires_known_sequence_source=True,
        contract_id="performance_kinase_library_fixed_15aa_centered_window",
    )

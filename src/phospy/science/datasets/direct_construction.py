"""Trusted direct dataset construction provenance helpers."""

from __future__ import annotations

import pandas as pd

from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_optional_table_strict
from phospy.provenance.models import (
    RunProvenance,
    TableFingerprint,
    TrustedDatasetConstructionAssertions,
)

DIRECT_CONSTRUCTION_WORKFLOW_NAME = "analysis_ready_dataset_direct_construction"
DIRECT_CONSTRUCTION_METHOD = "AnalysisReadyPhosphoDataset.__init__"
DIRECT_CONSTRUCTION_SOURCE = "direct_trusted_construction"
DIRECT_CONSTRUCTION_WARNING = (
    "Direct construction cannot prove biological correctness of caller-provided "
    "analysis-ready state."
)
DIRECT_CONSTRUCTION_DEPRECATION_WARNING = (
    "AnalysisReadyPhosphoDataset(...) direct construction is deprecated; use "
    "AnalysisReadyPhosphoDataset.from_trusted_tables(...) with "
    "TrustedDatasetConstructionAssertions for trusted construction, or "
    "AnalysisReadyDatasetBuilder.run(...) for ordinary dataset construction."
)
DIRECT_CONSTRUCTION_MISSING_ASSERTIONS_WARNING = (
    "No typed trusted construction assertion metadata was supplied; identity, "
    "intensity scale, quantitative meaning, aligned structure, localisation, "
    "sequence, and reference context are not recorded as user-asserted or "
    "waived."
)


def build_direct_construction_provenance(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    sample_metadata: pd.DataFrame | None,
    total: pd.DataFrame | None,
    comparisons: pd.DataFrame | None,
    imputation_observation_mask: pd.DataFrame | None,
    trusted_construction_assertions: TrustedDatasetConstructionAssertions | None,
) -> RunProvenance:
    """Build the minimal provenance marker for trusted direct construction."""

    assertions = (
        TrustedDatasetConstructionAssertions.missing()
        if trusted_construction_assertions is None
        else trusted_construction_assertions
    )
    table_entries = (
        ("dataset.phospho", phospho),
        ("dataset.site_metadata", site_metadata),
        ("dataset.sample_metadata", sample_metadata),
        ("dataset.total", total),
        ("dataset.comparisons", comparisons),
        ("dataset.imputation_observation_mask", imputation_observation_mask),
    )
    table_fingerprints = _fingerprint_direct_construction_tables(table_entries)
    construction_payload: dict[str, object] = {
        "method": DIRECT_CONSTRUCTION_METHOD,
        "source": DIRECT_CONSTRUCTION_SOURCE,
        "builder_used": False,
        "warning": DIRECT_CONSTRUCTION_WARNING,
        "trusted_construction_assertions": assertions.to_payload(),
        "trusted_construction_assertion_fingerprint": (
            assertions.assertion_fingerprint
        ),
        "trusted_assertion_metadata_provided": (assertions.assertion_metadata_provided),
        "missing_trusted_assertions": list(assertions.missing_assertions),
    }
    if not assertions.assertion_metadata_provided:
        construction_payload["assertion_warning"] = (
            DIRECT_CONSTRUCTION_MISSING_ASSERTIONS_WARNING
        )
    return RunProvenance(
        environment=collect_environment_provenance(),
        input_tables=table_fingerprints,
        preprocessing_stages=(),
        reference=None,
        workflow_name=DIRECT_CONSTRUCTION_WORKFLOW_NAME,
        workflow_parameters={"construction": construction_payload},
        random_state=None,
        random_seed_policy=None,
        output_tables=table_fingerprints,
    )


def _fingerprint_direct_construction_tables(
    entries: tuple[tuple[str, pd.DataFrame | None], ...],
) -> tuple[TableFingerprint, ...]:
    fingerprints: list[TableFingerprint] = []
    for name, table in entries:
        fingerprint = fingerprint_optional_table_strict(table, name=name)
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
    return tuple(fingerprints)


__all__ = [
    "DIRECT_CONSTRUCTION_METHOD",
    "DIRECT_CONSTRUCTION_DEPRECATION_WARNING",
    "DIRECT_CONSTRUCTION_MISSING_ASSERTIONS_WARNING",
    "DIRECT_CONSTRUCTION_SOURCE",
    "DIRECT_CONSTRUCTION_WARNING",
    "DIRECT_CONSTRUCTION_WORKFLOW_NAME",
    "build_direct_construction_provenance",
]

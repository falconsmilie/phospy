from __future__ import annotations

import inspect
from dataclasses import replace

import pandas as pd
import pytest

from phospy.api import KinaseWorkflowRequest, Organism, ReferenceBundle
from phospy.api.configs import KinaseScoringConfig
from phospy.errors.validation import PhosPyValidationError, WorkflowValidationError
from phospy.provenance.models import ReferenceProvenance
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.references.models import ReferenceContext
from phospy.validation.identity_contracts import (
    REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE,
    ReferenceContextCompatibilityWarning,
    validate_reference_context_compatibility,
)
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)


def _context(**overrides: object) -> ReferenceContext:
    values = {
        "organism": "rat",
        "protein_namespace": "gene_symbol",
        "source_name": "unit-reference",
        "source_version": "v1",
        "proteome_version": None,
        "reference_table_sha256": "a" * 64,
    }
    values.update(overrides)
    return ReferenceContext(**values)


def _dataset(context: ReferenceContext | None) -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;"]
    site_ids = site_key_index_from_display_ids(
        display_ids,
        protein_namespace="gene_symbol",
    )
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=pd.DataFrame({"sample_a": [1.0], "sample_b": [1.1]}, index=site_ids),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_ids.astype(str).tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_ids),
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["AAAAAAAYAAAAAAA"],
            },
            index=site_ids.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    if context is None:
        return dataset
    provenance = dataset.provenance
    if provenance is None:
        raise AssertionError(
            "analysis-ready dataset must carry construction provenance"
        )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=dataset.phospho,
        site_metadata=dataset.site_metadata,
        organism=dataset.organism,
        intensity_scale_state=dataset.intensity_scale_state,
        processing_state=dataset.processing_state,
        provenance=replace(provenance, reference_context=context),
    )


def _references(context: ReferenceContext | None) -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["MAP2K6"], "substrate_site": ["MAPK14;Y182;"]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["AAAAAAAYAAAAAAA"]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        provenance=ReferenceProvenance(
            source_type="explicit",
            organism=Organism.RAT.value,
            bundle_id=None,
            source_name="unit-reference",
            source_version=context.source_version,
            identifier_namespace="gene_symbol",
            table_fingerprints=(),
            reference_context=context,
        ),
    )


def test_matching_reference_contexts_pass() -> None:
    context = _context()

    assert (
        validate_reference_context_compatibility(
            context,
            _context(),
            operation="unit matching context",
        )
        is None
    )


def test_mismatched_source_version_fails() -> None:
    with pytest.raises(PhosPyValidationError) as exc_info:
        validate_reference_context_compatibility(
            _context(source_version="v1"),
            _context(source_version="v2"),
            operation="unit source version mismatch",
        )

    message = str(exc_info.value)
    assert "source_version" in message
    assert "unit source version mismatch" in message


def test_mismatched_namespace_fails() -> None:
    with pytest.raises(PhosPyValidationError) as exc_info:
        validate_reference_context_compatibility(
            _context(protein_namespace="gene_symbol"),
            _context(protein_namespace="uniprot_accession"),
            operation="unit namespace mismatch",
        )

    assert "protein_namespace" in str(exc_info.value)


def test_reference_context_validator_unknown_fails_by_default() -> None:
    with pytest.raises(PhosPyValidationError, match="unknown reference context"):
        validate_reference_context_compatibility(
            _context(),
            None,
            operation="unit unknown context",
        )


def test_unknown_context_passes_only_with_explicit_override_warning() -> None:
    warning = validate_reference_context_compatibility(
        _context(),
        None,
        operation="unit allowed unknown context",
        allow_unknown=True,
    )

    assert isinstance(warning, ReferenceContextCompatibilityWarning)
    assert warning.code == REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE
    assert warning.severity == "warning"
    assert warning.missing_contexts == ("right",)
    payload = warning.to_payload()
    assert payload["operation"] == "unit allowed unknown context"
    assert payload["missing_contexts"] == ["right"]


def test_kinase_workflow_validator_rejects_mismatched_reference_context() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(_context(source_version="dataset-v1")),
        references=_references(_context(source_version="reference-v2")),
        scoring_config=KinaseScoringConfig(min_substrates=2),
        activity_config=None,
    )

    with pytest.raises(WorkflowValidationError) as exc_info:
        KinaseWorkflowValidator().run(request)

    message = str(exc_info.value)
    assert "reference-context compatibility failed" in message
    assert "source_version" in message


def test_kinase_workflow_validator_uses_shared_reference_context_validator() -> None:
    source = inspect.getsource(KinaseWorkflowValidator.run)

    assert "validate_reference_context_compatibility(" in source

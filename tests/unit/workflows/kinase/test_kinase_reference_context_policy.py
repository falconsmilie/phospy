from __future__ import annotations

from dataclasses import replace

import pandas as pd
import pytest

from phospy.advanced.configs import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api import KinaseWorkflow, KinaseWorkflowRequest, Organism, ReferenceBundle
from phospy.errors.validation import WorkflowValidationError
from phospy.provenance.models import ReferenceProvenance
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.references.models import ReferenceContext
from phospy.validation.identity_contracts import REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE
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


def _window(display_id: str) -> str:
    residue = display_id.split(";")[1][0].upper()
    return ("A" * 15) + residue + ("A" * 15)


def _dataset(context: ReferenceContext | None) -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "MAPK1;T185;"]
    site_ids = site_key_index_from_display_ids(
        display_ids,
        protein_namespace="gene_symbol",
    )
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 1.5], "sample_b": [1.1, 1.7]},
            index=site_ids.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_ids.astype(str).tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_ids),
                "gene_symbol": ["MAPK14", "MAPK1"],
                "protein_id": ["MAPK14", "MAPK1"],
                "site": ["Y182", "T185"],
                "site_sequence": [_window(display_id) for display_id in display_ids],
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
    display_ids = ["MAPK14;Y182;", "MAPK1;T185;"]
    provenance = None
    if context is not None:
        provenance = ReferenceProvenance(
            source_type="explicit",
            organism=Organism.RAT.value,
            bundle_id=None,
            source_name="unit-reference",
            source_version=context.source_version,
            identifier_namespace="gene_symbol",
            table_fingerprints=(),
            reference_context=context,
        )
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["MAP2K6", "MAP2K6"], "substrate_site": display_ids}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": [_window(display_id) for display_id in display_ids]},
            index=pd.Index(display_ids, name="site_id"),
        ),
        provenance=provenance,
    )


def _request(
    *,
    dataset_context: ReferenceContext | None,
    reference_context: ReferenceContext | None,
    policy: ReferenceContextCompatibilityPolicy = (
        ReferenceContextCompatibilityPolicy.REQUIRE_KNOWN_MATCH
    ),
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=_dataset(dataset_context),
        references=_references(reference_context),
        scoring_config=KinaseScoringConfig(
            reliability_profile="custom",
            min_substrates=2,
            reference_context_compatibility_policy=policy,
        ),
        prediction_config=KinasePredictionConfig(
            top_k=1,
            deterministic_max_selected_kinases=1,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
    )


def _reference_context_caveats(result: object):
    return [
        caveat
        for caveat in result.caveats
        if caveat.code == REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE
    ]


def test_kinase_workflow_unknown_reference_context_fails_by_default() -> None:
    request = _request(dataset_context=_context(), reference_context=None)

    with pytest.raises(WorkflowValidationError, match="unknown reference context"):
        KinaseWorkflow().run(request)


def test_kinase_workflow_unknown_reference_context_passes_with_explicit_policy_and_caveat() -> (
    None
):
    result = KinaseWorkflow().run(
        _request(
            dataset_context=_context(),
            reference_context=None,
            policy=ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT,
        )
    )

    caveats = _reference_context_caveats(result)
    assert len(caveats) == 1
    caveat = caveats[0]
    assert caveat.severity == "warning"
    assert caveat.details["policy"] == "allow_unknown_with_caveat"
    assert caveat.details["workflow_scope"] == "kinase_scoring"
    assert caveat.details["missing_contexts"] == ("right",)
    assert caveat.details["operation"] == (
        "kinase workflow result dataset/reference bundle"
    )
    assert isinstance(caveat.details["left_reference_context_id"], str)
    assert caveat.details["right_reference_context_id"] is None


def test_kinase_workflow_mismatched_reference_context_fails_even_with_unknown_policy() -> (
    None
):
    request = _request(
        dataset_context=_context(source_version="dataset-v1"),
        reference_context=_context(source_version="reference-v2"),
        policy=ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT,
    )

    with pytest.raises(WorkflowValidationError, match="source_version"):
        KinaseWorkflow().run(request)


def test_kinase_workflow_known_reference_context_passes_without_caveat() -> None:
    context = _context()

    result = KinaseWorkflow().run(
        _request(dataset_context=context, reference_context=_context())
    )

    assert _reference_context_caveats(result) == []

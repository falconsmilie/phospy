from __future__ import annotations

from dataclasses import replace

import pandas as pd
import pytest

from phospy.api import (
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    SignalomeWorkflow,
    SignalomeWorkflowRequest,
)
from phospy.api.configs import ReferenceContextCompatibilityPolicy
from phospy.api.results import KinasePredictionResult, KinaseScoringResult
from phospy.errors.validation import WorkflowValidationError
from phospy.provenance.models import (
    EnvironmentProvenance,
    ReferenceProvenance,
    RunProvenance,
)
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
from tests.support.signalome_config import build_signalome_config
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


def _run_provenance(
    *,
    workflow_name: str,
    context: ReferenceContext | None,
) -> RunProvenance:
    return RunProvenance(
        environment=EnvironmentProvenance(
            package_name="phospy",
            package_version="test",
            python_version="3.13",
            dependency_versions={},
        ),
        input_tables=(),
        preprocessing_stages=(),
        reference=None,
        workflow_name=workflow_name,
        workflow_parameters={},
        random_state=None,
        random_seed_policy=None,
        output_tables=(),
        reference_context=context,
    )


def _window(display_id: str) -> str:
    residue = display_id.split(";")[1][0].upper()
    return ("A" * 15) + residue + ("A" * 15)


def _dataset(context: ReferenceContext | None) -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "MAPK1;T185;", "JUN;S63;"]
    site_ids = site_key_index_from_display_ids(
        display_ids,
        protein_namespace="gene_symbol",
    )
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=pd.DataFrame(
            {
                "sample_a": [1.0, 1.5, 2.0],
                "sample_b": [1.1, 1.7, 2.4],
            },
            index=site_ids.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_ids.astype(str).tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_ids),
                "gene_symbol": ["MAPK14", "MAPK1", "JUN"],
                "protein_id": ["MAPK14", "MAPK1", "JUN"],
                "site": ["Y182", "T185", "S63"],
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


def _references(
    *,
    context: ReferenceContext | None,
    display_ids: list[str],
) -> ReferenceBundle:
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
            {
                "kinase": ["MAP2K6", "MAP2K6", "MAPK8"],
                "substrate_site": display_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": [_window(display_id) for display_id in display_ids]},
            index=pd.Index(display_ids, name="site_id"),
        ),
        provenance=provenance,
    )


def _matrix(
    *,
    values: list[list[float]],
    site_ids: list[str],
) -> pd.DataFrame:
    return pd.DataFrame(
        values,
        index=pd.Index(site_ids, name="site_key"),
        columns=pd.Index(["MAP2K6", "MAPK8"], name="kinase"),
        dtype=float,
    )


def _kinase_result(
    *,
    dataset_context: ReferenceContext | None,
    result_context: ReferenceContext | None,
    reference_context: ReferenceContext | None,
) -> KinaseWorkflowResult:
    dataset = _dataset(dataset_context)
    site_ids = dataset.phospho.index.astype(str).tolist()
    display_ids = dataset.site_metadata.loc[:, "display_id"].astype(str).tolist()
    prediction_matrix = _matrix(
        values=[[0.9, 0.1], [0.3, 0.8], [0.7, 0.2]],
        site_ids=site_ids,
    )
    score_matrix = _matrix(
        values=[[1.0, 0.1], [0.2, 1.0], [0.8, 0.3]],
        site_ids=site_ids,
    )
    return KinaseWorkflowResult(
        dataset=dataset,
        references=_references(context=reference_context, display_ids=display_ids),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            rank_weighted_fusion_scores=score_matrix,
        ),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        provenance=(
            None
            if result_context is None
            else _run_provenance(
                workflow_name="kinase_workflow",
                context=result_context,
            )
        ),
        activity_result=None,
    )


def _request(
    *,
    dataset_context: ReferenceContext | None,
    result_context: ReferenceContext | None,
    reference_context: ReferenceContext | None,
    policy: ReferenceContextCompatibilityPolicy = (
        ReferenceContextCompatibilityPolicy.REQUIRE_KNOWN_MATCH
    ),
) -> SignalomeWorkflowRequest:
    return SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset_context=dataset_context,
            result_context=result_context,
            reference_context=reference_context,
        ),
        config=build_signalome_config(
            module_count=1,
            reference_context_compatibility_policy=policy,
        ),
    )


def _reference_context_caveats(result: object):
    return [
        caveat
        for caveat in result.caveats
        if caveat.code == REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE
    ]


def test_signalome_workflow_unknown_reference_context_fails_by_default() -> None:
    request = _request(
        dataset_context=_context(),
        result_context=None,
        reference_context=None,
    )

    with pytest.raises(WorkflowValidationError, match="unknown reference context"):
        SignalomeWorkflow().run(request)


def test_signalome_workflow_unknown_reference_context_passes_with_explicit_policy_and_caveat() -> (
    None
):
    result = SignalomeWorkflow().run(
        _request(
            dataset_context=_context(),
            result_context=None,
            reference_context=None,
            policy=ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT,
        )
    )

    caveats = _reference_context_caveats(result)
    assert len(caveats) == 2
    operations = {caveat.details["operation"] for caveat in caveats}
    assert operations == {
        "signalome workflow result dataset/upstream kinase result",
        "signalome workflow result dataset/upstream kinase reference",
    }
    for caveat in caveats:
        assert caveat.severity == "warning"
        assert caveat.details["policy"] == "allow_unknown_with_caveat"
        assert caveat.details["workflow_scope"] == "signalome"
        assert caveat.details["missing_contexts"] == ("right",)
        assert isinstance(caveat.details["left_reference_context_id"], str)
        assert caveat.details["right_reference_context_id"] is None


def test_signalome_workflow_mismatched_reference_context_fails_even_with_unknown_policy() -> (
    None
):
    request = _request(
        dataset_context=_context(source_version="dataset-v1"),
        result_context=_context(source_version="result-v2"),
        reference_context=None,
        policy=ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT,
    )

    with pytest.raises(WorkflowValidationError, match="source_version"):
        SignalomeWorkflow().run(request)


def test_signalome_workflow_known_reference_context_passes_without_caveat() -> None:
    context = _context()

    result = SignalomeWorkflow().run(
        _request(
            dataset_context=context,
            result_context=_context(),
            reference_context=_context(),
        )
    )

    assert _reference_context_caveats(result) == []

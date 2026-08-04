from __future__ import annotations

from dataclasses import replace

import pandas as pd
import pytest

from phospy.api import (
    KinasePredictionConfig,
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
)
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseExecutionConfig,
    ValidatedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.kinase.resolved_validator import (
    ResolvedKinaseEligibilityValidator,
    ResolvedKinaseInputs,
)
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


def _window(display_id: str) -> str:
    residue = display_id.split(";")[1][0].upper()
    return ("A" * 15) + residue + ("A" * 15)


def _validated_request(
    request: KinaseWorkflowRequest,
) -> ValidatedKinaseWorkflowRequest:
    return ValidatedKinaseWorkflowRequest(
        request=request,
        dataset_view=DatasetInternalView(request.dataset),
    )


def _dataset(
    display_ids: tuple[str, ...] = ("KIN1;S1;", "KIN2;T2;"),
) -> AnalysisReadyPhosphoDataset:
    site_index = site_key_index_from_display_ids(list(display_ids))
    return trusted_analysis_ready_dataset_from_tables(
        phospho=pd.DataFrame(
            {
                "sample_a": [1.0 + index for index, _ in enumerate(display_ids)],
                "sample_b": [2.0 + index for index, _ in enumerate(display_ids)],
            },
            index=site_index.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_index.astype(str).tolist(),
                "display_id": list(display_ids),
                **site_key_context_columns(site_index),
                "gene_symbol": [
                    display_id.split(";", 1)[0] for display_id in display_ids
                ],
                "protein_id": [
                    display_id.split(";", 1)[0] for display_id in display_ids
                ],
                "site": [display_id.split(";")[1] for display_id in display_ids],
                "site_sequence": [_window(display_id) for display_id in display_ids],
            },
            index=site_index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references(display_ids: tuple[str, ...] = ("KIN1;S1;", "KIN2;T2;")):
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K_REF" for _ in display_ids],
                "substrate_site": list(display_ids),
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": [_window(display_id) for display_id in display_ids]},
            index=pd.Index(list(display_ids), name="site_id"),
        ),
    )


def _site_identity_map(dataset: AnalysisReadyPhosphoDataset) -> pd.DataFrame:
    site_keys = dataset.phospho.index.astype(str).tolist()
    return pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": dataset.site_metadata.loc[:, "display_id"]
            .astype(str)
            .tolist(),
        },
        index=pd.Index(site_keys, name=dataset.phospho.index.name),
    )


def _site_sequences(
    dataset: AnalysisReadyPhosphoDataset,
    *,
    supported_site_keys: tuple[str, ...] | None = None,
) -> pd.DataFrame:
    identity_map = _site_identity_map(dataset)
    supported = list(supported_site_keys or tuple(identity_map.index.astype(str)))
    sequences = pd.DataFrame(
        {
            "site_sequence": [
                _window(str(identity_map.at[site_key, "display_id"]))
                for site_key in supported
            ],
            "display_id": [
                str(identity_map.at[site_key, "display_id"]) for site_key in supported
            ],
        },
        index=pd.Index(supported, name=dataset.phospho.index.name),
    )
    return sequences


def _execution_config(*, min_substrates: int = 2) -> ResolvedKinaseExecutionConfig:
    return ResolvedKinaseExecutionConfig(
        scoring_min_substrates=int(min_substrates),
        include_diagnostic_scoring_tables=False,
        profile_missing_value_strategy="strict",
        prediction_top_k=2,
        prediction_deterministic_max_selected_kinases=2,
        prediction_adaptive_ensemble_runs=2,
        prediction_mode="deterministic_ranking",
        prediction_adaptive_policy="stable",
        prediction_n_iterations=5,
        prediction_random_state=None,
    )


def _resolved_inputs(
    *,
    dataset: AnalysisReadyPhosphoDataset | None = None,
    kinase_substrate_map: pd.DataFrame | None = None,
    site_sequences: pd.DataFrame | None = None,
    min_substrates: int = 2,
    reference_site_count: int = 2,
) -> ResolvedKinaseInputs:
    dataset = dataset or _dataset()
    dataset_phospho = dataset.phospho
    site_sequences = (
        site_sequences if site_sequences is not None else _site_sequences(dataset)
    )
    scoring_site_index = site_sequences.index.copy()
    return ResolvedKinaseInputs(
        dataset=dataset,
        dataset_phospho=dataset_phospho,
        references=_references(),
        kinase_substrate_map=(
            kinase_substrate_map
            if kinase_substrate_map is not None
            else pd.DataFrame(
                {
                    "kinase": ["K_REF", "K_REF"],
                    "substrate_site": dataset_phospho.index.astype(str).tolist(),
                }
            )
        ),
        site_sequences=site_sequences,
        site_identity_map=_site_identity_map(dataset),
        scoring_site_index=scoring_site_index,
        activity_phospho_matrix=dataset_phospho.reindex(index=scoring_site_index),
        execution_config=_execution_config(min_substrates=min_substrates),
        reference_site_count=reference_site_count,
    )


def test_resolved_kinase_validator_rejects_no_reference_overlap() -> None:
    resolved_inputs = _resolved_inputs(
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["K_REF"], "substrate_site": ["unmatched_site_key"]}
        ),
        reference_site_count=1,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        ResolvedKinaseEligibilityValidator().run(resolved_inputs)

    assert exc_info.value.seam == "kinase.interpreter.reference_coverage"
    assert exc_info.value.details["overlap_sites"] == 0


def test_resolved_kinase_validator_rejects_no_eligible_kinases() -> None:
    dataset = _dataset()
    site_key = str(dataset.phospho.index[0])
    resolved_inputs = _resolved_inputs(
        dataset=dataset,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["K_REF"], "substrate_site": [site_key]}
        ),
        reference_site_count=1,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        ResolvedKinaseEligibilityValidator().run(resolved_inputs)

    assert exc_info.value.seam == "kinase.interpreter.eligible_kinases"
    assert exc_info.value.details["eligible_kinases"] == 0
    assert exc_info.value.details["max_quantified_sites_per_kinase"] == 1


def test_resolved_kinase_validator_rejects_no_sequence_supported_sites() -> None:
    dataset = _dataset()
    empty_sequences = pd.DataFrame(
        columns=["site_sequence", "display_id"],
        index=pd.Index([], name=dataset.phospho.index.name),
    )
    resolved_inputs = _resolved_inputs(
        dataset=dataset,
        site_sequences=empty_sequences,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        ResolvedKinaseEligibilityValidator().run(resolved_inputs)

    assert exc_info.value.seam == "kinase.interpreter.sequence_support"
    assert exc_info.value.details["sequence_supported_sites"] == 0


def test_resolved_kinase_validator_enforces_production_reliability_profile() -> None:
    resolved_inputs = _resolved_inputs()
    production_labelled_config = replace(
        resolved_inputs.execution_config,
        requested_reliability_profile=KinaseReliabilityProfile.PRODUCTION,
        effective_reliability_profile=KinaseReliabilityProfile.PRODUCTION,
    )
    resolved_inputs = replace(
        resolved_inputs,
        execution_config=production_labelled_config,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        ResolvedKinaseEligibilityValidator().run(resolved_inputs)

    assert exc_info.value.seam == "kinase.interpreter.production_reliability_profile"
    assert "min_substrates must be at least 5" in str(exc_info.value)
    assert exc_info.value.details["effective_reliability_profile"] == "production"
    assert exc_info.value.details["scoring_config_min_substrates"] == 2


def test_kinase_interpreter_delegates_resolved_validation() -> None:
    calls: list[ResolvedKinaseInputs] = []

    class _ResolvedValidatorSpy:
        def run(self, resolved_inputs: ResolvedKinaseInputs) -> None:
            calls.append(resolved_inputs)
            ResolvedKinaseEligibilityValidator().run(resolved_inputs)

    dataset = _dataset()
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=_references(),
        scoring_config=KinaseScoringConfig(
            reliability_profile="custom",
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
    )

    interpreted = KinaseWorkflowInterpreter(
        resolved_validator=_ResolvedValidatorSpy()  # type: ignore[arg-type]
    ).run(_validated_request(request))

    assert len(calls) == 1
    resolved_inputs = calls[0]
    assert isinstance(resolved_inputs, ResolvedKinaseInputs)
    assert resolved_inputs.references is interpreted.references
    assert interpreted.attrition_metrics is resolved_inputs.attrition_metrics
    assert not interpreted.kinase_substrate_map.empty

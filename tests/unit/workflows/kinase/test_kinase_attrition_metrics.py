from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import (
    KinaseAttritionPolicy,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
)
from phospy.errors.workflows import PhosPyWorkflowError, WorkflowBoundaryError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.workflows.kinase.attrition_metrics import (
    KinaseAttritionMetrics,
    build_kinase_attrition_provenance_payload,
)
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
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
from tests.support.unsafe_dataset_states import (
    unsafe_set_dataset_site_metadata_columns,
)


def _window(display_id: str) -> str:
    residue = display_id.split(";")[1][0].upper()
    return ("A" * 15) + residue + ("A" * 15)


def _dataset(
    display_ids: list[str],
    *,
    missing_sequence_display_ids: set[str] | None = None,
) -> AnalysisReadyPhosphoDataset:
    site_index = site_key_index_from_display_ids(display_ids)
    dataset = trusted_analysis_ready_dataset_from_tables(
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
                "display_id": display_ids,
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
    missing = missing_sequence_display_ids or set()
    if missing:
        unsafe_set_dataset_site_metadata_columns(
            dataset,
            {
                "site_sequence": [
                    "" if display_id in missing else _window(display_id)
                    for display_id in display_ids
                ]
            },
        )
    return dataset


def _references(
    *,
    substrate_display_ids: list[str],
    sequence_display_ids: list[str],
) -> ReferenceBundle:
    construction_sequence_ids = list(
        dict.fromkeys([*substrate_display_ids, *sequence_display_ids])
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["KINASE_A" for _ in substrate_display_ids],
                "substrate_site": substrate_display_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    _window(display_id) for display_id in construction_sequence_ids
                ]
            },
            index=pd.Index(construction_sequence_ids, name="site_id"),
        ),
    )
    # Keep public ReferenceBundle construction valid, then model degraded
    # interpreter-visible sequence support for focused attrition tests.
    object.__setattr__(
        references,
        "site_sequences",
        pd.DataFrame(
            {
                "site_sequence": [
                    _window(display_id) for display_id in sequence_display_ids
                ]
            },
            index=pd.Index(sequence_display_ids, name="site_id"),
        ),
    )
    return references


def _request(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    references: ReferenceBundle,
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
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


def test_kinase_attrition_metrics_no_attrition() -> None:
    display_ids = ["KIN1;S1;", "KIN2;T2;"]
    interpreted = KinaseWorkflowInterpreter().run(
        _request(
            dataset=_dataset(display_ids),
            references=_references(
                substrate_display_ids=display_ids,
                sequence_display_ids=display_ids,
            ),
        )
    )

    metrics = interpreted.attrition_metrics

    assert metrics is not None
    assert metrics.total_dataset_sites == 2
    assert metrics.reference_overlap_sites == 2
    assert metrics.sequence_supported_sites == 2
    assert metrics.scored_sites == 2
    assert metrics.reference_overlap_fraction == pytest.approx(1.0)
    assert metrics.sequence_supported_fraction == pytest.approx(1.0)
    assert metrics.scored_fraction == pytest.approx(1.0)


def test_kinase_attrition_metrics_reference_overlap_loss() -> None:
    display_ids = ["KIN1;S1;", "KIN2;T2;", "KIN3;Y3;"]
    reference_ids = display_ids[:2]
    interpreted = KinaseWorkflowInterpreter().run(
        _request(
            dataset=_dataset(display_ids),
            references=_references(
                substrate_display_ids=reference_ids,
                sequence_display_ids=reference_ids,
            ),
        )
    )

    metrics = interpreted.attrition_metrics

    assert metrics is not None
    assert metrics.total_dataset_sites == 3
    assert metrics.reference_overlap_sites == 2
    assert metrics.sequence_supported_sites == 3
    assert metrics.scored_sites == 2
    assert metrics.reference_overlap_fraction == pytest.approx(2 / 3)
    assert metrics.sequence_supported_fraction == pytest.approx(1.0)
    assert metrics.scored_fraction == pytest.approx(2 / 3)


def test_kinase_attrition_metrics_sequence_support_loss() -> None:
    display_ids = ["KIN1;S1;", "KIN2;T2;", "KIN3;Y3;"]
    sequence_supported_ids = display_ids[:2]
    interpreted = KinaseWorkflowInterpreter().run(
        _request(
            dataset=_dataset(
                display_ids,
                missing_sequence_display_ids={display_ids[2]},
            ),
            references=_references(
                substrate_display_ids=display_ids,
                sequence_display_ids=sequence_supported_ids,
            ),
        )
    )

    metrics = interpreted.attrition_metrics

    assert metrics is not None
    assert metrics.total_dataset_sites == 3
    assert metrics.reference_overlap_sites == 3
    assert metrics.sequence_supported_sites == 2
    assert metrics.scored_sites == 2
    assert metrics.reference_overlap_fraction == pytest.approx(1.0)
    assert metrics.sequence_supported_fraction == pytest.approx(2 / 3)
    assert metrics.scored_fraction == pytest.approx(2 / 3)


def test_kinase_attrition_metrics_zero_scored_sites_is_explicit() -> None:
    display_ids = ["KIN1;S1;", "KIN2;T2;", "KIN3;Y3;", "KIN4;S4;"]
    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflowInterpreter().run(
            _request(
                dataset=_dataset(
                    display_ids,
                    missing_sequence_display_ids=set(display_ids[:2]),
                ),
                references=_references(
                    substrate_display_ids=display_ids[:2],
                    sequence_display_ids=display_ids[2:],
                ),
            )
        )

    error = exc_info.value
    assert error.seam == "kinase.interpreter.scored_site_support"
    assert error.details["total_dataset_sites"] == 4
    assert error.details["reference_overlap_sites"] == 2
    assert error.details["sequence_supported_sites"] == 2
    assert error.details["scored_sites"] == 0
    assert error.details["scored_fraction"] == pytest.approx(0.0)


def test_kinase_attrition_metrics_zero_denominator_fails_explicitly() -> None:
    with pytest.raises(PhosPyWorkflowError, match="total_dataset_sites"):
        KinaseAttritionMetrics.from_counts(
            total_dataset_sites=0,
            reference_overlap_sites=0,
            sequence_supported_sites=0,
            scored_sites=0,
        )


def test_kinase_attrition_payload_preserves_calculated_metrics_and_policy() -> None:
    metrics = KinaseAttritionMetrics.from_counts(
        total_dataset_sites=4,
        reference_overlap_sites=3,
        sequence_supported_sites=2,
        scored_sites=1,
    )
    policy = KinaseAttritionPolicy(
        minimum_reference_overlap_fraction=0.25,
        minimum_sequence_supported_fraction=0.25,
        minimum_scored_fraction=0.25,
        on_violation="warn",
    )

    payload = build_kinase_attrition_provenance_payload(
        metrics=metrics,
        policy=policy,
        violations=(),
    )

    assert payload["metrics"] == metrics.to_payload()
    assert payload["policy"] == {
        "minimum_reference_overlap_fraction": 0.25,
        "minimum_sequence_supported_fraction": 0.25,
        "minimum_scored_fraction": 0.25,
        "on_violation": "warn",
    }
    assert payload["policy_outcome"] == "passed"

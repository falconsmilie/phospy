from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import pandas as pd
import pytest

import phospy.science.datasets.internal_frame_store as internal_frame_store_module
import phospy.validation.workflows.differential as differential_validation_module
import phospy.validation.workflows.differential_design_rules as design_rules_module
import phospy.validation.workflows.quantitative as quantitative_module
import phospy.workflows.differential.replicates as replicates_module
import phospy.workflows.differential.validator as validator_module
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.validation.workflows.differential import (
    DifferentialDatasetEligibilityValidator,
    ExperimentalDesignContractValidator,
)
from phospy.workflows.differential.models import (
    DifferentialDatasetEligibilityValidatorContract,
    DifferentialDesignValidatorContract,
    DifferentialTechnicalReplicatePlannerContract,
)
from phospy.workflows.differential.reliability import (
    resolved_minimum_condition_replicates,
)
from phospy.workflows.differential.replicates import (
    TechnicalReplicateAggregationPlanner,
)
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from tests.support.workflow_identity_coherence import (
    build_duplicate_display_differential_request,
)


@dataclass
class _DatasetViewEvents:
    views: list[DatasetInternalView] = field(default_factory=list)
    dataframe_snapshots: Counter[str] = field(default_factory=Counter)
    optional_dataframe_snapshots: Counter[str] = field(default_factory=Counter)


@pytest.fixture
def instrument_differential_dataset_view(
    monkeypatch: pytest.MonkeyPatch,
) -> _DatasetViewEvents:
    events = _DatasetViewEvents()
    original_view_class = DatasetInternalView
    original_dataframe_snapshot = (
        internal_frame_store_module.immutable_dataframe_snapshot
    )
    original_optional_dataframe_snapshot = (
        internal_frame_store_module.immutable_optional_dataframe_snapshot
    )

    class _CountingDatasetInternalView(original_view_class):
        def __init__(self, dataset: AnalysisReadyPhosphoDataset) -> None:
            events.views.append(self)
            super().__init__(dataset)

    def _counting_dataframe_snapshot(
        value: pd.DataFrame,
        *,
        field_name: str,
        error_type: type[Exception] = TypeError,
    ):
        events.dataframe_snapshots.update((field_name,))
        return original_dataframe_snapshot(
            value,
            field_name=field_name,
            error_type=error_type,
        )

    def _counting_optional_dataframe_snapshot(
        value: pd.DataFrame | None,
        *,
        field_name: str,
        error_type: type[Exception] = TypeError,
    ):
        if isinstance(value, pd.DataFrame):
            events.optional_dataframe_snapshots.update((field_name,))
        return original_optional_dataframe_snapshot(
            value,
            field_name=field_name,
            error_type=error_type,
        )

    for module in (
        validator_module,
        differential_validation_module,
        design_rules_module,
        quantitative_module,
        replicates_module,
    ):
        monkeypatch.setattr(
            module,
            "DatasetInternalView",
            _CountingDatasetInternalView,
        )
    monkeypatch.setattr(
        internal_frame_store_module,
        "immutable_dataframe_snapshot",
        _counting_dataframe_snapshot,
    )
    monkeypatch.setattr(
        internal_frame_store_module,
        "immutable_optional_dataframe_snapshot",
        _counting_optional_dataframe_snapshot,
    )
    return events


def test_production_differential_collaborators_satisfy_protocol_calls() -> None:
    request = build_duplicate_display_differential_request()
    config = request.config
    dataset_view = DatasetInternalView(request.dataset)

    eligibility_validator: DifferentialDatasetEligibilityValidatorContract = (
        DifferentialDatasetEligibilityValidator()
    )
    technical_replicate_planner: DifferentialTechnicalReplicatePlannerContract = (
        TechnicalReplicateAggregationPlanner()
    )
    design_validator: DifferentialDesignValidatorContract = (
        ExperimentalDesignContractValidator()
    )

    eligibility_validator.run(
        dataset=request.dataset,
        imputed_value_policy=config.imputed_value_policy,
        allow_suspicious_declared_input_scale=(
            config.allow_suspicious_declared_input_scale
        ),
        dataset_view=dataset_view,
    )
    technical_replicate_plan = technical_replicate_planner.run(
        dataset=request.dataset,
        design=request.design,
        technical_replicate_policy=config.technical_replicate_policy,
        dataset_view=dataset_view,
    )
    validated_design = design_validator.run(
        dataset=request.dataset,
        design=request.design,
        contrasts=request.contrasts,
        allow_design_subset=config.allow_design_subset,
        minimum_condition_replicates=resolved_minimum_condition_replicates(config),
        paired_design_policy=config.paired_design_policy,
        dataset_view=dataset_view,
    )

    assert not technical_replicate_plan.requires_aggregation
    assert validated_design.analysis_sample_ids == request.design.sample_ids()


def test_differential_validator_production_collaborators_reuse_one_dataset_view(
    instrument_differential_dataset_view: _DatasetViewEvents,
) -> None:
    validated = DifferentialAnalysisValidator().run(
        build_duplicate_display_differential_request()
    )

    assert instrument_differential_dataset_view.views == [validated.dataset_view]
    assert (
        instrument_differential_dataset_view.dataframe_snapshots[
            "dataset.phospho internal snapshot"
        ]
        == 1
    )
    assert (
        instrument_differential_dataset_view.dataframe_snapshots[
            "dataset.site_metadata internal snapshot"
        ]
        == 1
    )
    assert (
        instrument_differential_dataset_view.optional_dataframe_snapshots == Counter()
    )

from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
    TechnicalReplicatePolicy,
)
from phospy.errors import WorkflowValidationError
from phospy.provenance.derived_quantitative import (
    DerivedQuantitativeDataProvenance,
)
from phospy.provenance.hashing import fingerprint_table
from phospy.science.datasets.derived_quantitative import (
    DerivedAnalysisReadyPhosphoDataset,
)
from phospy.science.datasets.models import DatasetPreprocessingReport
from phospy.science.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialComputationExecutor,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as ComputationRequest,
)
from phospy.validation.workflows.differential import ExperimentalDesignContractValidator
from phospy.workflows.differential.executor import DifferentialAnalysisExecutor
from phospy.workflows.differential.interpreter import DifferentialAnalysisInterpreter
from phospy.workflows.differential.replicates import (
    TechnicalReplicateAggregationPlanner,
    TechnicalReplicateAggregator,
    TechnicalReplicateResolver,
)
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.processing_state import (
    imputed_processing_state as valid_imputed_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_GENES = ["MAPK14", "GSK3B"]
_SITES = ["Y182", "S9"]
_DISPLAY_IDS = ["MAPK14;Y182;", "GSK3B;S9;"]


def _site_index() -> pd.Index:
    return protein_site_key_index(protein_identifiers=_GENES, sites=_SITES)


def _dataset_with_technical_replicates() -> AnalysisReadyPhosphoDataset:
    site_index = _site_index()
    phospho = pd.DataFrame(
        {
            "A1_T1": [1.0, 10.0],
            "A1_T2": [3.0, 8.0],
            "A2_T1": [2.0, 2.0],
            "A2_T2": [4.0, 4.0],
            "B1_T1": [5.0, 20.0],
            "B1_T2": [7.0, 18.0],
            "B2_T1": [6.0, 6.0],
            "B2_T2": [8.0, 8.0],
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": _DISPLAY_IDS,
            **site_key_context_columns(site_index),
            "gene_symbol": _GENES,
            "site": _SITES,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in _SITES
            ],
            "protein_id": _GENES,
        },
        index=phospho.index.copy(),
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _repeated_design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A1_T1",
                condition="A",
                biological_replicate_id="A1",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="A1_T2",
                condition="A",
                biological_replicate_id="A1",
                technical_replicate_id="T2",
            ),
            SampleDesignRecord(
                sample_id="A2_T1",
                condition="A",
                biological_replicate_id="A2",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="A2_T2",
                condition="A",
                biological_replicate_id="A2",
                technical_replicate_id="T2",
            ),
            SampleDesignRecord(
                sample_id="B1_T1",
                condition="B",
                biological_replicate_id="B1",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="B1_T2",
                condition="B",
                biological_replicate_id="B1",
                technical_replicate_id="T2",
            ),
            SampleDesignRecord(
                sample_id="B2_T1",
                condition="B",
                biological_replicate_id="B2",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="B2_T2",
                condition="B",
                biological_replicate_id="B2",
                technical_replicate_id="T2",
            ),
        )
    )


def _dataset_with_technical_replicates_and_total() -> AnalysisReadyPhosphoDataset:
    phospho_only = _dataset_with_technical_replicates()
    total = pd.DataFrame(
        {
            "A1_T1": [1.0, 2.0],
            "A1_T2": [3.0, 4.0],
            "A2_T1": [2.0, 1.0],
            "A2_T2": [4.0, 3.0],
            "B1_T1": [5.0, 6.0],
            "B1_T2": [7.0, 8.0],
            "B2_T1": [6.0, 5.0],
            "B2_T2": [8.0, 7.0],
        },
        index=pd.Index(["MAPK14", "GSK3B"], name="protein_id"),
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho_only.phospho,
        site_metadata=phospho_only.site_metadata,
        total=total,
        organism=phospho_only.organism,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=True
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=True),
    )


def _independent_design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A1_T1",
                condition="A",
                biological_replicate_id="A1",
            ),
            SampleDesignRecord(
                sample_id="A1_T2",
                condition="A",
                biological_replicate_id="A2",
            ),
            SampleDesignRecord(
                sample_id="B1_T1",
                condition="B",
                biological_replicate_id="B1",
            ),
            SampleDesignRecord(
                sample_id="B1_T2",
                condition="B",
                biological_replicate_id="B2",
            ),
        )
    )


def _contrasts() -> tuple[Contrast, ...]:
    return (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )


def _request(
    *,
    dataset: AnalysisReadyPhosphoDataset | None = None,
    design: ExperimentalDesign | None = None,
    technical_replicate_policy: TechnicalReplicatePolicy = (
        TechnicalReplicatePolicy.REJECT
    ),
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=_dataset_with_technical_replicates() if dataset is None else dataset,
        design=_repeated_design() if design is None else design,
        contrasts=_contrasts(),
        config=DifferentialAnalysisConfig(
            technical_replicate_policy=technical_replicate_policy
        ),
    )


def _aggregate(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    design: ExperimentalDesign,
    policy: TechnicalReplicatePolicy,
):
    plan = TechnicalReplicateAggregationPlanner().run(
        dataset=dataset,
        design=design,
        technical_replicate_policy=policy,
    )
    return TechnicalReplicateAggregator().run(
        dataset=dataset,
        design=design,
        aggregation_plan=plan,
    )


def _fingerprint_by_name(fingerprints, name: str):
    for fingerprint in fingerprints:
        if fingerprint.name == name:
            return fingerprint
    raise AssertionError(f"missing fingerprint: {name}")


def test_independent_biological_replicates_pass_unchanged() -> None:
    dataset = _dataset_with_technical_replicates()
    request = DifferentialAnalysisRequest(
        dataset=dataset,
        design=_independent_design(),
        contrasts=_contrasts(),
        config=DifferentialAnalysisConfig(allow_design_subset=True),
    )
    validated = DifferentialAnalysisValidator().run(request)
    assert validated.analysis_sample_ids == ("A1_T1", "A1_T2", "B1_T1", "B1_T2")
    assert validated.workflow_provenance is None
    assert validated.technical_replicate_aggregation_plan is not None
    assert not validated.technical_replicate_aggregation_plan.requires_aggregation


def test_repeated_biological_replicates_fail_with_default_reject_policy() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="Technical replicates require explicit aggregation",
    ):
        DifferentialAnalysisValidator().run(_request())


def test_repeated_biological_replicates_mean_policy_produces_explicit_plan() -> None:
    validated = DifferentialAnalysisValidator().run(
        _request(technical_replicate_policy=TechnicalReplicatePolicy.MEAN)
    )
    assert validated.dataset.phospho.columns.tolist() == [
        "A1_T1",
        "A1_T2",
        "A2_T1",
        "A2_T2",
        "B1_T1",
        "B1_T2",
        "B2_T1",
        "B2_T2",
    ]
    plan = validated.technical_replicate_aggregation_plan
    assert plan is not None
    assert plan.requires_aggregation
    assert plan.technical_replicate_policy is TechnicalReplicatePolicy.MEAN
    assert [group.output_sample_id for group in plan.groups] == ["A1", "A2", "B1", "B2"]


def test_repeated_biological_replicates_median_policy_produces_explicit_plan() -> None:
    validated = DifferentialAnalysisValidator().run(
        _request(technical_replicate_policy=TechnicalReplicatePolicy.MEDIAN)
    )
    plan = validated.technical_replicate_aggregation_plan
    assert plan is not None
    assert plan.requires_aggregation
    assert plan.technical_replicate_policy is TechnicalReplicatePolicy.MEDIAN


def test_aggregation_groups_by_condition_plus_biological_replicate_id() -> None:
    site_index = _site_index()
    phospho = pd.DataFrame(
        {
            "A_R1_T1": [1.0, 2.0],
            "A_R1_T2": [3.0, 4.0],
            "B_R1_T1": [10.0, 20.0],
            "B_R1_T2": [30.0, 40.0],
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": _DISPLAY_IDS,
            **site_key_context_columns(site_index),
            "gene_symbol": _GENES,
            "site": _SITES,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in _SITES
            ],
            "protein_id": _GENES,
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_R1_T1",
                condition="A",
                biological_replicate_id="R1",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="A_R1_T2",
                condition="A",
                biological_replicate_id="R1",
                technical_replicate_id="T2",
            ),
            SampleDesignRecord(
                sample_id="B_R1_T1",
                condition="B",
                biological_replicate_id="R1",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="B_R1_T2",
                condition="B",
                biological_replicate_id="R1",
                technical_replicate_id="T2",
            ),
        )
    )
    validated = DifferentialAnalysisValidator().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=design,
            contrasts=_contrasts(),
            config=DifferentialAnalysisConfig(
                technical_replicate_policy=TechnicalReplicatePolicy.MEAN,
                minimum_condition_replicates=1,
            ),
        )
    )
    plan = validated.technical_replicate_aggregation_plan
    assert plan is not None
    assert [group.output_sample_id for group in plan.groups] == ["A__R1", "B__R1"]


def test_validator_does_not_mutate_dataset_or_design_for_aggregation_policy() -> None:
    dataset = _dataset_with_technical_replicates()
    before = dataset.phospho.copy(deep=True)
    before_design_sample_ids = tuple(
        record.sample_id for record in _repeated_design().samples
    )
    request = _request(
        dataset=dataset,
        technical_replicate_policy=TechnicalReplicatePolicy.MEAN,
    )
    DifferentialAnalysisValidator().run(request)
    pdt.assert_frame_equal(before, dataset.phospho)
    assert (
        tuple(record.sample_id for record in request.design.samples)
        == before_design_sample_ids
    )


def test_validator_keeps_phospho_and_total_frames_immutable() -> None:
    dataset = _dataset_with_technical_replicates_and_total()
    phospho_before = dataset.phospho.copy(deep=True)
    assert dataset.total is not None
    total_before = dataset.total.copy(deep=True)
    request = _request(
        dataset=dataset,
        technical_replicate_policy=TechnicalReplicatePolicy.MEAN,
    )
    validated = DifferentialAnalysisValidator().run(request)
    assert validated.dataset is dataset
    pdt.assert_frame_equal(phospho_before, dataset.phospho)
    assert total_before is not None
    assert dataset.total is not None
    pdt.assert_frame_equal(total_before, dataset.total)


def test_executor_receives_aggregated_matrix() -> None:
    observed_columns: list[str] = []

    class _ComputationExecutorSpy:
        def run(self, request: ComputationRequest):
            observed_columns.extend(request.matrix.columns.astype(str).tolist())
            return DifferentialComputationExecutor().run(request)

    workflow = DifferentialAnalysisWorkflow(
        executor=DifferentialAnalysisExecutor(
            computation_executor=_ComputationExecutorSpy()  # type: ignore[arg-type]
        )
    )
    workflow.run(_request(technical_replicate_policy=TechnicalReplicatePolicy.MEAN))
    assert observed_columns == ["A1", "A2", "B1", "B2"]


def test_phospho_only_aggregation_marks_total_matrix_as_not_aggregated() -> None:
    result = DifferentialAnalysisWorkflow().run(
        _request(technical_replicate_policy=TechnicalReplicatePolicy.MEAN)
    )
    assert result.workflow_provenance is not None
    matrices = result.workflow_provenance["matrices_aggregated"]
    assert matrices == {"phospho": True, "total_protein": False}
    assert result.workflow_provenance["both_phospho_and_total_aggregated"] is False


def test_phospho_and_total_matrices_are_both_aggregated_when_total_present() -> None:
    dataset = _dataset_with_technical_replicates_and_total()
    plan = TechnicalReplicateAggregationPlanner().run(
        dataset=dataset,
        design=_repeated_design(),
        technical_replicate_policy=TechnicalReplicatePolicy.MEAN,
    )
    resolved = TechnicalReplicateAggregator().run(
        dataset=dataset,
        design=_repeated_design(),
        aggregation_plan=plan,
    )
    assert resolved.dataset.phospho.columns.tolist() == ["A1", "A2", "B1", "B2"]
    assert resolved.dataset.total is not None
    assert resolved.dataset.total.columns.tolist() == ["A1", "A2", "B1", "B2"]
    assert resolved.workflow_provenance is not None
    assert resolved.workflow_provenance["matrices_aggregated"] == {
        "phospho": True,
        "total_protein": True,
    }
    assert resolved.workflow_provenance["both_phospho_and_total_aggregated"] is True


def test_aggregation_returns_derived_dataset_with_fresh_provenance() -> None:
    base_dataset = _dataset_with_technical_replicates()
    source_report = DatasetPreprocessingReport.from_rows()
    dataset = AnalysisReadyPhosphoDataset(
        phospho=base_dataset.phospho,
        site_metadata=base_dataset.site_metadata,
        sample_metadata=base_dataset.sample_metadata,
        total=base_dataset.total,
        comparisons=base_dataset.comparisons,
        organism=base_dataset.organism,
        intensity_scale_state=base_dataset.intensity_scale_state,
        processing_state=base_dataset.processing_state,
        preprocessing_report=source_report,
        provenance=base_dataset.provenance,
    )

    resolved = _aggregate(
        dataset=dataset,
        design=_repeated_design(),
        policy=TechnicalReplicatePolicy.MEAN,
    )

    assert isinstance(resolved.dataset, DerivedAnalysisReadyPhosphoDataset)
    assert resolved.dataset.provenance is not dataset.provenance
    assert resolved.dataset.provenance is not None
    assert (
        resolved.dataset.provenance.workflow_name == "technical_replicate_aggregation"
    )
    assert resolved.dataset.preprocessing_report is None
    assert dataset.preprocessing_report is source_report

    lineage = resolved.dataset.derived_lineage
    source_phospho = fingerprint_table(dataset.phospho, name="dataset.phospho")
    derived_phospho = fingerprint_table(
        resolved.dataset.phospho,
        name="dataset.phospho",
    )
    parent_recorded = _fingerprint_by_name(
        lineage.parent_dataset_fingerprints,
        "dataset.phospho",
    )
    derived_recorded = _fingerprint_by_name(
        lineage.derived_dataset_fingerprints,
        "dataset.phospho",
    )
    assert parent_recorded.exact_hash_value == source_phospho.exact_hash_value
    assert derived_recorded.exact_hash_value == derived_phospho.exact_hash_value
    assert derived_recorded.rows == resolved.dataset.phospho.shape[0]
    assert derived_recorded.columns == resolved.dataset.phospho.shape[1]
    assert parent_recorded.exact_hash_value != derived_recorded.exact_hash_value
    assert lineage.sample_count == resolved.dataset.phospho.shape[1]
    assert lineage.sample_groups()[0] == ("A1", ("A1_T1", "A1_T2"))


def test_mean_and_median_derived_provenance_differ() -> None:
    site_index = _site_index()
    phospho = pd.DataFrame(
        {
            "A1_T1": [1.0, 10.0],
            "A1_T2": [2.0, 20.0],
            "A1_T3": [100.0, 30.0],
            "B1_T1": [4.0, 40.0],
            "B1_T2": [5.0, 50.0],
            "B1_T3": [6.0, 500.0],
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": _DISPLAY_IDS,
            **site_key_context_columns(site_index),
            "gene_symbol": _GENES,
            "site": _SITES,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in _SITES
            ],
            "protein_id": _GENES,
        },
        index=site_index.copy(),
    )
    dataset = AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )
    design = ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=sample_id[0],
                biological_replicate_id=sample_id[:2],
                technical_replicate_id=sample_id[-2:],
            )
            for sample_id in phospho.columns
        )
    )

    mean = _aggregate(
        dataset=dataset,
        design=design,
        policy=TechnicalReplicatePolicy.MEAN,
    )
    median = _aggregate(
        dataset=dataset,
        design=design,
        policy=TechnicalReplicatePolicy.MEDIAN,
    )

    assert isinstance(mean.dataset, DerivedAnalysisReadyPhosphoDataset)
    assert isinstance(median.dataset, DerivedAnalysisReadyPhosphoDataset)
    assert mean.dataset.derived_lineage.aggregation_method == "mean"
    assert median.dataset.derived_lineage.aggregation_method == "median"
    assert (
        mean.dataset.derived_lineage.lineage_hash_value
        != median.dataset.derived_lineage.lineage_hash_value
    )
    assert (
        _fingerprint_by_name(
            mean.dataset.derived_lineage.derived_dataset_fingerprints,
            "dataset.phospho",
        ).exact_hash_value
        != _fingerprint_by_name(
            median.dataset.derived_lineage.derived_dataset_fingerprints,
            "dataset.phospho",
        ).exact_hash_value
    )


def test_total_protein_and_imputation_mask_lineage_are_recorded() -> None:
    base = _dataset_with_technical_replicates_and_total()
    mask = pd.DataFrame(
        True,
        index=base.phospho.index.copy(),
        columns=base.phospho.columns.copy(),
    )
    mask.loc[base.phospho.index[0], "A1_T2"] = False
    mask.loc[base.phospho.index[1], "B2_T1"] = False
    dataset = AnalysisReadyPhosphoDataset(
        phospho=base.phospho,
        site_metadata=base.site_metadata,
        sample_metadata=base.sample_metadata,
        total=base.total,
        organism=base.organism,
        intensity_scale_state=base.intensity_scale_state,
        processing_state=valid_imputed_processing_state(base.processing_state),
        imputation_observation_mask=mask,
        provenance=base.provenance,
    )

    resolved = _aggregate(
        dataset=dataset,
        design=_repeated_design(),
        policy=TechnicalReplicatePolicy.MEAN,
    )

    assert isinstance(resolved.dataset, DerivedAnalysisReadyPhosphoDataset)
    observed_mask = resolved.dataset.imputation_observed_mask_dataframe()
    assert observed_mask is not None
    assert bool(observed_mask.loc[base.phospho.index[0], "A1"]) is False
    assert bool(observed_mask.loc[base.phospho.index[1], "B2"]) is False
    lineage = resolved.dataset.derived_lineage
    assert lineage.matrices_transformed["total_protein"] is True
    assert lineage.matrices_transformed["imputation_observation_mask"] is True
    assert (
        _fingerprint_by_name(
            lineage.derived_dataset_fingerprints, "dataset.total"
        ).columns
        == 4
    )
    assert (
        _fingerprint_by_name(
            lineage.derived_dataset_fingerprints,
            "dataset.imputation_observation_mask",
        ).columns
        == 4
    )


def test_serialized_lineage_round_trips_and_verifies_actual_matrix() -> None:
    dataset = _dataset_with_technical_replicates()
    resolved = _aggregate(
        dataset=dataset,
        design=_repeated_design(),
        policy=TechnicalReplicatePolicy.MEAN,
    )

    assert isinstance(resolved.dataset, DerivedAnalysisReadyPhosphoDataset)
    lineage = resolved.dataset.derived_lineage
    payload = lineage.to_payload()
    restored = DerivedQuantitativeDataProvenance.from_payload(payload)
    assert restored.to_payload() == payload

    derived_phospho = fingerprint_table(
        resolved.dataset.phospho,
        name="dataset.phospho",
    )
    restored_phospho = _fingerprint_by_name(
        restored.derived_dataset_fingerprints,
        "dataset.phospho",
    )
    assert restored_phospho.exact_hash_value == derived_phospho.exact_hash_value
    replayed = pd.concat(
        [
            dataset.phospho.loc[:, list(input_ids)].mean(axis=1).rename(output_id)
            for output_id, input_ids in restored.sample_groups()
        ],
        axis=1,
    )
    replayed.columns = pd.Index([item[0] for item in restored.sample_groups()])
    pdt.assert_frame_equal(replayed, resolved.dataset.phospho)


def test_interpreter_consumes_derived_dataset_after_aggregation() -> None:
    validated = DifferentialAnalysisValidator().run(
        _request(technical_replicate_policy=TechnicalReplicatePolicy.MEAN)
    )
    seen_dataset_types: list[type[AnalysisReadyPhosphoDataset]] = []
    real_validator = ExperimentalDesignContractValidator()

    class _DesignValidatorSpy:
        def run(self, **kwargs):
            seen_dataset_types.append(type(kwargs["dataset"]))
            return real_validator.run(**kwargs)

    DifferentialAnalysisInterpreter(
        design_validator=_DesignValidatorSpy(),  # type: ignore[arg-type]
    ).run(validated)

    assert DerivedAnalysisReadyPhosphoDataset in seen_dataset_types


def test_technical_replicate_resolver_warns_and_preserves_wrapper_behaviour() -> None:
    with pytest.warns(
        DeprecationWarning,
        match="TechnicalReplicateResolver is deprecated",
    ) as construction_warnings:
        resolver = TechnicalReplicateResolver()
    assert "TechnicalReplicateAggregationPlanner" in str(
        construction_warnings[0].message
    )
    assert "TechnicalReplicateAggregator" in str(construction_warnings[0].message)

    dataset = _dataset_with_technical_replicates()
    design = _repeated_design()
    with pytest.warns(
        DeprecationWarning,
        match="TechnicalReplicateResolver is deprecated",
    ) as run_warnings:
        resolved = resolver.run(
            dataset=dataset,
            design=design,
            technical_replicate_policy=TechnicalReplicatePolicy.MEAN,
        )

    assert "TechnicalReplicateAggregationPlanner" in str(run_warnings[0].message)
    assert "TechnicalReplicateAggregator" in str(run_warnings[0].message)
    assert resolved.dataset.phospho.columns.tolist() == ["A1", "A2", "B1", "B2"]
    assert [record.sample_id for record in resolved.design.samples] == [
        "A1",
        "A2",
        "B1",
        "B2",
    ]
    assert resolved.workflow_provenance is not None
    assert resolved.workflow_provenance["aggregation_method"] == "mean"


def test_provenance_records_technical_replicate_lineage() -> None:
    result = DifferentialAnalysisWorkflow().run(
        _request(technical_replicate_policy=TechnicalReplicatePolicy.MEAN)
    )
    assert result.workflow_provenance is not None
    assert result.workflow_provenance["technical_replicate_policy"] == "mean"
    assert result.workflow_provenance["aggregation_policy"] == "mean"
    assert result.workflow_provenance["aggregation_method"] == "mean"
    assert result.workflow_provenance["grouped_samples"] == ["A1", "A2", "B1", "B2"]
    assert result.workflow_provenance["source_samples"] == [
        "A1_T1",
        "A1_T2",
        "A2_T1",
        "A2_T2",
        "B1_T1",
        "B1_T2",
        "B2_T1",
        "B2_T2",
    ]
    groups = result.workflow_provenance["groups"]
    assert isinstance(groups, list)
    a1_group = next(
        group
        for group in groups
        if group["condition"] == "A" and group["biological_replicate_id"] == "A1"
    )
    assert a1_group["output_sample_id"] == "A1"
    assert a1_group["input_sample_ids"] == ["A1_T1", "A1_T2"]
    assert a1_group["source_sample_ids"] == ["A1_T1", "A1_T2"]
    assert a1_group["technical_replicate_ids"] == ["T1", "T2"]
    assert a1_group["n_technical_replicates"] == 2
    assert a1_group["aggregation_method"] == "mean"
    assert result.policy_provenance is not None
    assert result.policy_provenance.replicates.technical_replicate_policy == "mean"
    structured_groups = result.policy_provenance.replicates.technical_replicate_groups
    assert len(structured_groups) == 4
    a1_structured = next(
        group
        for group in structured_groups
        if group.condition == "A" and group.biological_replicate_id == "A1"
    )
    assert a1_structured.output_sample_id == "A1"
    assert a1_structured.input_sample_ids == ("A1_T1", "A1_T2")
    assert a1_structured.technical_replicate_ids == ("T1", "T2")
    assert a1_structured.n_technical_replicates == 2


def test_invalid_technical_replicate_policy_fails() -> None:
    request = DifferentialAnalysisRequest(
        dataset=_dataset_with_technical_replicates(),
        design=_repeated_design(),
        contrasts=_contrasts(),
        config=DifferentialAnalysisConfig(
            technical_replicate_policy="invalid"  # type: ignore[arg-type]
        ),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="technical_replicate_policy must be TechnicalReplicatePolicy",
    ):
        DifferentialAnalysisValidator().run(request)

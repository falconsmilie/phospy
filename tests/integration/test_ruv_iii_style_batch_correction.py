from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd
import pytest

from phospy.advanced import (
    ControlSiteSet,
    ControlSiteSourceMetadata,
    CorrectionMissingnessPolicy,
    ObservationMask,
    OriginallyMissingCellTracking,
    SpsDiscoveryConfig,
    SpsDiscoveryRequest,
    SpsDiscoveryWorkflow,
    SpsReferenceDataset,
    SpsRuvBatchCorrectionConfig,
    TemporaryImputationMethod,
    TemporaryImputationPolicy,
)
from phospy.api import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    Organism,
)
from phospy.errors import PhosPyInputError
from phospy.provenance.serialization import (
    batch_correction_provenance_from_payload,
    batch_correction_provenance_to_payload,
)
from phospy.science.batch_correction import RuvIIIStyleExecutor
from phospy.validation.workflows.batch_correction.control_site_workflow import (
    BatchCorrectionWorkflowControlSiteValidator,
)
from phospy.validation.workflows.batch_correction.design import (
    BatchCorrectionWorkflowDesignValidator,
    BatchCorrectionWorkflowFactorFeasibilityValidator,
)
from phospy.validation.workflows.batch_correction.missingness import (
    BatchCorrectionWorkflowMissingnessValidator,
)
from phospy.workflows.batch_correction import (
    BatchCorrectionPlanInterpreter,
    BatchCorrectionWorkflowRequest,
)

pytestmark = pytest.mark.integration

_SAMPLES = tuple(f"sample_{position}" for position in range(1, 7))


def test_public_ruv_iii_style_executes_with_active_replicates_and_provenance() -> None:
    phospho = _phospho()
    dataset = _build(phospho=phospho, replicates=("r1", "r1", "r2", "r2", "r3", "r3"))

    report = dataset.preprocessing_report.batch_correction
    assert report is not None
    assert report.method == "ruv_iii_style"
    assert not dataset.phospho.equals(phospho)
    correction = next(
        stage.batch_correction_provenance
        for stage in dataset.provenance.preprocessing_stages
        if stage.stage == "batch_correction"
    )
    assert correction is not None
    assert correction.requested_method == "ruv_iii_style"
    assert correction.selected_site_key_rows == _control_keys()
    assert correction.replicate_metadata is not None
    assert correction.replicate_metadata["used_for_numerical_factor_estimation"] is True
    assert correction.replicate_metadata["ruv_iii_semantics_enabled"] is True
    config = correction.resolved_parameters["config"]
    assert isinstance(config, Mapping)
    assert config["n_unwanted_factors"] == 1
    interpreter_plan = correction.resolved_parameters["interpreter_plan"]
    assert isinstance(interpreter_plan, Mapping)
    assert interpreter_plan["n_unwanted_factors"] == 1
    executor_provenance = correction.resolved_parameters["executor"]
    assert isinstance(executor_provenance, Mapping)
    assert executor_provenance["k"] == 1
    executor = correction.diagnostics["executor"]
    assert isinstance(executor, Mapping)
    assert executor["method"] == "ruv_iii_style"
    assert executor["algorithm_id"] == "replicate_aware_ruv_iii_v1"
    assert executor["requested_unwanted_factors"] == 1
    assert executor["estimated_unwanted_factors"] == 1
    assert executor["missingness_imputation_summary"] == {
        "strategy": "no_internal_completion_required",
        "temporary_completion_applied": False,
        "temporary_imputation_method": "none",
        "temporary_imputation_parameters": {},
        "governed_missing_cell_count": 0,
        "restored_missing_cell_count": 0,
        "upstream_imputed_input_cell_count": 0,
        "imputed_values_are_observed_evidence": False,
        "output_policy": "no governed missing positions",
    }

    restored = batch_correction_provenance_from_payload(
        batch_correction_provenance_to_payload(correction)
    )
    assert restored == correction


def test_ruv_iii_style_replicate_definition_changes_numerical_result() -> None:
    phospho = _phospho()
    conditions = ("a",) * len(_SAMPLES)
    paired = _build(
        phospho=phospho,
        replicates=("r1", "r1", "r2", "r2", "r3", "r3"),
        conditions=conditions,
    )
    alternative = _build(
        phospho=phospho,
        replicates=("r1", "r2", "r1", "r2", "r3", "r3"),
        conditions=conditions,
    )

    assert not paired.phospho.equals(alternative.phospho)


def test_discovery_control_set_is_consumed_directly_without_reference_matrices() -> (
    None
):
    controls = _discovered_controls()
    config = _config(replicate_column="replicate", control_site_set=controls)

    assert config.control_site_set is controls
    dataset = _build(
        phospho=_phospho(),
        replicates=("r1", "r1", "r2", "r2", "r3", "r3"),
        config=config,
    )

    correction = next(
        stage.batch_correction_provenance
        for stage in dataset.provenance.preprocessing_stages
        if stage.stage == "batch_correction"
    )
    assert correction is not None
    assert correction.control_site_source["control_site_set_source_type"] == (
        "sps_discovery"
    )
    assert "reference_datasets" not in correction.control_site_source
    assert "intensities" not in str(correction.resolved_parameters)


def test_ruv_iii_executor_temporarily_completes_and_restores_missing_positions() -> (
    None
):
    phospho = _resolved_phospho()
    missing_site = str(phospho.index[2])
    phospho.loc[missing_site, "sample_1"] = np.nan
    mask = ObservationMask(
        feature_ids=tuple(phospho.index.astype(str)),
        sample_ids=tuple(phospho.columns.astype(str)),
        originally_missing_cells=((missing_site, "sample_1"),),
    )
    policy = CorrectionMissingnessPolicy(
        temporary_imputation=TemporaryImputationPolicy(
            allowed=True,
            method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
            method_parameters=(("min_observed_values", 2),),
        ),
        originally_missing_cells_tracked_by=(
            OriginallyMissingCellTracking.OBSERVATION_MASK
        ),
        observation_mask=mask,
    )
    config = _config(replicate_column="replicate", missingness_policy=policy)
    request = BatchCorrectionWorkflowRequest(
        phospho=phospho,
        config=config.to_internal_request(),
        sample_metadata=_resolved_sample_metadata(),
        control_site_set=config.control_site_set,
        missingness_policy=policy,
        dataset_organism=Organism.RAT,
    )
    metadata = BatchCorrectionWorkflowDesignValidator().run(request=request)
    mapping = BatchCorrectionWorkflowControlSiteValidator().run(request=request)
    resolved_policy = BatchCorrectionWorkflowMissingnessValidator().run(request=request)
    BatchCorrectionWorkflowFactorFeasibilityValidator().run(
        request=request,
        dataset_metadata=metadata,
        control_site_mapping=mapping,
        missingness_policy=resolved_policy,
    )
    plan = BatchCorrectionPlanInterpreter().run(
        config=request.config,
        dataset_metadata=metadata,
        control_site_mapping=mapping,
        missingness_policy=resolved_policy,
    )

    result = RuvIIIStyleExecutor().run(phospho=phospho, plan=plan)

    assert pd.isna(result.corrected_matrix.loc[missing_site, "sample_1"])
    assert not bool(result.output_observation_mask.loc[missing_site, "sample_1"])
    assert result.corrected_preprocessing_output is None
    summary = result.diagnostics.missingness_imputation_summary
    assert summary["temporary_completion_applied"] is True
    assert summary["restored_missing_cell_count"] == 1
    assert summary["imputed_values_are_observed_evidence"] is False


def test_ruv_iii_style_requires_replicate_column_and_non_singleton_sets() -> None:
    with pytest.raises(PhosPyInputError, match="replicate_column is required"):
        _config(replicate_column=None)

    with pytest.raises(
        PhosPyInputError,
        match="RUV-III replicate sets must each contain at least two samples",
    ):
        _build(
            phospho=_phospho(),
            replicates=("r1", "r2", "r3", "r4", "r5", "r6"),
        )


def test_public_ruv_iii_style_rejects_missing_configured_replicate_metadata() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="sample_metadata is missing required column 'missing_replicate'",
    ):
        _build(
            phospho=_phospho(),
            replicates=("r1", "r1", "r2", "r2", "r3", "r3"),
            config=_config(replicate_column="missing_replicate"),
        )


def test_public_ruv_iii_style_rejects_invalid_and_infeasible_k() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="n_unwanted_factors must be greater than or equal to 1",
    ):
        _config(replicate_column="replicate", n_unwanted_factors=0)

    controls = ControlSiteSet.from_site_keys(
        _control_keys(count=5),
        source_metadata=ControlSiteSourceMetadata(
            organism="rat",
            identifier_namespace="site_key",
            source_name="manual RUV-III controls",
            source_version="v1",
            license="caller local use",
            redistribution="not redistributed",
        ),
    )
    with pytest.raises(
        PhosPyInputError,
        match=("requested k=4 exceeds replicate-residual degrees of freedom=3"),
    ):
        _build(
            phospho=_phospho(),
            replicates=("r1", "r1", "r2", "r2", "r3", "r3"),
            config=_config(
                replicate_column="replicate",
                control_site_set=controls,
                n_unwanted_factors=4,
            ),
        )


def test_public_ruv_iii_style_enforces_protected_condition_strata() -> None:
    aligned = _build(
        phospho=_phospho(),
        replicates=("r1", "r1", "r2", "r2", "r3", "r3"),
    )
    assert aligned.preprocessing_report.batch_correction is not None

    with pytest.raises(
        PhosPyInputError,
        match=(
            "RUV-III replicate sets must not cross protected condition strata; "
            "replicate set 'r1' contains conditions 'a', 'b'"
        ),
    ):
        _build(
            phospho=_phospho(),
            replicates=("r1", "r2", "r1", "r2", "r3", "r3"),
        )


def _build(
    *,
    phospho: pd.DataFrame,
    replicates: tuple[str, ...],
    conditions: tuple[str, ...] = ("a", "a", "b", "b", "c", "c"),
    config: SpsRuvBatchCorrectionConfig | None = None,
):
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho),
            sample_metadata=pd.DataFrame(
                {
                    "batch": ("run_1", "run_2", "run_1", "run_2", "run_1", "run_2"),
                    "condition": conditions,
                    "replicate": replicates,
                },
                index=phospho.columns.copy(),
            ),
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                batch_correction=(
                    _config(replicate_column="replicate") if config is None else config
                )
            ),
        )
    )


def _config(
    *,
    replicate_column: str | None,
    control_site_set: ControlSiteSet | None = None,
    missingness_policy: CorrectionMissingnessPolicy | None = None,
    n_unwanted_factors: int = 1,
) -> SpsRuvBatchCorrectionConfig:
    return SpsRuvBatchCorrectionConfig(
        control_site_set=(
            ControlSiteSet.from_site_keys(
                _control_keys(),
                source_metadata=ControlSiteSourceMetadata(
                    organism="rat",
                    identifier_namespace="site_key",
                    source_name="manual RUV-III controls",
                    source_version="v1",
                    license="caller local use",
                    redistribution="not redistributed",
                ),
            )
            if control_site_set is None
            else control_site_set
        ),
        batch_column="batch",
        condition_columns=("condition",),
        replicate_column=replicate_column,
        missingness_policy=(
            CorrectionMissingnessPolicy()
            if missingness_policy is None
            else missingness_policy
        ),
        n_unwanted_factors=n_unwanted_factors,
        method="ruv_iii_style",
    )


def _phospho() -> pd.DataFrame:
    unwanted = np.asarray((-2.0, 2.0, -1.0, 1.0, -3.0, 3.0))
    biological = np.asarray((-4.0, -4.0, 1.0, 1.0, 5.0, 5.0))
    values = np.vstack(
        (
            15.0 + unwanted,
            12.0 - (2.0 * unwanted),
            18.0 + biological + (0.5 * unwanted),
            12.0 + biological,
            np.repeat(7.0, len(_SAMPLES)),
        )
    )
    return pd.DataFrame(
        values,
        index=pd.Index(
            ("P0001;S1;", "P0002;S2;", "P0003;S3;", "P0004;S4;", "P0005;S5;"),
            name="site_id",
        ),
        columns=_SAMPLES,
    )


def _site_metadata(phospho: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": tuple(f"G{position}" for position in range(1, 6)),
            "protein_id": tuple(f"P{position:04d}" for position in range(1, 6)),
            "site": tuple(f"S{position}" for position in range(1, 6)),
            "site_sequence": tuple(("A" * 15) + "S" + ("A" * 15) for _ in range(5)),
            "localisation_confidence": (0.99,) * 5,
        },
        index=phospho.index.copy(),
    )


def _control_keys(*, count: int = 2) -> tuple[str, ...]:
    return tuple(
        "phospy:v1|organism=rat|protein_namespace=protein_id|"
        f"protein_identifier=P{position:04d}|residue=S|position={position}"
        for position in range(1, count + 1)
    )


def _resolved_phospho() -> pd.DataFrame:
    resolved = _phospho().copy(deep=True)
    resolved.index = pd.Index(
        (
            *_control_keys(),
            "phospy:v1|organism=rat|protein_namespace=protein_id|"
            "protein_identifier=P0003|residue=S|position=3",
            "phospy:v1|organism=rat|protein_namespace=protein_id|"
            "protein_identifier=P0004|residue=S|position=4",
            "phospy:v1|organism=rat|protein_namespace=protein_id|"
            "protein_identifier=P0005|residue=S|position=5",
        ),
        name="site_key",
    )
    return resolved


def _resolved_sample_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ("run_1", "run_2", "run_1", "run_2", "run_1", "run_2"),
            "condition": ("a", "a", "b", "b", "c", "c"),
            "replicate": ("r1", "r1", "r2", "r2", "r3", "r3"),
        },
        index=pd.Index(_SAMPLES),
    )


def _discovered_controls() -> ControlSiteSet:
    site_keys = pd.Index(
        (
            *_control_keys(),
            "phospy:v1|organism=rat|protein_namespace=protein_id|"
            "protein_identifier=P0003|residue=S|position=3",
        ),
        name="site_key",
    )
    condition_by_sample = {
        "a_1": "a",
        "a_2": "a",
        "b_1": "b",
        "b_2": "b",
    }
    references = tuple(
        SpsReferenceDataset.from_condition_relative_log2(
            dataset_id=dataset_id,
            intensities=pd.DataFrame(values, index=site_keys),
            condition_by_sample=condition_by_sample,
            log2_scale_established_by=f"{dataset_id} governed log2 transform",
            baseline_centering_established_by=(
                f"{dataset_id} governed condition-relative centering"
            ),
            source_name=f"governed-{dataset_id}",
            source_version="2026-09",
        )
        for dataset_id, values in (
            (
                "reference_a",
                {
                    "a_1": (0.0, 0.0, 0.0),
                    "a_2": (0.0, 0.0, 0.0),
                    "b_1": (0.1, 0.2, 2.0),
                    "b_2": (0.1, 0.2, 2.0),
                },
            ),
            (
                "reference_b",
                {
                    "a_1": (0.0, 0.0, 0.0),
                    "a_2": (0.0, 0.0, 0.0),
                    "b_1": (0.2, 0.3, 3.0),
                    "b_2": (0.2, 0.3, 3.0),
                },
            ),
        )
    )
    return (
        SpsDiscoveryWorkflow()
        .run(
            SpsDiscoveryRequest(
                reference_datasets=references,
                config=SpsDiscoveryConfig(top_n=2, minimum_shared_sites=2),
            )
        )
        .control_site_set
    )

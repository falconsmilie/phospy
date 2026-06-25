from __future__ import annotations

from collections.abc import Mapping
from typing import cast

import numpy as np
import pandas as pd
import pytest

from phospy.api import (
    AnalysisReadyDatasetBuilder,
    ControlSiteSet,
    DatasetBuildRequest,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    Organism,
    SpsRuvBatchCorrectionConfig,
)
from phospy.contracts.configs.preprocessing import (
    CorrectionMaskPolicy,
    CorrectionMissingnessPolicy,
    InternalBatchCorrectionControlSiteMode,
    InternalBatchCorrectionControlSiteSource,
    InternalBatchCorrectionImputationPolicy,
    InternalBatchCorrectionMethod,
    InternalBatchCorrectionMissingValuePolicy,
    InternalBatchCorrectionRequest,
    InternalBatchCorrectionStageOrder,
    OriginallyMissingCellTracking,
    TemporaryImputationMethod,
    TemporaryImputationPolicy,
)
from phospy.errors import PhosPyInputError
from phospy.provenance.models import BatchCorrectionProvenance
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteSourceMetadata,
)
from phospy.workflows.batch_correction import (
    BatchCorrectionWorkflow,
    BatchCorrectionWorkflowRequest,
)
from phospy.workflows.batch_correction.contracts import BatchCorrectionExecutorContract

_MAPK14_SITE_KEY = (
    "phospy:v1|organism=rat|protein_namespace=protein_id|"
    "protein_identifier=MAPK14|residue=Y|position=182"
)
_AKT1_SITE_KEY = (
    "phospy:v1|organism=rat|protein_namespace=protein_id|"
    "protein_identifier=AKT1|residue=T|position=308"
)
_GSK3B_SITE_KEY = (
    "phospy:v1|organism=rat|protein_namespace=protein_id|"
    "protein_identifier=GSK3B|residue=S|position=9"
)


def test_dataset_sps_ruv_preserves_upstream_imputation_observation_mask() -> None:
    phospho = _phospho()
    phospho.loc["AKT1;T308;", "sample_2"] = np.nan

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            sample_metadata=_sample_metadata(),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=_preprocessing_config(),
        )
    )

    observed_mask = built.imputation_observed_mask_dataframe()
    assert observed_mask is not None
    assert bool(observed_mask.loc[_AKT1_SITE_KEY, "sample_2"]) is False
    assert np.isfinite(built.phospho.to_numpy(dtype="float64")).all()

    provenance = _batch_correction_provenance(built)
    executor_diagnostics = cast(
        Mapping[str, object],
        provenance.diagnostics["executor"],
    )
    assert provenance.phospy_version
    assert provenance.phospy_version != "unknown"
    assert provenance.python_version
    assert provenance.python_version != "unknown"
    assert {"numpy", "pandas"}.issubset(set(provenance.dependency_versions))
    missingness_summary = cast(
        Mapping[str, object],
        executor_diagnostics["missingness_imputation_summary"],
    )
    assert executor_diagnostics["originally_missing_cell_count"] == 1
    assert missingness_summary["upstream_imputed_input_cell_count"] == 1
    provenance_mask = cast(
        Mapping[str, object],
        provenance.missing_value_policy["observation_mask"],
    )
    assert provenance_mask["originally_missing_cell_count"] == 1
    assert "batch_correction.workflow.upstream_observation_mask" in {
        fingerprint.name for fingerprint in provenance.observation_masks
    }
    observation_summary = built.imputation_observation_summary_dataframe(
        feature_ids=[_AKT1_SITE_KEY],
        sample_ids=["sample_2"],
    )
    assert observation_summary is not None
    assert observation_summary["observed_cell_count"].tolist() == [0]


def test_dataset_sps_ruv_without_upstream_missing_values_keeps_all_true_mask() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(_phospho().index),
            sample_metadata=_sample_metadata(),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=_preprocessing_config(),
        )
    )

    observed_mask = built.imputation_observed_mask_dataframe()
    assert observed_mask is not None
    assert bool(observed_mask.to_numpy(dtype=bool).all())

    provenance = _batch_correction_provenance(built)
    executor_diagnostics = cast(
        Mapping[str, object],
        provenance.diagnostics["executor"],
    )
    assert executor_diagnostics["originally_missing_cell_count"] == 0
    assert not any("originally missing" in warning for warning in provenance.warnings)


def test_dataset_sps_ruv_records_replicate_metadata_for_provenance_only() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(_phospho().index),
            sample_metadata=_sample_metadata(),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=_preprocessing_config(),
        )
    )

    provenance = _batch_correction_provenance(built)
    replicate_metadata = cast(Mapping[str, object], provenance.replicate_metadata)

    assert replicate_metadata == {
        "replicate_by_sample": {
            "sample_1": "r1",
            "sample_2": "r1",
            "sample_3": "r2",
            "sample_4": "r2",
        },
        "replicate_labels": ["r1", "r1", "r2", "r2"],
    }
    assert provenance.requested_method == "sps_ruv_style"
    config_payload = cast(
        Mapping[str, object], provenance.resolved_parameters["config"]
    )
    assert config_payload["method"] == "sps_ruv_style"
    assert config_payload["replicate_column"] == "replicate"


def test_workflow_rejects_misaligned_upstream_mask_before_executor() -> None:
    executor = _SpyExecutor()
    request = _workflow_request(
        upstream_observation_mask=_upstream_mask().loc[
            :, ["sample_2", "sample_1", "sample_3", "sample_4"]
        ],
        missingness_policy=_missingness_policy(),
    )

    with pytest.raises(PhosPyInputError, match="mask alignment"):
        BatchCorrectionWorkflow(
            executor=cast(BatchCorrectionExecutorContract, executor)
        ).run(request)

    assert executor.call_count == 0


def test_workflow_rejects_unsupported_upstream_imputed_policy_before_executor() -> None:
    executor = _SpyExecutor()
    request = _workflow_request(
        upstream_observation_mask=_upstream_mask(),
        missingness_policy=CorrectionMissingnessPolicy(),
    )

    with pytest.raises(PhosPyInputError, match="upstream-imputed.*observed evidence"):
        BatchCorrectionWorkflow(
            executor=cast(BatchCorrectionExecutorContract, executor)
        ).run(request)

    assert executor.call_count == 0


class _SpyExecutor:
    call_count = 0

    def run(self, **_: object) -> object:
        self.call_count += 1
        raise AssertionError("executor should not be called")


def _preprocessing_config() -> DatasetPreprocessingConfig:
    return DatasetPreprocessingConfig(
        missing_data=DatasetMissingDataConfig(
            policy="impute_row_median",
            min_observed_values=2,
        ),
        batch_correction=SpsRuvBatchCorrectionConfig(
            control_site_set=ControlSiteSet.from_site_keys(
                (_MAPK14_SITE_KEY, _GSK3B_SITE_KEY),
                source_metadata=_control_source_metadata(),
            ),
            batch_column="batch",
            condition_columns=("condition",),
            replicate_column="replicate",
            missingness_policy=_missingness_policy(),
            n_unwanted_factors=1,
            diagnostics_enabled=True,
            provenance_enabled=True,
        ),
    )


def _missingness_policy() -> CorrectionMissingnessPolicy:
    return CorrectionMissingnessPolicy(
        temporary_imputation=TemporaryImputationPolicy(
            allowed=True,
            method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
            method_parameters={"min_observed_values": 2},  # type: ignore[arg-type]
        ),
        originally_missing_cells_tracked_by=(
            OriginallyMissingCellTracking.OBSERVATION_MASK
        ),
        correction_mask_policy=CorrectionMaskPolicy(),
    )


def _workflow_request(
    *,
    upstream_observation_mask: pd.DataFrame,
    missingness_policy: CorrectionMissingnessPolicy,
) -> BatchCorrectionWorkflowRequest:
    return BatchCorrectionWorkflowRequest(
        phospho=_phospho(),
        config=_internal_config(),
        sample_metadata=_sample_metadata(),
        control_site_set=ControlSiteSet.from_site_keys(
            ("MAPK14;Y182;", "GSK3B;S9;"),
            source_metadata=_control_source_metadata(),
        ),
        missingness_policy=missingness_policy,
        upstream_observation_mask=upstream_observation_mask,
    )


def _control_source_metadata() -> ControlSiteSourceMetadata:
    return ControlSiteSourceMetadata(
        organism="rat",
        identifier_namespace="site_key",
        source_name="manual-curated-controls",
        source_version="manual-v1",
        license="caller local use",
        redistribution="not redistributed",
    )


def _internal_config() -> InternalBatchCorrectionRequest:
    return InternalBatchCorrectionRequest(
        method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
        batch_column="batch",
        condition_columns=("condition",),
        replicate_column="replicate",
        control_site_source=InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED,
        control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
        missing_value_policy=(
            InternalBatchCorrectionMissingValuePolicy.ALLOW_TEMPORARY_IMPUTATION
        ),
        imputation_policy=InternalBatchCorrectionImputationPolicy.ROW_MEDIAN_TEMPORARY,
        n_unwanted_factors=1,
        stage_order=InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM,
        diagnostics_enabled=True,
    )


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [10.0, 5.0, 20.0],
            "sample_2": [10.0, 9.0, 20.0],
            "sample_3": [14.0, 8.0, 28.0],
            "sample_4": [14.0, 12.0, 28.0],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"],
            name="site_id",
        ),
    )


def _upstream_mask() -> pd.DataFrame:
    mask = pd.DataFrame(
        True,
        index=_phospho().index.copy(),
        columns=_phospho().columns.copy(),
    )
    mask.loc["AKT1;T308;", "sample_2"] = False
    return mask


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
            "protein_id": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "site_sequence": [
                ("A" * 15) + "Y" + ("A" * 15),
                ("A" * 15) + "T" + ("A" * 15),
                ("A" * 15) + "S" + ("A" * 15),
            ],
            "localisation_confidence": [0.95, 0.92, 0.9],
        },
        index=index.copy(),
    )


def _sample_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ("run_1", "run_1", "run_2", "run_2"),
            "condition": ("control", "treated", "control", "treated"),
            "replicate": ("r1", "r1", "r2", "r2"),
        },
        index=("sample_1", "sample_2", "sample_3", "sample_4"),
    )


def _batch_correction_provenance(
    built: AnalysisReadyPhosphoDataset,
) -> BatchCorrectionProvenance:
    assert built.provenance is not None
    for item in built.provenance.preprocessing_stages:
        if item.stage != "batch_correction":
            continue
        provenance = item.batch_correction_provenance
        assert provenance is not None
        return provenance
    raise AssertionError("batch_correction provenance was not recorded")

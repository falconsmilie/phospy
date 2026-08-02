from __future__ import annotations

import ast
import importlib
import inspect
import os
import subprocess
import sys
import warnings
from collections.abc import Callable, Iterator, Mapping
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
    SampleDesignRecord,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    DifferentialAnalysisResult,
    KinasePredictionResult,
    KinaseWorkflowResult,
    PhosphositeImportResult,
    SignalomeWorkflowResult,
)
from phospy.errors.validation import DatasetValidationError, PhosPyValidationError
from phospy.frames.ownership import (
    _borrow_dataframe,
    _borrow_series,
    export_dataframe,
    export_optional_dataframe,
    export_optional_series,
    export_series,
    own_dataframe,
    own_optional_dataframe,
    own_optional_series,
    own_series,
)
from phospy.provenance.hashing import fingerprint_table
from phospy.science.activities.models import (
    KinaseActivityInputs,
    KinaseActivityResult,
    PredMatOverlapSummary,
)
from phospy.science.activities.semantics import ActivityInputMatrix
from phospy.science.datasets.builders.executor import DatasetBuildExecutor
from phospy.science.datasets.builders.interpreter import DatasetBuildRequestInterpreter
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import DatasetPreprocessingReport
from phospy.science.datasets.preprocessing.report_schema import (
    ComparisonGroupStatsRow,
    ComparisonPairStatsRow,
    DuplicateSiteResolutionRow,
    MetadataConflictRow,
    PreprocessingOperationRow,
    PreprocessingRowAuditRow,
    PreprocessingRowCountRow,
)
from phospy.science.differential.models.diagnostics import (
    EmpiricalBayesPriorDiagnostics,
)
from phospy.science.prediction.models import KinaseScoringResult
from phospy.science.signalomes.constants import (
    CORRELATION_COLUMN,
    DISPLAY_ID_COLUMN,
    GENE_SYMBOL_COLUMN,
    ISOFORM_ID_COLUMN,
    MODULE_ID_COLUMN,
    MODULE_TOP_KINASE_CANDIDATES_COLUMN,
    MODULE_TOP_KINASE_COLUMN,
    MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN,
    MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN,
    MODULE_TOP_KINASE_TIE_COUNT_COLUMN,
    PROTEIN_ACCESSION_COLUMN,
    PROTEIN_COLUMN,
    SITE_COLUMN,
    SITE_KEY_COLUMN,
    SOURCE_KINASE_COLUMN,
    TARGET_KINASE_COLUMN,
    TOP_KINASE_CANDIDATES_COLUMN,
    TOP_KINASE_COLUMN,
    TOP_KINASE_IS_AMBIGUOUS_COLUMN,
    TOP_KINASE_SELECTION_POLICY_COLUMN,
    TOP_KINASE_TIE_COUNT_COLUMN,
    TOP_KINASE_WEIGHTS_COLUMN,
    TOP_SCORE_COLUMN,
)
from phospy.science.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)
from phospy.tables.base import TableSchema
from phospy.tables.datasets import PhosphoIntensityMatrix
from phospy.tables.kinase import KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS
from phospy.workflows.signalome.interpreter import SignalomeWorkflowInterpreter
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import (
    protein_site_key,
    protein_site_key_index,
    site_key_context_columns,
)

_DISPLAY_IDS = ["MAPK14;Y182;", "GSK3B;S9;"]
_GENE_SYMBOLS = ["MAPK14", "GSK3B"]
_SITES = ["Y182", "S9"]
_SITE_INDEX = protein_site_key_index(
    protein_identifiers=_GENE_SYMBOLS,
    sites=_SITES,
)
_SITE_KEYS = _SITE_INDEX.astype(str).tolist()
_AKT1_T308_KEY = protein_site_key(protein_identifier="AKT1", site="T308")
_KINASE_DISPLAY_IDS = [*_DISPLAY_IDS, "AKT1;T308;"]
_KINASE_SITE_INDEX = pd.Index([*_SITE_KEYS, _AKT1_T308_KEY], name="site_key")
_KINASE_SITE_KEYS = _KINASE_SITE_INDEX.astype(str).tolist()
_ALLOW_UNKNOWN_REFERENCE_CONTEXT = (
    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
)
_REPO_ROOT = Path(__file__).resolve().parents[2]
_NATIVE_PANDAS_COPY_ON_WRITE = int(str(pd.__version__).split(".", maxsplit=1)[0]) >= 3
_PANDAS_COPY_ON_WRITE_OPTION_IS_MUTABLE = not _NATIVE_PANDAS_COPY_ON_WRITE
_OBJECT_PAYLOAD_COLUMN = "ownership_payload"
_OBJECT_PAYLOAD_STATE = {
    "list": ("list-start",),
    "dict": ("dict-start",),
    "array": (1.0, 2.0),
    "set": ("set-start",),
    "nested_array": (3.0, 4.0),
    "nested_set": ("nested-set-start",),
    "nested_list": ("nested-list-start",),
}


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [2.0, 1.0],
        },
        index=_SITE_INDEX.copy(),
    )


def _site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": _GENE_SYMBOLS,
            "site": _SITES,
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            ],
            "protein_id": _GENE_SYMBOLS,
            "localisation_confidence": [0.95, 0.9],
            "site_key": _SITE_KEYS,
            "display_id": _DISPLAY_IDS,
            **site_key_context_columns(_SITE_INDEX),
        },
        index=_SITE_INDEX.copy(),
    )


def _mixed_numeric_phospho() -> pd.DataFrame:
    frame = _phospho()
    frame.loc[:, "sample_c"] = np.asarray([3, 4], dtype=np.int64)
    return frame


def _mutable_object_payload() -> dict[str, object]:
    return {
        "list": ["list-start"],
        "dict": {"inner": ["dict-start"]},
        "array": np.asarray([1.0, 2.0]),
        "set": {"set-start"},
        "nested": [
            {"array": np.asarray([3.0, 4.0])},
            {"set": {"nested-set-start"}},
            ["nested-list-start"],
        ],
    }


def _mutate_object_payload(payload: object, marker: str) -> None:
    assert isinstance(payload, dict)
    list_value = payload["list"]
    dict_value = payload["dict"]
    array_value = payload["array"]
    set_value = payload["set"]
    nested_value = payload["nested"]
    assert isinstance(list_value, list)
    assert isinstance(dict_value, dict)
    assert isinstance(array_value, np.ndarray)
    assert isinstance(set_value, set)
    assert isinstance(nested_value, list)
    nested_array_mapping = nested_value[0]
    nested_set_mapping = nested_value[1]
    nested_list = nested_value[2]
    assert isinstance(nested_array_mapping, dict)
    assert isinstance(nested_set_mapping, dict)
    assert isinstance(nested_list, list)
    nested_array = nested_array_mapping["array"]
    nested_set = nested_set_mapping["set"]
    assert isinstance(nested_array, np.ndarray)
    assert isinstance(nested_set, set)

    list_value.append(f"{marker}-list")
    dict_inner = dict_value["inner"]
    assert isinstance(dict_inner, list)
    dict_inner.append(f"{marker}-dict")
    array_value[0] = 99.0
    set_value.add(f"{marker}-set")
    nested_array[0] = 88.0
    nested_set.add(f"{marker}-nested-set")
    nested_list.append(f"{marker}-nested-list")


def _object_payload_state(payload: object) -> dict[str, tuple[object, ...]]:
    assert isinstance(payload, dict)
    nested_value = payload["nested"]
    assert isinstance(nested_value, list)
    nested_array_mapping = nested_value[0]
    nested_set_mapping = nested_value[1]
    nested_list = nested_value[2]
    assert isinstance(nested_array_mapping, dict)
    assert isinstance(nested_set_mapping, dict)
    assert isinstance(nested_list, list)
    array_value = payload["array"]
    nested_array = nested_array_mapping["array"]
    set_value = payload["set"]
    nested_set = nested_set_mapping["set"]
    dict_value = payload["dict"]
    assert isinstance(array_value, np.ndarray)
    assert isinstance(nested_array, np.ndarray)
    assert isinstance(set_value, set)
    assert isinstance(nested_set, set)
    assert isinstance(dict_value, dict)
    dict_inner = dict_value["inner"]
    assert isinstance(dict_inner, list)
    list_value = payload["list"]
    assert isinstance(list_value, list)
    return {
        "list": tuple(list_value),
        "dict": tuple(dict_inner),
        "array": tuple(float(value) for value in array_value.tolist()),
        "set": tuple(sorted(str(value) for value in set_value)),
        "nested_array": tuple(float(value) for value in nested_array.tolist()),
        "nested_set": tuple(sorted(str(value) for value in nested_set)),
        "nested_list": tuple(nested_list),
    }


def _immutable_object_payload_state(payload: object) -> dict[str, tuple[object, ...]]:
    assert isinstance(payload, Mapping)
    nested_value = payload["nested"]
    assert isinstance(nested_value, tuple)
    nested_array_mapping = nested_value[0]
    nested_set_mapping = nested_value[1]
    nested_list = nested_value[2]
    assert isinstance(nested_array_mapping, Mapping)
    assert isinstance(nested_set_mapping, Mapping)
    assert isinstance(nested_list, tuple)
    array_value = payload["array"]
    nested_array = nested_array_mapping["array"]
    set_value = payload["set"]
    nested_set = nested_set_mapping["set"]
    dict_value = payload["dict"]
    assert isinstance(array_value, np.ndarray)
    assert isinstance(nested_array, np.ndarray)
    assert isinstance(set_value, frozenset)
    assert isinstance(nested_set, frozenset)
    assert isinstance(dict_value, Mapping)
    dict_inner = dict_value["inner"]
    assert isinstance(dict_inner, tuple)
    list_value = payload["list"]
    assert isinstance(list_value, tuple)
    return {
        "list": tuple(list_value),
        "dict": tuple(dict_inner),
        "array": tuple(float(value) for value in array_value.tolist()),
        "set": tuple(sorted(str(value) for value in set_value)),
        "nested_array": tuple(float(value) for value in nested_array.tolist()),
        "nested_set": tuple(sorted(str(value) for value in nested_set)),
        "nested_list": tuple(nested_list),
    }


def _object_payload_frame_from_site_metadata(payload: object) -> pd.DataFrame:
    site_metadata = _site_metadata()
    site_metadata.loc[:, _OBJECT_PAYLOAD_COLUMN] = pd.Series(
        [payload, _mutable_object_payload()],
        index=site_metadata.index,
        dtype=object,
    )
    return site_metadata


def _references() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                    "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                ]
            },
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
        ),
    )


def _analysis_ready_dataset() -> AnalysisReadyPhosphoDataset:
    return trusted_analysis_ready_dataset_from_tables(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _differential_workflow_request(
    *,
    n_sites: int = 8,
) -> DifferentialAnalysisRequest:
    genes = [f"GENE{i}" for i in range(n_sites)]
    sites = [f"S{i + 1}" for i in range(n_sites)]
    site_index = protein_site_key_index(
        protein_identifiers=genes,
        sites=sites,
    )
    baseline = np.arange(1, n_sites + 1, dtype=float)
    phospho = pd.DataFrame(
        {
            "A_1": baseline,
            "A_2": baseline + 0.2,
            "B_1": baseline + 1.0,
            "B_2": baseline + 1.2,
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": [
                f"{gene};{site};" for gene, site in zip(genes, sites, strict=True)
            ],
            **site_key_context_columns(site_index),
            "gene_symbol": genes,
            "site": sites,
            "site_sequence": [
                ("A" * 15) + site.strip().upper()[0] + ("A" * 15) for site in sites
            ],
            "protein_id": genes,
        },
        index=site_index.copy(),
    )
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=ExperimentalDesign(
            samples=(
                SampleDesignRecord(
                    sample_id="A_1",
                    condition="A",
                    biological_replicate_id="A_r1",
                ),
                SampleDesignRecord(
                    sample_id="A_2",
                    condition="A",
                    biological_replicate_id="A_r2",
                ),
                SampleDesignRecord(
                    sample_id="B_1",
                    condition="B",
                    biological_replicate_id="B_r1",
                ),
                SampleDesignRecord(
                    sample_id="B_2",
                    condition="B",
                    biological_replicate_id="B_r2",
                ),
            ),
        ),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
    )


def _with_numeric_payload_column(frame: pd.DataFrame) -> pd.DataFrame:
    copied = frame.copy(deep=True)
    copied.loc[:, "numeric_payload"] = np.arange(
        1,
        len(copied.index) + 1,
        dtype=float,
    )
    return copied


def _with_object_payload_column(
    frame: pd.DataFrame,
    payload: object,
) -> pd.DataFrame:
    copied = frame.copy(deep=True)
    payloads = [payload, *[_mutable_object_payload() for _ in range(len(copied) - 1)]]
    copied.loc[:, _OBJECT_PAYLOAD_COLUMN] = pd.Series(
        payloads,
        index=copied.index,
        dtype=object,
    )
    return copied


def _first_numeric_column(frame: pd.DataFrame) -> str:
    for column in frame.columns:
        if pd.api.types.is_numeric_dtype(frame.loc[:, column]):
            return str(column)
    raise AssertionError("registry frame has no numeric column")


def _mutate_first_numeric_cell(frame: pd.DataFrame, value: float) -> None:
    column = _first_numeric_column(frame)
    frame.loc[frame.index[0], column] = value


def _first_object_payload(frame: pd.DataFrame) -> object:
    return frame.loc[frame.index[0], _OBJECT_PAYLOAD_COLUMN]


def _differential_result_table() -> pd.DataFrame:
    table = _site_metadata().copy(deep=True)
    table.loc[:, "logFC"] = [1.0, -1.0]
    table.loc[:, "t"] = [2.0, -2.0]
    table.loc[:, "P.Value"] = [0.05, 0.10]
    table.loc[:, "adj.P.Val"] = [0.10, 0.10]
    table.loc[:, "result_status"] = ["tested", "tested"]
    table.loc[:, "result_status_reason"] = ["", ""]
    return _with_numeric_payload_column(table)


def _differential_prior_diagnostics(
    index: pd.Index,
) -> EmpiricalBayesPriorDiagnostics:
    return EmpiricalBayesPriorDiagnostics(
        method="standard",
        robust=False,
        trend=False,
        winsor_tail_p=(0.05, 0.1),
        base_prior_variance=1.0,
        base_prior_degrees_of_freedom=10.0,
        robust_outlier_count=0,
        robust_outlier_fraction=0.0,
        winsorized_low_count=0,
        winsorized_high_count=0,
        prior_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="prior_residual_variance",
        ),
        prior_degrees_of_freedom=pd.Series(
            np.full(index.size, 10.0),
            index=index.copy(),
            name="prior_degrees_of_freedom",
        ),
    )


def _differential_result_from_table(table: pd.DataFrame) -> DifferentialAnalysisResult:
    index = table.index.copy()
    return DifferentialAnalysisResult(
        residual_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="residual_variance",
        ),
        posterior_residual_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="posterior_residual_variance",
        ),
        prior_residual_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="prior_residual_variance",
        ),
        prior_degrees_of_freedom_series_value=pd.Series(
            np.full(index.size, 10.0),
            index=index.copy(),
            name="prior_degrees_of_freedom",
        ),
        prior_variance=1.0,
        prior_degrees_of_freedom=10.0,
        residual_degrees_of_freedom=4.0,
        empirical_bayes_method="standard",
        empirical_bayes_robust=False,
        empirical_bayes_trend=False,
        prior_diagnostics=_differential_prior_diagnostics(index),
        mean_variance_trend_diagnostics=None,
        contrast_tables={"B_vs_A": table},
    )


def _kinase_score_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "MAP2K6": [0.8, 0.2],
            "AKT1": [0.2, 0.8],
        },
        index=_SITE_INDEX.copy(),
    )


def _prediction_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [0.2, 0.1],
        },
        index=_SITE_INDEX.copy(),
    )


def _activity_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {"MAP2K6": [1.0, 2.0], "AKT1": [0.5, 1.5]},
        index=pd.Index(["sample_a", "sample_b"]),
    )


def _activity_target_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_id": ["MAPK14;Y182;", "GSK3B;S9;"],
            "site_key": _SITE_KEYS,
            "display_id": _DISPLAY_IDS,
            "kinase": ["MAP2K6", "AKT1"],
            "score": [0.9, 0.8],
        }
    )


def _activity_result_from_matrix(matrix: pd.DataFrame) -> KinaseActivityResult:
    activity_input = ActivityInputMatrix.sample_level_abundance(matrix)
    return KinaseActivityResult(
        activity_matrix=matrix,
        thresholded_substrate_mean_activity=_activity_matrix(),
        thresholded_substrate_counts=pd.Series(
            [2, 2],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_substrates",
        ),
        target_counts=pd.Series(
            [1, 1],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_targets",
        ),
        target_table=_activity_target_table(),
        input_semantics=activity_input.semantics,
        profile_metadata=activity_input.profile_metadata,
    )


def _activity_result_from_target_table(table: pd.DataFrame) -> KinaseActivityResult:
    activity_matrix = _activity_matrix()
    activity_input = ActivityInputMatrix.sample_level_abundance(activity_matrix)
    return KinaseActivityResult(
        activity_matrix=activity_matrix,
        thresholded_substrate_mean_activity=_activity_matrix(),
        thresholded_substrate_counts=pd.Series(
            [2, 2],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_substrates",
        ),
        target_counts=pd.Series(
            [1, 1],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_targets",
        ),
        target_table=table,
        input_semantics=activity_input.semantics,
        profile_metadata=activity_input.profile_metadata,
    )


def _activity_result_from_thresholded_counts(
    series: pd.Series,
) -> KinaseActivityResult:
    activity_matrix = _activity_matrix()
    activity_input = ActivityInputMatrix.sample_level_abundance(activity_matrix)
    return KinaseActivityResult(
        activity_matrix=activity_matrix,
        thresholded_substrate_mean_activity=_activity_matrix(),
        thresholded_substrate_counts=series,
        target_counts=pd.Series(
            [1, 1],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_targets",
        ),
        target_table=_activity_target_table(),
        input_semantics=activity_input.semantics,
        profile_metadata=activity_input.profile_metadata,
    )


def _substrate_list() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "kinase": ["MAP2K6", "AKT1"],
            "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            "site_key": _SITE_KEYS,
            "display_id": _DISPLAY_IDS,
            "score": [0.9, 0.8],
            "rank": [1, 1],
        }
    )


def _kinase_substrate_contribution_table() -> pd.DataFrame:
    return pd.DataFrame.from_records(
        [
            {
                "kinase": "MAP2K6",
                "substrate_site": "MAPK14;Y182;",
                "substrate_identifier": "MAPK14;Y182;",
                "value_used_in_scoring": 0.8,
                "score_component": "rank_weighted_fusion_scores",
                "score_source": "profile_only_motif_missing_or_constant",
                "reference_source_name": "fixture",
                "reference_source_version": "v1",
                "reference_bundle_id": "fixture_bundle",
                "reference_identifier_namespace": "display_id",
                "status": "included",
                "exclusion_reason": None,
                "ambiguous": False,
            }
        ],
        columns=pd.Index(KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS),
    )


def _kinase_workflow_result_from_contributions(
    table: pd.DataFrame,
) -> KinaseWorkflowResult:
    return KinaseWorkflowResult(
        dataset=_analysis_ready_dataset(),
        references=_references(),
        scoring_result=KinaseScoringResult(profile_scores=_kinase_score_matrix()),
        prediction_result=KinasePredictionResult(pred_mat=_prediction_matrix()),
        substrate_contributions=table,
    )


def _empty_signalome_assignments_table() -> pd.DataFrame:
    columns = (
        SITE_KEY_COLUMN,
        DISPLAY_ID_COLUMN,
        GENE_SYMBOL_COLUMN,
        SITE_COLUMN,
        PROTEIN_COLUMN,
        PROTEIN_ACCESSION_COLUMN,
        ISOFORM_ID_COLUMN,
        MODULE_ID_COLUMN,
        TOP_KINASE_COLUMN,
        TOP_SCORE_COLUMN,
        TOP_KINASE_CANDIDATES_COLUMN,
        TOP_KINASE_WEIGHTS_COLUMN,
        TOP_KINASE_TIE_COUNT_COLUMN,
        TOP_KINASE_IS_AMBIGUOUS_COLUMN,
        TOP_KINASE_SELECTION_POLICY_COLUMN,
        MODULE_TOP_KINASE_COLUMN,
        MODULE_TOP_KINASE_CANDIDATES_COLUMN,
        MODULE_TOP_KINASE_TIE_COUNT_COLUMN,
        MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN,
        MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN,
    )
    return pd.DataFrame(columns=columns, index=_SITE_INDEX[:0].copy())


def _empty_kinase_network_edges() -> pd.DataFrame:
    return pd.DataFrame(
        columns=(
            SOURCE_KINASE_COLUMN,
            TARGET_KINASE_COLUMN,
            CORRELATION_COLUMN,
            "valid_observations",
        )
    )


def _empty_signalome_modules_table() -> pd.DataFrame:
    return pd.DataFrame({"MAP2K6": pd.Series(dtype="float64")})


def _signalome_expanded_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            SITE_KEY_COLUMN: _SITE_KEYS,
            DISPLAY_ID_COLUMN: _DISPLAY_IDS,
            "numeric_payload": [1.0, 2.0],
        },
        index=_SITE_INDEX.copy(),
    )


def _signalome_result_from_expanded_table(
    table: pd.DataFrame,
) -> SignalomeWorkflowResult:
    return SignalomeWorkflowResult(
        dataset=_analysis_ready_dataset(),
        kinase_result=_kinase_workflow_result_from_contributions(
            _kinase_substrate_contribution_table()
        ),
        module_assignments=SignalomeAssignments._from_owned(
            table=_empty_signalome_assignments_table()
        ),
        signalome_modules=SignalomeModules._from_owned(
            table=_empty_signalome_modules_table()
        ),
        kinase_network=KinaseNetwork._from_owned(edges=_empty_kinase_network_edges()),
        expanded_signalome=table,
    )


def _reference_kinase_map() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "kinase": ["MAP2K6", "AKT1"],
            "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            "numeric_payload": [1.0, 2.0],
        }
    )


@dataclass(frozen=True, slots=True)
class _PublicFrameOwnerCase:
    name: str
    make_numeric_source: Callable[[], pd.DataFrame]
    construct_from_numeric: Callable[[pd.DataFrame], object]
    observe_numeric: Callable[[object], pd.DataFrame]
    make_object_source: Callable[[object], pd.DataFrame] | None = None
    construct_from_object: Callable[[pd.DataFrame], object] | None = None
    observe_object_payload: Callable[[object], object] | None = None


def _public_frame_owner_cases() -> tuple[_PublicFrameOwnerCase, ...]:
    return (
        _PublicFrameOwnerCase(
            name="analysis-ready-dataset",
            make_numeric_source=_phospho,
            construct_from_numeric=lambda frame: (
                trusted_analysis_ready_dataset_from_tables(
                    phospho=frame,
                    site_metadata=_site_metadata(),
                    organism=Organism.RAT,
                    intensity_scale_state=supported_linear_intensity_scale_state(
                        has_total_matrix=False
                    ),
                    processing_state=supported_linear_processing_state(
                        has_total_matrix=False
                    ),
                )
            ),
            observe_numeric=lambda owner: owner.to_dataframe(),
            make_object_source=lambda payload: _object_payload_frame_from_site_metadata(
                payload
            ),
            construct_from_object=lambda frame: (
                trusted_analysis_ready_dataset_from_tables(
                    phospho=_phospho(),
                    site_metadata=frame,
                    organism=Organism.RAT,
                    intensity_scale_state=supported_linear_intensity_scale_state(
                        has_total_matrix=False
                    ),
                    processing_state=supported_linear_processing_state(
                        has_total_matrix=False
                    ),
                )
            ),
            observe_object_payload=lambda owner: _first_object_payload(
                owner.site_metadata_dataframe()
            ),
        ),
        _PublicFrameOwnerCase(
            name="phosphosite-import-result",
            make_numeric_source=_phospho,
            construct_from_numeric=lambda frame: PhosphositeImportResult(
                phospho_matrix_candidate=frame,
                site_metadata_candidate=_site_metadata(),
                sample_column_mapping={"sample_a": "sample_a"},
            ),
            observe_numeric=lambda owner: owner.phospho_matrix_candidate,
            make_object_source=lambda payload: _with_object_payload_column(
                _site_metadata(),
                payload,
            ),
            construct_from_object=lambda frame: PhosphositeImportResult(
                phospho_matrix_candidate=_phospho(),
                site_metadata_candidate=frame,
                sample_column_mapping={"sample_a": "sample_a"},
            ),
            observe_object_payload=lambda owner: _first_object_payload(
                owner.site_metadata_candidate
            ),
        ),
        _PublicFrameOwnerCase(
            name="differential-result",
            make_numeric_source=_differential_result_table,
            construct_from_numeric=_differential_result_from_table,
            observe_numeric=lambda owner: owner.table_for("B_vs_A"),
            make_object_source=lambda payload: _with_object_payload_column(
                _differential_result_table(),
                payload,
            ),
            construct_from_object=_differential_result_from_table,
            observe_object_payload=lambda owner: _first_object_payload(
                owner.table_for("B_vs_A")
            ),
        ),
        _PublicFrameOwnerCase(
            name="kinase-activity-result",
            make_numeric_source=_activity_matrix,
            construct_from_numeric=_activity_result_from_matrix,
            observe_numeric=lambda owner: owner.activity_matrix,
            make_object_source=lambda payload: _with_object_payload_column(
                _activity_target_table(),
                payload,
            ),
            construct_from_object=_activity_result_from_target_table,
            observe_object_payload=lambda owner: _first_object_payload(
                owner.target_table
            ),
        ),
        _PublicFrameOwnerCase(
            name="kinase-scoring-result",
            make_numeric_source=_kinase_score_matrix,
            construct_from_numeric=lambda frame: KinaseScoringResult(
                profile_scores=frame
            ),
            observe_numeric=lambda owner: owner.profile_scores,
        ),
        _PublicFrameOwnerCase(
            name="kinase-prediction-result",
            make_numeric_source=_prediction_matrix,
            construct_from_numeric=lambda frame: KinasePredictionResult(pred_mat=frame),
            observe_numeric=lambda owner: owner.pred_mat,
            make_object_source=lambda payload: _with_object_payload_column(
                _substrate_list(),
                payload,
            ),
            construct_from_object=lambda frame: KinasePredictionResult(
                pred_mat=_prediction_matrix(),
                substrate_list=frame,
            ),
            observe_object_payload=lambda owner: _first_object_payload(
                owner.substrate_list
            ),
        ),
        _PublicFrameOwnerCase(
            name="kinase-workflow-result",
            make_numeric_source=_kinase_substrate_contribution_table,
            construct_from_numeric=_kinase_workflow_result_from_contributions,
            observe_numeric=lambda owner: owner.substrate_contributions,
        ),
        _PublicFrameOwnerCase(
            name="signalome-workflow-result",
            make_numeric_source=_signalome_expanded_table,
            construct_from_numeric=_signalome_result_from_expanded_table,
            observe_numeric=lambda owner: owner.to_dataframe(),
            make_object_source=lambda payload: _with_object_payload_column(
                _signalome_expanded_table(),
                payload,
            ),
            construct_from_object=_signalome_result_from_expanded_table,
            observe_object_payload=lambda owner: _first_object_payload(
                owner.to_dataframe()
            ),
        ),
        _PublicFrameOwnerCase(
            name="reference-bundle",
            make_numeric_source=_reference_kinase_map,
            construct_from_numeric=lambda frame: ReferenceBundle(
                organism=Organism.RAT,
                kinase_substrate_map=frame,
                site_sequences=_references().site_sequences_dataframe(),
            ),
            observe_numeric=lambda owner: owner.kinase_substrate_map_dataframe(),
            make_object_source=lambda payload: _with_object_payload_column(
                _reference_kinase_map(),
                payload,
            ),
            construct_from_object=lambda frame: ReferenceBundle(
                organism=Organism.RAT,
                kinase_substrate_map=frame,
                site_sequences=_references().site_sequences_dataframe(),
            ),
            observe_object_payload=lambda owner: _first_object_payload(
                owner.kinase_substrate_map_dataframe()
            ),
        ),
    )


@dataclass(frozen=True, slots=True)
class _PublicSeriesOwnerCase:
    name: str
    make_source: Callable[[], pd.Series]
    construct: Callable[[pd.Series], object]
    observe: Callable[[object], pd.Series]


def _public_series_owner_cases() -> tuple[_PublicSeriesOwnerCase, ...]:
    return (
        _PublicSeriesOwnerCase(
            name="differential-residual-variance",
            make_source=lambda: pd.Series(
                [1.0, 2.0],
                index=_SITE_INDEX.copy(),
                name="residual_variance",
            ),
            construct=lambda series: DifferentialAnalysisResult(
                residual_variance=series,
                posterior_residual_variance=pd.Series(
                    [1.0, 2.0],
                    index=_SITE_INDEX.copy(),
                    name="posterior_residual_variance",
                ),
                prior_residual_variance=pd.Series(
                    [1.0, 2.0],
                    index=_SITE_INDEX.copy(),
                    name="prior_residual_variance",
                ),
                prior_degrees_of_freedom_series_value=pd.Series(
                    [10.0, 10.0],
                    index=_SITE_INDEX.copy(),
                    name="prior_degrees_of_freedom",
                ),
                prior_variance=1.0,
                prior_degrees_of_freedom=10.0,
                residual_degrees_of_freedom=4.0,
                empirical_bayes_method="standard",
                empirical_bayes_robust=False,
                empirical_bayes_trend=False,
                prior_diagnostics=_differential_prior_diagnostics(_SITE_INDEX.copy()),
                mean_variance_trend_diagnostics=None,
                contrast_tables={"B_vs_A": _differential_result_table()},
            ),
            observe=lambda owner: owner.residual_variance_series(),
        ),
        _PublicSeriesOwnerCase(
            name="kinase-activity-thresholded-counts",
            make_source=lambda: pd.Series(
                [2, 2],
                index=pd.Index(["MAP2K6", "AKT1"]),
                name="n_substrates",
            ),
            construct=_activity_result_from_thresholded_counts,
            observe=lambda owner: owner.thresholded_substrate_counts,
        ),
    )


@pytest.mark.parametrize(
    "case",
    _public_frame_owner_cases(),
    ids=lambda case: case.name,
)
def test_public_dataframe_bearing_types_isolate_caller_inputs(
    case: _PublicFrameOwnerCase,
) -> None:
    source = case.make_numeric_source()
    owner = case.construct_from_numeric(source)
    before = case.observe_numeric(owner).copy(deep=True)

    _mutate_first_numeric_cell(source, 999.0)

    pd.testing.assert_frame_equal(case.observe_numeric(owner), before)


@pytest.mark.parametrize(
    "case",
    _public_frame_owner_cases(),
    ids=lambda case: case.name,
)
def test_public_dataframe_bearing_types_isolate_caller_numpy_views(
    case: _PublicFrameOwnerCase,
) -> None:
    source = case.make_numeric_source()
    owner = case.construct_from_numeric(source)
    before = case.observe_numeric(owner).copy(deep=True)
    numeric_column = _first_numeric_column(source)

    assert _force_numeric_array_and_bases_writeable_and_mutate(
        source.loc[:, [numeric_column]].to_numpy(copy=False),
        value=888.0,
    )

    pd.testing.assert_frame_equal(case.observe_numeric(owner), before)


@pytest.mark.parametrize(
    "case",
    _public_frame_owner_cases(),
    ids=lambda case: case.name,
)
def test_public_dataframe_exports_do_not_mutate_owned_state(
    case: _PublicFrameOwnerCase,
) -> None:
    source = case.make_numeric_source()
    owner = case.construct_from_numeric(source)
    before = case.observe_numeric(owner).copy(deep=True)
    exported = case.observe_numeric(owner)

    _mutate_first_numeric_cell(exported, 777.0)

    pd.testing.assert_frame_equal(case.observe_numeric(owner), before)


@pytest.mark.parametrize(
    "case",
    [
        case
        for case in _public_frame_owner_cases()
        if case.make_object_source is not None
        and case.construct_from_object is not None
        and case.observe_object_payload is not None
    ],
    ids=lambda case: case.name,
)
def test_public_dataframe_bearing_types_isolate_nested_object_cells(
    case: _PublicFrameOwnerCase,
) -> None:
    assert case.make_object_source is not None
    assert case.construct_from_object is not None
    assert case.observe_object_payload is not None
    payload = _mutable_object_payload()
    source = case.make_object_source(payload)
    owner = case.construct_from_object(source)

    _mutate_object_payload(payload, "caller")
    observed_payload = case.observe_object_payload(owner)

    assert _object_payload_state(observed_payload) == _OBJECT_PAYLOAD_STATE

    exported_payload = case.observe_object_payload(owner)
    _mutate_object_payload(exported_payload, "export")

    assert (
        _object_payload_state(case.observe_object_payload(owner))
        == _OBJECT_PAYLOAD_STATE
    )


@pytest.mark.parametrize(
    "case",
    _public_series_owner_cases(),
    ids=lambda case: case.name,
)
def test_public_series_bearing_types_isolate_caller_inputs(
    case: _PublicSeriesOwnerCase,
) -> None:
    source = case.make_source()
    owner = case.construct(source)
    before = case.observe(owner).copy(deep=True)

    source.iloc[0] = 999

    pd.testing.assert_series_equal(case.observe(owner), before)


@pytest.mark.parametrize(
    "case",
    _public_series_owner_cases(),
    ids=lambda case: case.name,
)
def test_public_series_bearing_types_isolate_caller_numpy_views(
    case: _PublicSeriesOwnerCase,
) -> None:
    source = case.make_source()
    owner = case.construct(source)
    before = case.observe(owner).copy(deep=True)

    assert _force_numeric_array_and_bases_writeable_and_mutate(
        source.to_numpy(copy=False),
        value=888.0,
    )

    pd.testing.assert_series_equal(case.observe(owner), before)


@pytest.mark.parametrize(
    "case",
    _public_series_owner_cases(),
    ids=lambda case: case.name,
)
def test_public_series_exports_do_not_mutate_owned_state(
    case: _PublicSeriesOwnerCase,
) -> None:
    source = case.make_source()
    owner = case.construct(source)
    before = case.observe(owner).copy(deep=True)
    exported = case.observe(owner)

    exported.iloc[0] = 777

    pd.testing.assert_series_equal(case.observe(owner), before)


def _kinase_result():
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 4.0],
            "sample_b": [2.0, 4.0, 1.0],
        },
        index=_KINASE_SITE_INDEX.copy(),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "S9", "T308"]
            ],
            "protein_id": ["MAPK14", "GSK3B", "AKT1"],
            "localisation_confidence": [0.95, 0.9, 0.92],
            "site_key": _KINASE_SITE_KEYS,
            "display_id": _KINASE_DISPLAY_IDS,
            **site_key_context_columns(_KINASE_SITE_INDEX),
        },
        index=_KINASE_SITE_INDEX.copy(),
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6", "AKT1", "AKT1"],
                "substrate_site": [
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                    "GSK3B;S9;",
                    "AKT1;T308;",
                ],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31, "B" * 31, "C" * 31]},
            index=pd.Index(_KINASE_DISPLAY_IDS, name="site_id"),
        ),
    )
    return KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=trusted_analysis_ready_dataset_from_tables(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                intensity_scale_state=supported_linear_intensity_scale_state(
                    has_total_matrix=False
                ),
                processing_state=supported_linear_processing_state(
                    has_total_matrix=False
                ),
            ),
            references=references,
            scoring_config=KinaseScoringConfig(
                reliability_profile="custom",
                min_substrates=2,
                reference_context_compatibility_policy=(
                    _ALLOW_UNKNOWN_REFERENCE_CONTEXT
                ),
            ),
            prediction_config=KinasePredictionConfig(
                top_k=1,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=None,
        )
    )


def _signalome_request_for_read_path_mutation_checks() -> SignalomeWorkflowRequest:
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    prediction_matrix = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [0.1, 0.2],
        },
        index=dataset._borrow_phospho_frame().index.copy(),
    )
    score_matrix = pd.DataFrame(
        {
            "MAP2K6": [1.5, 1.2],
            "AKT1": [0.6, 0.7],
        },
        index=dataset._borrow_phospho_frame().index.copy(),
    )
    return SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=dataset,
            references=_references(),
            scoring_result=KinaseScoringResult(
                profile_scores=score_matrix,
                rank_weighted_fusion_scores=score_matrix,
            ),
            prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
            activity_result=None,
        ),
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            reference_context_compatibility_policy=(_ALLOW_UNKNOWN_REFERENCE_CONTEXT),
        ),
    )


@dataclass(slots=True)
class _CopyCounts:
    dataframe_deep: int = 0


@dataclass(slots=True)
class _FullMatrixCopyCounts:
    full_matrix_deep: int = 0


@contextmanager
def _count_dataframe_deep_copies() -> Iterator[_CopyCounts]:
    counts = _CopyCounts()
    original_copy = pd.DataFrame.copy

    def wrapped_copy(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        deep = kwargs.get("deep", args[0] if args else True)
        if bool(deep):
            counts.dataframe_deep += 1
        return original_copy(self, *args, **kwargs)

    pd.DataFrame.copy = wrapped_copy
    try:
        yield counts
    finally:
        pd.DataFrame.copy = original_copy


@contextmanager
def _count_full_matrix_deep_copies(
    *,
    shape: tuple[int, int],
    columns: tuple[object, ...],
) -> Iterator[_FullMatrixCopyCounts]:
    counts = _FullMatrixCopyCounts()
    original_copy = pd.DataFrame.copy

    def wrapped_copy(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        deep = kwargs.get("deep", args[0] if args else True)
        if bool(deep) and self.shape == shape and tuple(self.columns) == columns:
            counts.full_matrix_deep += 1
        return original_copy(self, *args, **kwargs)

    pd.DataFrame.copy = wrapped_copy
    try:
        yield counts
    finally:
        pd.DataFrame.copy = original_copy


def _assert_numpy_blocks_readonly(frame: pd.DataFrame | pd.Series) -> None:
    blocks = frame._mgr.blocks
    for block in blocks:
        values = block.values
        flags = getattr(values, "flags", None)
        if flags is not None and hasattr(flags, "writeable"):
            assert flags.writeable is False


def _mutate_first_frame_cell(frame: pd.DataFrame) -> None:
    if pd.api.types.is_numeric_dtype(frame.dtypes.iloc[0]):
        frame.iloc[0, 0] = float(frame.iloc[0, 0]) + 1.0
    else:
        frame.iloc[0, 0] = f"{frame.iloc[0, 0]}_changed"


def _mutate_existing_borrowed_frame_cell(
    *,
    owner: pd.DataFrame,
    borrowed: pd.DataFrame,
    value: object,
) -> None:
    before_value = owner.iloc[0, 0]
    before_shape = owner.shape
    before_columns = tuple(owner.columns)
    try:
        borrowed.iloc[0, 0] = value
    except (TypeError, ValueError):
        pass
    assert owner.iloc[0, 0] == before_value
    assert owner.shape == before_shape
    assert tuple(owner.columns) == before_columns


def _mutate_existing_borrowed_series_value(
    *,
    owner: pd.Series,
    borrowed: pd.Series,
    value: object,
) -> None:
    before_value = owner.iloc[0]
    before_shape = owner.shape
    try:
        borrowed.iloc[0] = value
    except (TypeError, ValueError):
        pass
    assert owner.iloc[0] == before_value
    assert owner.shape == before_shape


def _numeric_block_count(frame: pd.DataFrame) -> int:
    manager = getattr(frame, "_mgr", None)
    blocks = getattr(manager, "blocks", ())
    count = 0
    for block in blocks:
        values = getattr(block, "values", None)
        if isinstance(values, np.ndarray) and np.issubdtype(values.dtype, np.number):
            count += 1
    return count


def _force_numeric_array_writeable_and_mutate(
    array: object,
    *,
    value: float,
) -> bool:
    if not isinstance(array, np.ndarray):
        return False
    if not np.issubdtype(array.dtype, np.number):
        return False
    try:
        array.setflags(write=True)
    except ValueError:
        pass
    try:
        array[...] = value
    except (TypeError, ValueError, OverflowError):
        return False
    return True


def _force_numeric_array_and_bases_writeable_and_mutate(
    array: object,
    *,
    value: float,
) -> bool:
    mutated = _force_numeric_array_writeable_and_mutate(array, value=value)
    seen: set[int] = set()
    base = getattr(array, "base", None)
    while isinstance(base, np.ndarray) and id(base) not in seen:
        seen.add(id(base))
        mutated = (
            _force_numeric_array_writeable_and_mutate(base, value=value) or mutated
        )
        base = getattr(base, "base", None)
    return mutated


def _mutate_numeric_slice_values(frame: pd.DataFrame) -> bool:
    numeric_columns = [
        column
        for column in frame.columns
        if pd.api.types.is_numeric_dtype(frame.loc[:, column])
    ]
    assert numeric_columns
    slice_frame = frame.loc[:, numeric_columns[:1]]
    return _force_numeric_array_and_bases_writeable_and_mutate(
        slice_frame.values,
        value=444.0,
    )


def _mutate_numeric_values(frame: pd.DataFrame) -> bool:
    return _force_numeric_array_and_bases_writeable_and_mutate(
        frame.values,
        value=555.0,
    )


def _mutate_numeric_to_numpy_copy_false(frame: pd.DataFrame) -> bool:
    return _force_numeric_array_and_bases_writeable_and_mutate(
        frame.to_numpy(copy=False),
        value=666.0,
    )


def _mutate_numeric_blocks_and_bases(frame: pd.DataFrame) -> bool:
    expected_count = _numeric_block_count(frame)
    assert expected_count >= 1
    manager = getattr(frame, "_mgr", None)
    blocks = getattr(manager, "blocks", ())
    mutated_count = 0
    for block in blocks:
        values = getattr(block, "values", None)
        if _force_numeric_array_and_bases_writeable_and_mutate(
            values,
            value=777.0,
        ):
            mutated_count += 1
    assert mutated_count == expected_count
    return True


def _try_borrowed_column_write(frame: pd.DataFrame) -> None:
    try:
        frame.loc[:, "borrowed_only"] = [1.0, 2.0]
    except (TypeError, ValueError):
        pass


def _attribute_path(node: ast.AST) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parent = _attribute_path(node.value)
        if parent is None:
            return None
        return f"{parent}.{node.attr}"
    return None


def _assignment_targets(node: ast.AST) -> tuple[ast.AST, ...]:
    if isinstance(node, ast.Assign):
        return tuple(node.targets)
    if isinstance(node, ast.AnnAssign | ast.AugAssign):
        return (node.target,)
    return ()


def _pandas4_warning_type() -> type[Warning]:
    return getattr(pd.errors, "Pandas4Warning", Warning)


def _pandas_global_option_snapshot() -> dict[str, object]:
    option_names = ["mode.chained_assignment"]
    if _PANDAS_COPY_ON_WRITE_OPTION_IS_MUTABLE:
        option_names.append("mode.copy_on_write")
    return {name: pd.get_option(name) for name in option_names}


def _assert_dataframe_getter_defensive_snapshot(
    getter: Callable[[], pd.DataFrame],
) -> None:
    exported = getter()
    _mutate_first_frame_cell(exported)
    reread = getter()
    assert exported is not reread
    assert exported.iloc[0, 0] != reread.iloc[0, 0]


def _assert_optional_dataframe_getter_defensive_snapshot(
    getter: Callable[[], pd.DataFrame | None],
) -> None:
    exported = getter()
    assert exported is not None
    _mutate_first_frame_cell(exported)
    reread = getter()
    assert reread is not None
    assert exported is not reread
    assert exported.iloc[0, 0] != reread.iloc[0, 0]


def _assert_optional_dataframe_getter_defensive_snapshot_when_present(
    getter: Callable[[], pd.DataFrame | None],
) -> None:
    exported = getter()
    if exported is None:
        return

    _mutate_first_frame_cell(exported)
    reread = getter()
    assert reread is not None
    assert exported is not reread
    assert exported.iloc[0, 0] != reread.iloc[0, 0]


def _assert_copy_keyword_rejected(
    export: Callable[..., object],
) -> None:
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        export(copy=False)  # type: ignore[call-arg]


@pytest.mark.parametrize(
    (
        "builder",
        "mutated_phospho",
        "mutated_site_metadata",
        "expected_phospho",
        "expected_gene",
    ),
    [
        pytest.param(
            lambda p, s: trusted_analysis_ready_dataset_from_tables(
                phospho=p,
                site_metadata=s,
                organism=Organism.RAT,
                intensity_scale_state=supported_linear_intensity_scale_state(
                    has_total_matrix=False
                ),
                processing_state=supported_linear_processing_state(
                    has_total_matrix=False
                ),
            ),
            (0, 0, 999.0),
            (0, 0, "CHANGED"),
            1.0,
            "MAPK14",
            id="public-constructor-copies-caller-inputs",
        ),
        pytest.param(
            lambda p, s: AnalysisReadyDatasetBuilder().run(
                DatasetBuildRequest(
                    phospho=p,
                    site_metadata=s,
                    organism=Organism.RAT,
                    input_intensity_scale="linear",
                )
            ),
            (1, 1, 777.0),
            (1, 0, "CHANGED"),
            1.0,
            "GSK3B",
            id="builder-result-copies-caller-inputs",
        ),
    ],
)
def test_public_constructor_copy_contract_matrix(
    builder: Callable[[pd.DataFrame, pd.DataFrame], AnalysisReadyPhosphoDataset],
    mutated_phospho: tuple[int, int, float],
    mutated_site_metadata: tuple[int, int, str],
    expected_phospho: float,
    expected_gene: str,
) -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    built = builder(phospho, site_metadata)

    phospho.iloc[mutated_phospho[0], mutated_phospho[1]] = mutated_phospho[2]
    site_metadata.iloc[mutated_site_metadata[0], mutated_site_metadata[1]] = (
        mutated_site_metadata[2]
    )

    assert float(built.phospho.iloc[mutated_phospho[0], mutated_phospho[1]]) == (
        expected_phospho
    )
    assert (
        str(
            built.site_metadata.iloc[mutated_site_metadata[0], mutated_site_metadata[1]]
        )
        == expected_gene
    )


@pytest.mark.parametrize(
    ("case_name", "owner_factory"),
    [
        pytest.param(
            "dataset-site-metadata",
            lambda payload: (
                lambda dataset: dataset.site_metadata,
                trusted_analysis_ready_dataset_from_tables(
                    phospho=_phospho(),
                    site_metadata=_object_payload_frame_from_site_metadata(payload),
                    organism=Organism.RAT,
                    intensity_scale_state=supported_linear_intensity_scale_state(
                        has_total_matrix=False
                    ),
                    processing_state=supported_linear_processing_state(
                        has_total_matrix=False
                    ),
                ),
                _SITE_INDEX[0],
            ),
            id="dataset-site-metadata",
        ),
        pytest.param(
            "prediction-result-substrate-list",
            lambda payload: (
                lambda result: result.substrate_list,
                KinasePredictionResult(
                    pred_mat=pd.DataFrame(
                        {"MAP2K6": [0.9, 0.8]},
                        index=_SITE_INDEX.copy(),
                    ),
                    substrate_list=pd.DataFrame(
                        {
                            _OBJECT_PAYLOAD_COLUMN: pd.Series(
                                [payload],
                                dtype=object,
                            )
                        }
                    ),
                ),
                0,
            ),
            id="prediction-result-substrate-list",
        ),
    ],
)
def test_public_construction_and_exports_recursively_isolate_mutable_object_cells(
    case_name: str,
    owner_factory: Callable[
        [dict[str, object]],
        tuple[Callable[[object], pd.DataFrame | None], object, object],
    ],
) -> None:
    payload = _mutable_object_payload()
    frame_getter, owner, row_label = owner_factory(payload)

    _mutate_object_payload(payload, "caller")
    exported_after_caller_mutation = frame_getter(owner)
    assert exported_after_caller_mutation is not None
    owner_payload = exported_after_caller_mutation.loc[
        row_label,
        _OBJECT_PAYLOAD_COLUMN,
    ]
    assert _object_payload_state(owner_payload) == _OBJECT_PAYLOAD_STATE, case_name

    export_one = frame_getter(owner)
    export_two = frame_getter(owner)
    assert export_one is not None
    assert export_two is not None
    export_one_payload = export_one.loc[row_label, _OBJECT_PAYLOAD_COLUMN]
    export_two_payload = export_two.loc[row_label, _OBJECT_PAYLOAD_COLUMN]

    _mutate_object_payload(export_one_payload, "export-one")
    assert _object_payload_state(export_two_payload) == _OBJECT_PAYLOAD_STATE
    owner_reread = frame_getter(owner)
    assert owner_reread is not None
    assert (
        _object_payload_state(owner_reread.loc[row_label, _OBJECT_PAYLOAD_COLUMN])
        == _OBJECT_PAYLOAD_STATE
    )

    _mutate_object_payload(export_two_payload, "export-two")
    assert "export-two-list" not in _object_payload_state(export_one_payload)["list"]


class _UnsupportedMutableObject:
    def __init__(self) -> None:
        self.values: list[str] = []

    def mutate(self) -> None:
        self.values.append("changed")


@pytest.mark.parametrize(
    ("factory", "error_type"),
    [
        pytest.param(
            lambda payload: trusted_analysis_ready_dataset_from_tables(
                phospho=_phospho(),
                site_metadata=_object_payload_frame_from_site_metadata(payload),
                organism=Organism.RAT,
                intensity_scale_state=supported_linear_intensity_scale_state(
                    has_total_matrix=False
                ),
                processing_state=supported_linear_processing_state(
                    has_total_matrix=False
                ),
            ),
            DatasetValidationError,
            id="dataset",
        ),
        pytest.param(
            lambda payload: KinasePredictionResult(
                pred_mat=pd.DataFrame(
                    {"MAP2K6": [0.9, 0.8]},
                    index=_SITE_INDEX.copy(),
                ),
                substrate_list=pd.DataFrame(
                    {
                        _OBJECT_PAYLOAD_COLUMN: pd.Series(
                            [payload],
                            dtype=object,
                        )
                    }
                ),
            ),
            PhosPyValidationError,
            id="prediction-result",
        ),
    ],
)
def test_public_construction_rejects_unsupported_mutable_object_cells(
    factory: Callable[[object], object],
    error_type: type[Exception],
) -> None:
    with pytest.raises(
        error_type,
        match="unsupported mutable object.*ownership_payload",
    ):
        factory(_UnsupportedMutableObject())


def test_builder_stage_handoff_transfers_owned_frames_without_recopies() -> None:
    request = DatasetBuildRequest(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        input_intensity_scale="linear",
    )
    interpreted = DatasetBuildRequestInterpreter().run(request)
    built = DatasetBuildExecutor().run(interpreted)

    assert built._phospho is interpreted.phospho
    assert built._site_metadata is interpreted.site_metadata


def test_builder_dataframe_copy_churn_regression_budget() -> None:
    request = DatasetBuildRequest(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        input_intensity_scale="linear",
    )

    with _count_dataframe_deep_copies() as counts:
        AnalysisReadyDatasetBuilder().run(request)

    assert counts.dataframe_deep == 2


def test_internal_activity_inputs_alias_owned_frames() -> None:
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [0.2, 0.1],
        },
        index=_SITE_INDEX.copy(),
    )
    phospho_matrix = _phospho()
    overlap_summary = PredMatOverlapSummary(
        overlap_count=2,
        pred_mat_rows=2,
        phospho_rows=2,
    )

    inputs = KinaseActivityInputs(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=0.5,
        min_substrates=1,
        top_n_substrates=2,
        overlap_summary=overlap_summary,
        activity_input=ActivityInputMatrix.sample_level_abundance(
            phospho_matrix,
            _assume_owned=True,
        ),
    )

    assert inputs.pred_mat is pred_mat
    assert inputs.phospho_matrix is phospho_matrix


def test_prediction_result_boundary_copy_and_owned_transfer_modes() -> None:
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [0.2, 0.1],
        },
        index=_SITE_INDEX.copy(),
    )
    substrate_list = pd.DataFrame(
        {
            "kinase": ["MAP2K6"],
            "substrate_site": ["MAPK14;Y182;"],
            "score": [0.75],
            "rank": [1],
        }
    )

    copied_result = KinasePredictionResult(
        pred_mat=pred_mat,
        substrate_list=substrate_list,
    )
    owned_result = KinasePredictionResult._from_owned(
        pred_mat=pred_mat,
        substrate_list=substrate_list,
    )

    assert copied_result._pred_mat is not pred_mat
    assert copied_result._substrate_list is not substrate_list
    assert owned_result._pred_mat is pred_mat
    assert owned_result._substrate_list is substrate_list

    pred_mat.iloc[0, 0] = 999.0
    substrate_list.iloc[0, 0] = "CHANGED"
    assert float(copied_result.pred_mat.iloc[0, 0]) == 0.9
    assert str(copied_result.substrate_list.iloc[0, 0]) == "MAP2K6"


def test_public_dataframe_accessors_do_not_accept_copy_keyword() -> None:
    public_accessors = (
        (AnalysisReadyPhosphoDataset, "to_dataframe"),
        (AnalysisReadyPhosphoDataset, "site_metadata_dataframe"),
        (AnalysisReadyPhosphoDataset, "sample_metadata_dataframe"),
        (AnalysisReadyPhosphoDataset, "total_dataframe"),
        (AnalysisReadyPhosphoDataset, "comparisons_dataframe"),
        (AnalysisReadyPhosphoDataset, "imputation_feature_metadata_dataframe"),
        (AnalysisReadyPhosphoDataset, "imputation_observed_mask_dataframe"),
        (DatasetPreprocessingReport, "row_counts_dataframe"),
        (DatasetPreprocessingReport, "operations_dataframe"),
        (DatasetPreprocessingReport, "row_audit_dataframe"),
        (DatasetPreprocessingReport, "duplicate_site_resolution_dataframe"),
        (DatasetPreprocessingReport, "metadata_conflicts_dataframe"),
        (DatasetPreprocessingReport, "comparison_group_stats_dataframe"),
        (DatasetPreprocessingReport, "comparison_pair_stats_dataframe"),
        (KinaseScoringResult, "to_dataframe"),
        (KinaseScoringResult, "motif_scores_dataframe"),
        (KinaseScoringResult, "rank_weighted_fusion_scores_dataframe"),
        (KinaseScoringResult, "score_fusion_weights_dataframe"),
        (KinasePredictionResult, "to_dataframe"),
        (KinasePredictionResult, "substrate_list_dataframe"),
        (KinaseActivityResult, "to_dataframe"),
        (KinaseActivityResult, "thresholded_substrate_mean_activity_dataframe"),
        (KinaseActivityResult, "target_table_dataframe"),
        (ReferenceBundle, "kinase_substrate_map_dataframe"),
        (ReferenceBundle, "site_sequences_dataframe"),
        (SignalomeWorkflowResult, "to_dataframe"),
        (SignalomeWorkflowResult, "site_membership_dataframe"),
        (SignalomeWorkflowResult, "protein_site_context_dataframe"),
        (SignalomeAssignments, "to_pandas"),
        (SignalomeModules, "to_pandas"),
        (KinaseNetwork, "to_pandas"),
        (KinaseNetwork, "nodes_dataframe"),
        (KinaseNetwork, "candidate_correlations_dataframe"),
        (TableSchema, "to_pandas"),
    )

    for owner, method_name in public_accessors:
        signature = inspect.signature(getattr(owner, method_name))
        assert "copy" not in signature.parameters


def test_source_tree_has_no_pandas_global_option_assignments() -> None:
    offenders: list[str] = []
    for path in sorted((_REPO_ROOT / "src" / "phospy").rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            for target in _assignment_targets(node):
                target_path = _attribute_path(target)
                if target_path is None:
                    continue
                if target_path.startswith(("pd.options", "pandas.options")):
                    line_number = getattr(node, "lineno", "?")
                    offenders.append(
                        f"{path.relative_to(_REPO_ROOT)}:{line_number} {target_path}"
                    )

    assert offenders == []


def test_importing_phospy_does_not_emit_pandas_copy_on_write_warning() -> None:
    env = os.environ.copy()
    python_path = str(_REPO_ROOT / "src")
    if env.get("PYTHONPATH"):
        python_path = f"{python_path}{os.pathsep}{env['PYTHONPATH']}"
    env["PYTHONPATH"] = python_path
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import warnings; import pandas as pd; "
                "warning_type = getattr(pd.errors, 'Pandas4Warning', Warning); "
                "warnings.simplefilter('error', warning_type); "
                "import phospy"
            ),
        ],
        cwd=_REPO_ROOT,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr


@pytest.mark.skipif(
    not _PANDAS_COPY_ON_WRITE_OPTION_IS_MUTABLE,
    reason="pandas >=3 deprecates mode.copy_on_write and always enables CoW",
)
def test_importing_phospy_preserves_mutable_pandas_copy_on_write_option() -> None:
    env = os.environ.copy()
    python_path = str(_REPO_ROOT / "src")
    if env.get("PYTHONPATH"):
        python_path = f"{python_path}{os.pathsep}{env['PYTHONPATH']}"
    env["PYTHONPATH"] = python_path
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import pandas as pd; "
                "pd.set_option('mode.copy_on_write', False); "
                "import phospy; "
                "raise SystemExit(0 if pd.get_option('mode.copy_on_write') "
                "is False else 1)"
            ),
        ],
        cwd=_REPO_ROOT,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr


@pytest.mark.skipif(
    not _PANDAS_COPY_ON_WRITE_OPTION_IS_MUTABLE,
    reason="pandas >=3 deprecates mode.copy_on_write and always enables CoW",
)
@pytest.mark.parametrize("copy_on_write", [False, True])
def test_phospy_frame_borrowing_preserves_pandas_copy_on_write_option(
    copy_on_write: bool,
) -> None:
    with pd.option_context("mode.copy_on_write", copy_on_write):
        dataset = trusted_analysis_ready_dataset_from_tables(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        )

        public_snapshot = dataset.to_dataframe()
        public_snapshot.iloc[0, 0] = 123.0
        borrowed = dataset._borrow_phospho_frame()
        _mutate_existing_borrowed_frame_cell(
            owner=dataset._phospho,
            borrowed=borrowed,
            value=999.0,
        )

        assert pd.get_option("mode.copy_on_write") == copy_on_write


@pytest.mark.skipif(
    _PANDAS_COPY_ON_WRITE_OPTION_IS_MUTABLE,
    reason="pandas <3 has a mutable copy_on_write option covered separately",
)
def test_phospy_frame_borrowing_does_not_touch_deprecated_copy_on_write_option() -> (
    None
):
    with warnings.catch_warnings():
        warnings.simplefilter("error", _pandas4_warning_type())
        dataset = trusted_analysis_ready_dataset_from_tables(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        )

        borrowed = dataset._borrow_phospho_frame()
        _mutate_existing_borrowed_frame_cell(
            owner=dataset._phospho,
            borrowed=borrowed,
            value=999.0,
        )


def test_representative_kinase_workflow_preserves_pandas_global_options() -> None:
    before = _pandas_global_option_snapshot()

    _kinase_result()

    assert _pandas_global_option_snapshot() == before


def test_concurrent_borrowed_access_preserves_pandas_copy_on_write_option() -> None:
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )

    def worker(value: float) -> float:
        borrowed = dataset._borrow_phospho_frame()
        try:
            borrowed.iloc[0, 0] = value
        except ValueError:
            pass
        return float(dataset._phospho.iloc[0, 0])

    if _PANDAS_COPY_ON_WRITE_OPTION_IS_MUTABLE:
        context = pd.option_context("mode.copy_on_write", False)
    else:
        context = warnings.catch_warnings()

    with context:
        if not _PANDAS_COPY_ON_WRITE_OPTION_IS_MUTABLE:
            warnings.simplefilter("error", _pandas4_warning_type())
        with ThreadPoolExecutor(max_workers=4) as executor:
            owner_values = tuple(executor.map(worker, range(16)))

        assert owner_values == (1.0,) * 16
        if _PANDAS_COPY_ON_WRITE_OPTION_IS_MUTABLE:
            assert pd.get_option("mode.copy_on_write") is False


def test_internal_borrowed_dataset_access_is_detached_snapshot_contract() -> None:
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=_mixed_numeric_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )

    with _count_dataframe_deep_copies() as counts:
        borrowed = dataset._borrow_phospho_frame()
        _mutate_existing_borrowed_frame_cell(
            owner=dataset._phospho,
            borrowed=borrowed,
            value=999.0,
        )
        _try_borrowed_column_write(borrowed)

    assert borrowed is not dataset._phospho
    assert counts.dataframe_deep == 1
    assert not hasattr(dataset, "borrow_phospho_frame")
    assert float(dataset._phospho.iloc[0, 0]) == 1.0
    assert "borrowed_only" not in dataset._phospho.columns


def test_dataset_internal_view_reuses_one_immutable_numeric_snapshot() -> None:
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=_mixed_numeric_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    view = DatasetInternalView(dataset)

    with _count_dataframe_deep_copies() as counts:
        first = view.phospho
        second = view.phospho
        _assert_numpy_blocks_readonly(first)
        _assert_numpy_blocks_readonly(second)
        _mutate_existing_borrowed_frame_cell(
            owner=dataset._phospho,
            borrowed=first,
            value=999.0,
        )
        second.loc[:, "borrowed_only"] = [1.0, 2.0]
        third = view.phospho

    assert counts.dataframe_deep == 1
    assert first is not second
    assert second is not third
    assert float(dataset._phospho.iloc[0, 0]) == 1.0
    assert "borrowed_only" not in third.columns


def test_dataset_internal_view_freezes_metadata_object_cells_once() -> None:
    original_payload = _mutable_object_payload()
    site_metadata = _object_payload_frame_from_site_metadata(original_payload)
    site_metadata.loc[:, "extension_label"] = pd.Series(
        ["alpha", "beta"],
        index=site_metadata.index,
        dtype="string",
    )
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=_phospho(),
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    view = DatasetInternalView(dataset)

    with _count_dataframe_deep_copies() as counts:
        first = view.site_metadata
        second = view.site_metadata

    row_label = first.index[0]
    first_payload = first.loc[row_label, _OBJECT_PAYLOAD_COLUMN]
    second_payload = second.loc[row_label, _OBJECT_PAYLOAD_COLUMN]
    _mutate_object_payload(original_payload, "external")

    assert counts.dataframe_deep == 1
    assert first_payload is second_payload
    assert _immutable_object_payload_state(first_payload) == _OBJECT_PAYLOAD_STATE
    assert _immutable_object_payload_state(second_payload) == _OBJECT_PAYLOAD_STATE
    _assert_numpy_blocks_readonly(first)


def test_dataset_internal_view_extension_dtype_metadata_is_owner_safe() -> None:
    site_metadata = _site_metadata()
    site_metadata.loc[:, "extension_label"] = pd.Series(
        ["alpha", "beta"],
        index=site_metadata.index,
        dtype="string",
    )
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=_phospho(),
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    view = DatasetInternalView(dataset)

    first = view.site_metadata
    first.loc[first.index[0], "extension_label"] = "changed"
    second = view.site_metadata

    assert str(second.loc[second.index[0], "extension_label"]) == "alpha"
    assert (
        str(
            dataset._site_metadata.loc[
                dataset._site_metadata.index[0], "extension_label"
            ]
        )
        == "alpha"
    )


def test_independent_dataset_internal_views_do_not_share_snapshot_cache() -> None:
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=_mixed_numeric_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )

    with _count_dataframe_deep_copies() as counts:
        first_run = DatasetInternalView(dataset).phospho
        second_run = DatasetInternalView(dataset).phospho

    assert counts.dataframe_deep == 2
    assert first_run is not second_run
    _assert_numpy_blocks_readonly(first_run)
    _assert_numpy_blocks_readonly(second_run)


def test_representative_differential_workflow_bounds_full_matrix_copies() -> None:
    request = _differential_workflow_request(n_sites=8)

    with _count_full_matrix_deep_copies(
        shape=request.dataset._phospho.shape,
        columns=tuple(request.dataset._phospho.columns),
    ) as counts:
        result = DifferentialAnalysisWorkflow().run(request)

    assert isinstance(result, DifferentialAnalysisResult)
    assert counts.full_matrix_deep <= 3


@pytest.mark.parametrize(
    ("surface_name", "mutator"),
    [
        pytest.param("slice-values", _mutate_numeric_slice_values, id="slice-values"),
        pytest.param("values", _mutate_numeric_values, id="values"),
        pytest.param(
            "to-numpy-copy-false",
            _mutate_numeric_to_numpy_copy_false,
            id="to-numpy-copy-false",
        ),
        pytest.param(
            "blocks-and-bases",
            _mutate_numeric_blocks_and_bases,
            id="blocks-and-bases",
        ),
    ],
)
@pytest.mark.parametrize(
    ("access_name", "accessor"),
    [
        pytest.param(
            "internal-borrow",
            lambda dataset: dataset._borrow_phospho_frame(),
            id="internal-borrow",
        ),
        pytest.param(
            "public-export",
            lambda dataset: dataset.to_dataframe(),
            id="public-export",
        ),
    ],
)
def test_restoring_writeability_on_dataset_numeric_surfaces_cannot_mutate_owner(
    access_name: str,
    accessor: Callable[[AnalysisReadyPhosphoDataset], pd.DataFrame],
    surface_name: str,
    mutator: Callable[[pd.DataFrame], bool],
) -> None:
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=_mixed_numeric_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )

    before = dataset._phospho.copy(deep=True)
    exposed = accessor(dataset)

    assert mutator(exposed), f"{access_name}:{surface_name}"
    pd.testing.assert_frame_equal(dataset._phospho, before)


def test_internal_borrowed_prediction_and_scoring_access_is_mutation_isolated() -> None:
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [0.2, 0.1],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
    profile_scores = pd.DataFrame(
        {"MAP2K6": [0.8, 0.2]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
    )
    rank_weighted = pd.DataFrame(
        {"MAP2K6": [0.75, 0.15]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
    )
    prediction_result = KinasePredictionResult._from_owned(pred_mat=pred_mat)
    scoring_result = KinaseScoringResult._from_owned(
        profile_scores=profile_scores,
        rank_weighted_fusion_scores=rank_weighted,
    )

    with _count_dataframe_deep_copies() as counts:
        borrowed_pred = prediction_result._borrow_pred_mat_frame()
        borrowed_profile = scoring_result._borrow_profile_scores_frame()
        borrowed_rank_weighted = (
            scoring_result._borrow_rank_weighted_fusion_scores_frame()
        )
        _mutate_existing_borrowed_frame_cell(
            owner=prediction_result._pred_mat,
            borrowed=borrowed_pred,
            value=99.0,
        )
        _mutate_existing_borrowed_frame_cell(
            owner=scoring_result._profile_scores,
            borrowed=borrowed_profile,
            value=88.0,
        )
        assert borrowed_rank_weighted is not None
        assert scoring_result._rank_weighted_fusion_scores is not None
        _mutate_existing_borrowed_frame_cell(
            owner=scoring_result._rank_weighted_fusion_scores,
            borrowed=borrowed_rank_weighted,
            value=77.0,
        )

    assert counts.dataframe_deep == 3
    assert borrowed_pred is not prediction_result._pred_mat
    assert borrowed_profile is not scoring_result._profile_scores
    assert borrowed_rank_weighted is not scoring_result._rank_weighted_fusion_scores
    assert float(prediction_result._pred_mat.iloc[0, 0]) == 0.9
    assert float(scoring_result._profile_scores.iloc[0, 0]) == 0.8
    assert scoring_result._rank_weighted_fusion_scores is not None
    assert float(scoring_result._rank_weighted_fusion_scores.iloc[0, 0]) == 0.75
    assert not hasattr(prediction_result, "borrow_pred_mat_frame")
    assert not hasattr(scoring_result, "borrow_profile_scores_frame")


def test_borrowed_extension_array_frame_falls_back_to_deep_copy_for_isolation() -> None:
    frame = pd.DataFrame({"label": pd.Series(["a", "b"], dtype="string")})

    with _count_dataframe_deep_copies() as counts:
        borrowed = _borrow_dataframe(frame)
        borrowed.iloc[0, 0] = "changed"

    assert str(frame.iloc[0, 0]) == "a"
    assert str(borrowed.iloc[0, 0]) == "changed"
    assert counts.dataframe_deep == 1


def test_borrowed_series_access_is_detached_and_mutation_isolated() -> None:
    series = pd.Series([1.0, 2.0], index=["a", "b"], name="values")

    with _count_dataframe_deep_copies() as counts:
        borrowed = _borrow_series(series)
        _mutate_existing_borrowed_series_value(
            owner=series,
            borrowed=borrowed,
            value=999.0,
        )

    assert borrowed is not series
    assert counts.dataframe_deep == 0
    assert float(series.iloc[0]) == 1.0
    series.iloc[0] = 321.0
    assert float(series.iloc[0]) == 321.0


def test_signalome_validator_read_path_does_not_mutate_internal_frames() -> None:
    request = _signalome_request_for_read_path_mutation_checks()
    dataset = request.kinase_result.dataset
    prediction_result = request.kinase_result.prediction_result
    scoring_result = request.kinase_result.scoring_result

    dataset_phospho_before = fingerprint_table(
        dataset._borrow_phospho_frame(),
        name="dataset.phospho",
    )
    dataset_site_metadata_before = fingerprint_table(
        dataset._borrow_site_metadata_frame(),
        name="dataset.site_metadata",
    )
    prediction_before = fingerprint_table(
        prediction_result._borrow_pred_mat_frame(),
        name="prediction_result.pred_mat",
    )
    score_before = fingerprint_table(
        scoring_result._borrow_profile_scores_frame(),
        name="scoring_result.profile_scores",
    )
    assert scoring_result._borrow_rank_weighted_fusion_scores_frame() is not None
    rank_weighted_before = fingerprint_table(
        scoring_result._borrow_rank_weighted_fusion_scores_frame(),
        name="scoring_result.rank_weighted_fusion_scores",
    )

    validated = SignalomeWorkflowValidator().run(request)
    assert validated is request

    dataset_phospho_after = fingerprint_table(
        dataset._borrow_phospho_frame(),
        name="dataset.phospho",
    )
    dataset_site_metadata_after = fingerprint_table(
        dataset._borrow_site_metadata_frame(),
        name="dataset.site_metadata",
    )
    prediction_after = fingerprint_table(
        prediction_result._borrow_pred_mat_frame(),
        name="prediction_result.pred_mat",
    )
    score_after = fingerprint_table(
        scoring_result._borrow_profile_scores_frame(),
        name="scoring_result.profile_scores",
    )
    assert scoring_result._borrow_rank_weighted_fusion_scores_frame() is not None
    rank_weighted_after = fingerprint_table(
        scoring_result._borrow_rank_weighted_fusion_scores_frame(),
        name="scoring_result.rank_weighted_fusion_scores",
    )

    assert (
        dataset_phospho_before.tolerance_hash_value
        == dataset_phospho_after.tolerance_hash_value
    )
    assert (
        dataset_site_metadata_before.tolerance_hash_value
        == dataset_site_metadata_after.tolerance_hash_value
    )
    assert (
        prediction_before.tolerance_hash_value == prediction_after.tolerance_hash_value
    )
    assert score_before.tolerance_hash_value == score_after.tolerance_hash_value
    assert (
        rank_weighted_before.tolerance_hash_value
        == rank_weighted_after.tolerance_hash_value
    )


def test_signalome_interpreter_read_path_does_not_mutate_dataset_frames() -> None:
    request = _signalome_request_for_read_path_mutation_checks()
    dataset = request.kinase_result.dataset

    phospho_before = fingerprint_table(
        dataset._borrow_phospho_frame(),
        name="dataset.phospho",
    )
    site_metadata_before = fingerprint_table(
        dataset._borrow_site_metadata_frame(),
        name="dataset.site_metadata",
    )

    SignalomeWorkflowInterpreter().run(request)

    phospho_after = fingerprint_table(
        dataset._borrow_phospho_frame(),
        name="dataset.phospho",
    )
    site_metadata_after = fingerprint_table(
        dataset._borrow_site_metadata_frame(),
        name="dataset.site_metadata",
    )

    assert phospho_before.tolerance_hash_value == phospho_after.tolerance_hash_value
    assert (
        site_metadata_before.tolerance_hash_value
        == site_metadata_after.tolerance_hash_value
    )


def test_owned_construction_frames_can_be_mutated_after_owned_transfer() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    provenance_source = trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    assert provenance_source.provenance is not None
    builder_provenance = replace(
        provenance_source.provenance,
        workflow_name="dataset_builder",
        workflow_parameters={
            "construction": {
                "method": "AnalysisReadyDatasetBuilder.run",
                "processing_state_establishment": {
                    "source": "test builder-shaped provenance"
                },
            }
        },
    )

    dataset = AnalysisReadyPhosphoDataset._from_builder_output(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
        provenance=builder_provenance,
    )

    phospho.iloc[0, 0] = 321.0
    site_metadata.iloc[0, 0] = "UPDATED_GENE"
    assert float(dataset._borrow_phospho_frame().iloc[0, 0]) == 321.0
    assert str(dataset._borrow_site_metadata_frame().iloc[0, 0]) == "UPDATED_GENE"

    public_snapshot = dataset.phospho
    public_snapshot.iloc[0, 0] = 123.0
    assert float(dataset._borrow_phospho_frame().iloc[0, 0]) == 321.0


def test_internal_borrowed_accessors_are_not_public_api_exports() -> None:
    import phospy
    import phospy.frames.ownership as frame_ownership
    import phospy.science.datasets as datasets

    assert not any("borrow" in name for name in phospy.__all__)
    assert not any("borrow" in name for name in datasets.__all__)
    assert not any("borrow" in name for name in frame_ownership.__all__)
    assert not hasattr(AnalysisReadyPhosphoDataset, "borrow_phospho_frame")
    assert not hasattr(AnalysisReadyPhosphoDataset, "borrow_site_metadata_frame")
    assert not hasattr(DatasetPreprocessingReport, "borrow_row_counts_frame")


def test_frame_ownership_helpers_are_owned_by_frames_ownership_module() -> None:
    assert own_dataframe.__module__ == "phospy.frames.ownership"
    assert own_optional_dataframe.__module__ == "phospy.frames.ownership"
    assert export_dataframe.__module__ == "phospy.frames.ownership"
    assert export_optional_dataframe.__module__ == "phospy.frames.ownership"
    assert own_series.__module__ == "phospy.frames.ownership"
    assert own_optional_series.__module__ == "phospy.frames.ownership"
    assert export_series.__module__ == "phospy.frames.ownership"
    assert export_optional_series.__module__ == "phospy.frames.ownership"
    assert _borrow_dataframe.__module__ == "phospy.frames.ownership"
    assert _borrow_series.__module__ == "phospy.frames.ownership"


def test_root_frame_ownership_module_is_removed() -> None:
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("phospy._frame_ownership")


def test_frame_ownership_helper_constructor_and_export_behaviour_is_preserved() -> None:
    original = _phospho()
    owned = own_dataframe(original, field_name="x")
    original.iloc[0, 0] = 999.0
    assert float(owned.iloc[0, 0]) == 1.0

    exported = export_dataframe(owned)
    exported.iloc[0, 0] = 111.0
    assert float(owned.iloc[0, 0]) == 1.0


def test_frame_ownership_helper_policy_matrix_documents_export_and_borrow_modes() -> (
    None
):
    frame_owner = _phospho()
    frame_cases: tuple[
        tuple[str, Callable[[pd.DataFrame], pd.DataFrame], bool],
        ...,
    ] = (
        ("safe_public_copy", export_dataframe, True),
        ("borrowed_internal_snapshot", _borrow_dataframe, False),
    )

    for category, accessor, writable_export in frame_cases:
        snapshot = accessor(frame_owner)
        assert snapshot is not frame_owner
        try:
            snapshot.iloc[0, 0] = 222.0
        except (TypeError, ValueError):
            assert not writable_export, category
        else:
            if writable_export:
                assert float(snapshot.iloc[0, 0]) == 222.0
        assert float(frame_owner.iloc[0, 0]) == 1.0

    series_owner = pd.Series([1.0, 2.0], index=["a", "b"], name="values")
    series_cases: tuple[
        tuple[str, Callable[[pd.Series], pd.Series], bool],
        ...,
    ] = (
        ("safe_public_copy", export_series, True),
        ("borrowed_internal_snapshot", _borrow_series, False),
    )

    for category, accessor, writable_export in series_cases:
        snapshot = accessor(series_owner)
        assert snapshot is not series_owner
        try:
            snapshot.iloc[0] = 222.0
        except (TypeError, ValueError):
            assert not writable_export, category
        else:
            if writable_export:
                assert float(snapshot.iloc[0]) == 222.0
        assert float(series_owner.iloc[0]) == 1.0


def test_frame_ownership_optional_helpers_return_none_for_none() -> None:
    assert own_optional_dataframe(None, field_name="x") is None
    assert export_optional_dataframe(None) is None
    assert own_optional_series(None, field_name="x") is None
    assert export_optional_series(None) is None


def test_frame_ownership_series_helper_copy_behaviour_is_preserved() -> None:
    original = pd.Series([1.0, 2.0], index=["a", "b"], name="values")
    owned = own_series(original, field_name="x")
    original.iloc[0] = 999.0
    assert float(owned.iloc[0]) == 1.0

    exported = export_series(owned)
    exported.iloc[0] = 111.0
    assert float(owned.iloc[0]) == 1.0


def test_prediction_result_public_export_copy_default_is_safe() -> None:
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [0.2, 0.1],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
    result = KinasePredictionResult._from_owned(pred_mat=pred_mat)

    exported = result.to_dataframe()
    exported.iloc[0, 0] = 0.0

    assert float(result.pred_mat.iloc[0, 0]) == 0.9


def test_safe_public_export_does_not_change_owned_provenance_state() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    fingerprint_before = fingerprint_table(built.phospho, name="dataset.phospho")

    safe_copy = built.to_dataframe()
    safe_copy.iloc[0, 0] = 999.0

    fingerprint_after = fingerprint_table(built.phospho, name="dataset.phospho")
    assert (
        fingerprint_before.tolerance_hash_value
        == fingerprint_after.tolerance_hash_value
    )


@pytest.mark.parametrize(
    ("export_factory",),
    [
        pytest.param(
            lambda: (
                trusted_analysis_ready_dataset_from_tables(
                    phospho=_phospho(),
                    site_metadata=_site_metadata(),
                    organism=Organism.RAT,
                    intensity_scale_state=supported_linear_intensity_scale_state(
                        has_total_matrix=False
                    ),
                    processing_state=supported_linear_processing_state(
                        has_total_matrix=False
                    ),
                ).to_dataframe
            ),
            id="dataset-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: (
                AnalysisReadyDatasetBuilder()
                .run(
                    DatasetBuildRequest(
                        phospho=_phospho(),
                        site_metadata=_site_metadata(),
                        organism=Organism.RAT,
                        input_intensity_scale="linear",
                    )
                )
                .to_dataframe
            ),
            id="builder-output-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: (
                KinasePredictionResult._from_owned(
                    pred_mat=pd.DataFrame(
                        {
                            "MAP2K6": [0.9, 0.8],
                            "AKT1": [0.2, 0.1],
                        },
                        index=["MAPK14;Y182;", "GSK3B;S9;"],
                    )
                ).to_dataframe
            ),
            id="prediction-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: (
                SignalomeWorkflow()
                .run(
                    SignalomeWorkflowRequest(
                        kinase_result=_kinase_result(),
                        config=build_signalome_config(
                            substrate_support_cutoff=0.5,
                            reference_context_compatibility_policy=(
                                _ALLOW_UNKNOWN_REFERENCE_CONTEXT
                            ),
                        ),
                    )
                )
                .to_dataframe
            ),
            id="signalome-result-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: (
                SignalomeWorkflow()
                .run(
                    SignalomeWorkflowRequest(
                        kinase_result=_kinase_result(),
                        config=build_signalome_config(
                            substrate_support_cutoff=0.5,
                            reference_context_compatibility_policy=(
                                _ALLOW_UNKNOWN_REFERENCE_CONTEXT
                            ),
                        ),
                    )
                )
                .module_assignments.to_pandas
            ),
            id="signalome-module-assignments-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: (
                SignalomeWorkflow()
                .run(
                    SignalomeWorkflowRequest(
                        kinase_result=_kinase_result(),
                        config=build_signalome_config(
                            substrate_support_cutoff=0.5,
                            reference_context_compatibility_policy=(
                                _ALLOW_UNKNOWN_REFERENCE_CONTEXT
                            ),
                        ),
                    )
                )
                .signalome_modules.to_pandas
            ),
            id="signalome-modules-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: (
                SignalomeWorkflow()
                .run(
                    SignalomeWorkflowRequest(
                        kinase_result=_kinase_result(),
                        config=build_signalome_config(
                            substrate_support_cutoff=0.5,
                            reference_context_compatibility_policy=(
                                _ALLOW_UNKNOWN_REFERENCE_CONTEXT
                            ),
                        ),
                    )
                )
                .kinase_network.to_pandas
            ),
            id="signalome-network-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: PhosphoIntensityMatrix(frame=_phospho()).to_pandas,
            id="table-schema-export-rejects-copy-keyword",
        ),
    ],
)
def test_public_export_copy_keyword_rejection_contract_matrix(
    export_factory: Callable[[], Callable[..., object]],
) -> None:
    _assert_copy_keyword_rejected(export_factory())


def test_public_signalome_exports_isolated_from_mutation() -> None:
    result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=_kinase_result(),
            config=build_signalome_config(
                substrate_support_cutoff=0.5,
                network_min_paired_finite_observations=3,
                reference_context_compatibility_policy=(
                    _ALLOW_UNKNOWN_REFERENCE_CONTEXT
                ),
            ),
        )
    )

    # Public exports must be defensive snapshots across representative result types.
    for getter in (
        result.module_assignments.to_pandas,
        result.signalome_modules.to_pandas,
        result.kinase_network.to_pandas,
    ):
        _assert_dataframe_getter_defensive_snapshot(getter)


@pytest.mark.parametrize(
    ("getter_factory",),
    [
        pytest.param(
            lambda: (
                trusted_analysis_ready_dataset_from_tables(
                    phospho=_phospho(),
                    site_metadata=_site_metadata(),
                    organism=Organism.RAT,
                    intensity_scale_state=supported_linear_intensity_scale_state(
                        has_total_matrix=False
                    ),
                    processing_state=supported_linear_processing_state(
                        has_total_matrix=False
                    ),
                ).to_dataframe
            ),
            id="dataset-export-snapshot",
        ),
        pytest.param(
            lambda: (
                KinasePredictionResult._from_owned(
                    pred_mat=pd.DataFrame(
                        {
                            "MAP2K6": [0.9, 0.8],
                            "AKT1": [0.2, 0.1],
                        },
                        index=["MAPK14;Y182;", "GSK3B;S9;"],
                    )
                ).to_dataframe
            ),
            id="prediction-export-snapshot",
        ),
        pytest.param(
            lambda: _references().kinase_substrate_map_dataframe,
            id="reference-export-snapshot",
        ),
    ],
)
def test_public_export_snapshot_contract_matrix(
    getter_factory: Callable[[], Callable[[], pd.DataFrame]],
) -> None:
    _assert_dataframe_getter_defensive_snapshot(getter_factory())


def test_dataset_dataframe_properties_are_defensive_snapshots() -> None:
    sample_metadata = pd.DataFrame(
        {"batch": [1, 2]},
        index=pd.Index(["sample_a", "sample_b"]),
    )
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        sample_metadata=sample_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )

    _assert_dataframe_getter_defensive_snapshot(lambda: dataset.phospho)
    _assert_dataframe_getter_defensive_snapshot(lambda: dataset.site_metadata)
    _assert_optional_dataframe_getter_defensive_snapshot(
        lambda: dataset.sample_metadata
    )


def test_preprocessing_report_dataframe_properties_are_defensive_snapshots() -> None:
    report = DatasetPreprocessingReport.from_rows(
        row_count_rows=(
            PreprocessingRowCountRow(
                stage="missing_data",
                input_rows=2,
                output_rows=2,
                dropped_rows=0,
            ),
        ),
        operation_rows=(
            PreprocessingOperationRow(
                step_order=1,
                stage="missing_data",
                operation="forbid",
                parameters={},
                input_rows=2,
                output_rows=2,
                notes=None,
            ),
        ),
        row_audit_rows=(
            PreprocessingRowAuditRow(
                stage="missing_data",
                action="retained",
                reason="complete",
                source_row_id="MAPK14;Y182;",
                site_id="MAPK14;Y182;",
                retained=True,
                retained_row_id="MAPK14;Y182;",
                source_rows=None,
                retained_row=None,
                parameter_snapshot={},
            ),
        ),
        duplicate_site_resolution_rows=(
            DuplicateSiteResolutionRow(
                site_key="MAPK14;Y182;",
                display_id="MAPK14;Y182;",
                site_id="MAPK14;Y182;",
                source_row_id="MAPK14;Y182;",
                retained=True,
                resolution_policy="max_mean_signal",
                aggregation_method="max_mean_signal",
                missing_value_policy=None,
                metadata_resolution_policy="retain_row_ranked_by_observed_values_then_mean_signal_then_input_order",
                retained_reason=None,
                dropped_reason=None,
                observed_values=None,
                mean_signal=1.0,
                n_source_rows=1,
                n_aggregated_rows=1,
                source_protein_id="MAPK14",
                source_gene_symbol="MAPK14",
                source_site="Y182",
                source_site_sequence="A",
                metadata_conflict_detected=False,
            ),
        ),
        metadata_conflict_rows=(
            MetadataConflictRow(
                site_key="MAPK14;Y182;",
                display_id="MAPK14;Y182;",
                site_id="MAPK14;Y182;",
                field="protein_id",
                values=["MAPK14"],
                n_distinct_values=1,
                source_row_ids=["MAPK14;Y182;"],
            ),
        ),
        comparison_group_stats_rows=(
            ComparisonGroupStatsRow(
                site_id="MAPK14;Y182;",
                group="group_a",
                n=1,
                mean=1.0,
                sd=None,
                sem=None,
                median=1.0,
                min=1.0,
                max=1.0,
                sample_ids=["sample_a"],
            ),
        ),
        comparison_pair_stats_rows=(
            ComparisonPairStatsRow(
                site_id="MAPK14;Y182;",
                comparison="p_group_a_group_b",
                left_group="group_a",
                right_group="group_b",
                left_n=1,
                right_n=1,
                left_mean=1.0,
                right_mean=2.0,
                left_sd=None,
                right_sd=None,
                left_sem=None,
                right_sem=None,
                effect_size=-1.0,
                left_median=1.0,
                right_median=2.0,
                left_min=1.0,
                right_min=2.0,
                left_max=1.0,
                right_max=2.0,
            ),
        ),
    )

    for getter in (
        lambda: report.row_counts,
        lambda: report.operations,
        lambda: report.row_audit,
    ):
        _assert_dataframe_getter_defensive_snapshot(getter)

    for optional_getter in (
        lambda: report.duplicate_site_resolution,
        lambda: report.metadata_conflicts,
        lambda: report.comparison_group_stats,
        lambda: report.comparison_pair_stats,
    ):
        _assert_optional_dataframe_getter_defensive_snapshot(optional_getter)


def test_kinase_result_table_properties_are_defensive_snapshots() -> None:
    scoring_result = KinaseScoringResult._from_owned(
        profile_scores=pd.DataFrame(
            {"MAP2K6": [0.8, 0.2]},
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
        ),
        motif_scores=pd.DataFrame(
            {"MAP2K6": [0.7, 0.1]},
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
        ),
        rank_weighted_fusion_scores=pd.DataFrame(
            {"MAP2K6": [0.75, 0.15]},
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
        ),
        score_fusion_weights=pd.DataFrame(
            {"MAP2K6": [1.0, 1.0]},
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
        ),
    )
    _assert_dataframe_getter_defensive_snapshot(lambda: scoring_result.profile_scores)
    for optional_getter in (
        lambda: scoring_result.motif_scores,
        lambda: scoring_result.rank_weighted_fusion_scores,
        lambda: scoring_result.score_fusion_weights,
    ):
        _assert_optional_dataframe_getter_defensive_snapshot(optional_getter)

    prediction_result = KinasePredictionResult._from_owned(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.9, 0.8]},
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
        ),
        substrate_list=pd.DataFrame(
            {
                "kinase": ["MAP2K6"],
                "substrate_site": ["MAPK14;Y182;"],
                "score": [0.9],
                "rank": [1],
            }
        ),
    )
    _assert_dataframe_getter_defensive_snapshot(lambda: prediction_result.pred_mat)
    _assert_optional_dataframe_getter_defensive_snapshot(
        lambda: prediction_result.substrate_list
    )

    contribution_table = pd.DataFrame.from_records(
        [
            {
                "kinase": "MAP2K6",
                "substrate_site": "MAPK14;Y182;",
                "substrate_identifier": "MAPK14;Y182;",
                "value_used_in_scoring": 0.8,
                "score_component": "rank_weighted_fusion_scores",
                "score_source": "profile_only_motif_missing_or_constant",
                "reference_source_name": "fixture",
                "reference_source_version": "v1",
                "reference_bundle_id": "fixture_bundle",
                "reference_identifier_namespace": "display_id",
                "status": "included",
                "exclusion_reason": None,
                "ambiguous": False,
            }
        ],
        columns=pd.Index(KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS),
    )
    workflow_result = KinaseWorkflowResult(
        dataset=trusted_analysis_ready_dataset_from_tables(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        ),
        references=_references(),
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        substrate_contributions=contribution_table,
    )
    _assert_optional_dataframe_getter_defensive_snapshot(
        lambda: workflow_result.substrate_contributions
    )

    activity_matrix = pd.DataFrame(
        {"MAP2K6": [1.0, 2.0]},
        index=pd.Index(["sample_a", "sample_b"]),
    )
    activity_input = ActivityInputMatrix.sample_level_abundance(activity_matrix)
    activity_result = KinaseActivityResult._from_owned(
        activity_matrix=activity_matrix,
        thresholded_substrate_mean_activity=pd.DataFrame(
            {"MAP2K6": [0.5, 1.5]},
            index=pd.Index(["sample_a", "sample_b"]),
        ),
        thresholded_substrate_counts=pd.Series(
            [2, 2],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_substrates",
        ),
        target_counts=pd.Series(
            [1, 1],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_targets",
        ),
        target_table=pd.DataFrame(
            {
                "site_id": ["MAPK14;Y182;"],
                "kinase": ["MAP2K6"],
                "score": [0.9],
            }
        ),
        input_semantics=activity_input.semantics,
        profile_metadata=activity_input.profile_metadata,
    )
    for getter in (
        lambda: activity_result.activity_matrix,
        lambda: activity_result.thresholded_substrate_mean_activity,
        lambda: activity_result.target_table,
    ):
        _assert_dataframe_getter_defensive_snapshot(getter)


def test_kinase_activity_result_series_properties_are_defensive_snapshots() -> None:
    activity_matrix = pd.DataFrame(
        {"MAP2K6": [1.0, 2.0]},
        index=pd.Index(["sample_a", "sample_b"]),
    )
    activity_input = ActivityInputMatrix.sample_level_abundance(activity_matrix)
    activity_result = KinaseActivityResult._from_owned(
        activity_matrix=activity_matrix,
        thresholded_substrate_mean_activity=pd.DataFrame(
            {"MAP2K6": [0.5, 1.5]},
            index=pd.Index(["sample_a", "sample_b"]),
        ),
        thresholded_substrate_counts=pd.Series(
            [2, 2],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_substrates",
        ),
        target_counts=pd.Series(
            [1, 1],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_targets",
        ),
        target_table=pd.DataFrame(
            {
                "site_id": ["MAPK14;Y182;"],
                "kinase": ["MAP2K6"],
                "score": [0.9],
            }
        ),
        input_semantics=activity_input.semantics,
        profile_metadata=activity_input.profile_metadata,
    )

    assert hasattr(activity_result, "thresholded_substrate_counts")
    assert hasattr(activity_result, "target_counts")

    thresholded_before = fingerprint_table(
        activity_result.thresholded_substrate_counts.to_frame(name="n_substrates"),
        name="outputs.activity.thresholded_substrate_counts",
    )
    exported_thresholded = activity_result.thresholded_substrate_counts
    exported_thresholded.iloc[0] = 999
    reread_thresholded = activity_result.thresholded_substrate_counts
    assert exported_thresholded is not reread_thresholded
    assert reread_thresholded.to_dict() == {"MAP2K6": 2, "AKT1": 2}
    thresholded_after = fingerprint_table(
        activity_result.thresholded_substrate_counts.to_frame(name="n_substrates"),
        name="outputs.activity.thresholded_substrate_counts",
    )
    assert (
        thresholded_before.tolerance_hash_value
        == thresholded_after.tolerance_hash_value
    )

    target_before = fingerprint_table(
        activity_result.target_counts.to_frame(name="n_targets"),
        name="outputs.activity.target_counts",
    )
    exported_target = activity_result.target_counts
    exported_target.iloc[0] = 999
    reread_target = activity_result.target_counts
    assert exported_target is not reread_target
    assert reread_target.to_dict() == {"MAP2K6": 1, "AKT1": 1}
    target_after = fingerprint_table(
        activity_result.target_counts.to_frame(name="n_targets"),
        name="outputs.activity.target_counts",
    )
    assert target_before.tolerance_hash_value == target_after.tolerance_hash_value


def test_signalome_result_table_properties_are_defensive_snapshots() -> None:
    signalome_result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=_kinase_result(),
            config=build_signalome_config(
                substrate_support_cutoff=0.5,
                network_min_paired_finite_observations=3,
                reference_context_compatibility_policy=(
                    _ALLOW_UNKNOWN_REFERENCE_CONTEXT
                ),
            ),
        )
    )

    for getter in (
        lambda: signalome_result.module_assignments.table,
        lambda: signalome_result.signalome_modules.table,
        lambda: signalome_result.kinase_network.edges,
    ):
        _assert_dataframe_getter_defensive_snapshot(getter)

    for optional_getter in (
        lambda: signalome_result.kinase_network.nodes,
        lambda: signalome_result.kinase_network.candidate_correlations,
        lambda: signalome_result.expanded_signalome,
        lambda: signalome_result.site_membership,
        lambda: signalome_result.protein_site_context,
    ):
        _assert_optional_dataframe_getter_defensive_snapshot_when_present(
            optional_getter
        )

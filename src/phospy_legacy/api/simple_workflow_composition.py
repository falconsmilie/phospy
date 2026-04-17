from __future__ import annotations

import inspect
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable

import pandas as pd

from ..activities.analysis import KinaseActivityAnalyzer
from ..activities.results import KinaseActivityResult
from ..datasets.models import AnalysisReadyPhosphoDataset
from ..datasets.schema import DatasetSchema
from ..internal.constants import ComparisonSpec
from ..internal.defaults import DEFAULT_MOTIF_FLANK_SIZE
from ..internal.types import PREDICTION_SVM_MODE_DEFAULT, PredictionSvmMode
from ..prediction.engines import KinaseWorkflowExecutionResult, KinaseWorkflowExecutor
from ..prediction.results import PredMatResult
from ..preprocessing.core import CorePreprocessingConfig
from ..preprocessing.modes import AnalysisReadyDatasetBuilder
from ..references import BundledReferenceProvider, ReferenceBundle
from ..validation.requests.workflow import WorkflowInputs

__all__ = [
    "ActivityAnalyzerProtocol",
    "AnalysisReadyBuilderProtocol",
    "ReferenceProviderProtocol",
    "SimpleKinaseExecutionGraph",
    "WorkflowExecutorProtocol",
    "create_default_simple_kinase_execution_graph",
]


@runtime_checkable
class AnalysisReadyBuilderProtocol(Protocol):
    """Protocol for analysis-ready preprocessing collaborators."""

    def build(
        self,
        *,
        phospho: pd.DataFrame | str | Path,
        preprocessing_config: CorePreprocessingConfig,
        total: pd.DataFrame | str | Path | None = None,
        phospho_encoding: str | None = None,
        schema: DatasetSchema | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
        source: str = "analysis ready dataset builder",
        phospho_only_source: str = "analysis ready dataset builder (phospho only)",
    ) -> AnalysisReadyPhosphoDataset: ...


@runtime_checkable
class ReferenceProviderProtocol(Protocol):
    """Protocol for reference-bundle resolution collaborators."""

    def resolve(
        self,
        *,
        species: str,
        reference: str = "auto",
    ) -> ReferenceBundle: ...


@runtime_checkable
class ActivityAnalyzerProtocol(Protocol):
    """Protocol for kinase activity analysis collaborators."""

    def run(
        self,
        *,
        pred_mat: pd.DataFrame | PredMatResult,
        phospho_matrix: pd.DataFrame,
        threshold: float = 0.0,
        min_substrates: int = 1,
        top_n_substrates: int = 1,
    ) -> KinaseActivityResult: ...


@runtime_checkable
class WorkflowExecutorProtocol(Protocol):
    """Protocol for validated workflow execution collaborators."""

    def validate_request(
        self,
        *,
        phospho_matrix: pd.DataFrame,
        site_sequences: Mapping[str, str] | pd.Series | None = None,
        reference_bundle: ReferenceBundle | None = None,
        prediction_config: object | None = None,
    ) -> WorkflowInputs: ...

    def execute_validated_request(
        self,
        request: WorkflowInputs,
    ) -> KinaseWorkflowExecutionResult: ...


def _resolve_callable(
    *,
    collaborator_name: str,
    collaborator: object,
    method_name: str,
) -> object:
    method = getattr(collaborator, method_name, None)
    if callable(method):
        return method
    msg = (
        f"SimpleKinaseExecutionGraph collaborator '{collaborator_name}' must define "
        f"a callable '{method_name}(...)' method."
    )
    raise TypeError(msg)


def _validate_keyword_invocation_shape(
    *,
    collaborator_name: str,
    collaborator: object,
    method_name: str,
    call_kwargs: tuple[str, ...],
) -> None:
    method = _resolve_callable(
        collaborator_name=collaborator_name,
        collaborator=collaborator,
        method_name=method_name,
    )
    signature = inspect.signature(method)
    placeholder_kwargs = {name: object() for name in call_kwargs}
    try:
        signature.bind(**placeholder_kwargs)
    except TypeError as error:
        joined_kwargs = ", ".join(call_kwargs)
        msg = (
            f"SimpleKinaseExecutionGraph collaborator '{collaborator_name}' has an "
            f"incompatible '{method_name}(...)' signature. Expected invocation with "
            f"keyword arguments: {joined_kwargs}. ({error})"
        )
        raise TypeError(msg) from error


def _validate_positional_invocation_shape(
    *,
    collaborator_name: str,
    collaborator: object,
    method_name: str,
    positional_count: int,
) -> None:
    method = _resolve_callable(
        collaborator_name=collaborator_name,
        collaborator=collaborator,
        method_name=method_name,
    )
    signature = inspect.signature(method)
    placeholder_args = tuple(object() for _ in range(positional_count))
    try:
        signature.bind(*placeholder_args)
    except TypeError as error:
        msg = (
            f"SimpleKinaseExecutionGraph collaborator '{collaborator_name}' has an "
            f"incompatible '{method_name}(...)' signature. Expected invocation with "
            f"{positional_count} positional argument(s). ({error})"
        )
        raise TypeError(msg) from error


@dataclass(frozen=True, slots=True)
class SimpleKinaseExecutionGraph:
    """Advanced collaborator graph for `SimpleKinaseWorkflow`.

    This is an extension seam for advanced users who need to replace default
    workflow collaborators as a unit. The graph is treated as a stable public
    composition contract, while concrete collaborator implementations may evolve
    independently.
    """

    analysis_ready_builder: AnalysisReadyBuilderProtocol
    reference_provider: ReferenceProviderProtocol
    activity_analyzer: ActivityAnalyzerProtocol
    workflow_executor: WorkflowExecutorProtocol

    def __post_init__(self) -> None:
        required_collaborators = (
            (
                "analysis_ready_builder",
                self.analysis_ready_builder,
                AnalysisReadyBuilderProtocol,
            ),
            (
                "reference_provider",
                self.reference_provider,
                ReferenceProviderProtocol,
            ),
            ("activity_analyzer", self.activity_analyzer, ActivityAnalyzerProtocol),
            ("workflow_executor", self.workflow_executor, WorkflowExecutorProtocol),
        )
        for collaborator_name, collaborator, protocol in required_collaborators:
            if collaborator is None:
                msg = (
                    f"SimpleKinaseExecutionGraph collaborator "
                    f"'{collaborator_name}' cannot be None."
                )
                raise TypeError(msg)
            if not isinstance(collaborator, protocol):
                msg = (
                    f"SimpleKinaseExecutionGraph collaborator '{collaborator_name}' "
                    f"must satisfy {protocol.__name__}."
                )
                raise TypeError(msg)

        _validate_keyword_invocation_shape(
            collaborator_name="analysis_ready_builder",
            collaborator=self.analysis_ready_builder,
            method_name="build",
            call_kwargs=(
                "phospho",
                "total",
                "phospho_encoding",
                "schema",
                "comparisons",
                "preprocessing_config",
                "source",
                "phospho_only_source",
            ),
        )
        _validate_keyword_invocation_shape(
            collaborator_name="reference_provider",
            collaborator=self.reference_provider,
            method_name="resolve",
            call_kwargs=("species", "reference"),
        )
        _validate_keyword_invocation_shape(
            collaborator_name="activity_analyzer",
            collaborator=self.activity_analyzer,
            method_name="run",
            call_kwargs=(
                "pred_mat",
                "phospho_matrix",
                "threshold",
                "min_substrates",
                "top_n_substrates",
            ),
        )
        _validate_keyword_invocation_shape(
            collaborator_name="workflow_executor",
            collaborator=self.workflow_executor,
            method_name="validate_request",
            call_kwargs=(
                "phospho_matrix",
                "site_sequences",
                "reference_bundle",
                "prediction_config",
            ),
        )
        _validate_positional_invocation_shape(
            collaborator_name="workflow_executor",
            collaborator=self.workflow_executor,
            method_name="execute_validated_request",
            positional_count=1,
        )


def create_default_simple_kinase_execution_graph(
    *,
    flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
    kernel: str = "rbf",
    svm_mode: PredictionSvmMode = PREDICTION_SVM_MODE_DEFAULT,
) -> SimpleKinaseExecutionGraph:
    """Build the default `SimpleKinaseExecutionGraph` used by public workflows."""

    return SimpleKinaseExecutionGraph(
        analysis_ready_builder=AnalysisReadyDatasetBuilder(),
        reference_provider=BundledReferenceProvider(),
        activity_analyzer=KinaseActivityAnalyzer(),
        workflow_executor=KinaseWorkflowExecutor(
            flank_size=flank_size,
            kernel=kernel,
            svm_mode=svm_mode,
        ),
    )

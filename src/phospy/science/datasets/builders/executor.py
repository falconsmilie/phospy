"""Internal executor for analysis-ready dataset construction."""

from __future__ import annotations

import inspect
from dataclasses import replace

import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.provenance.models import RunProvenance
from phospy.science.datasets.builders.contracts import (
    DatasetIntensityScaleResolverContract,
    DatasetPreprocessorContract,
    InterpretedDatasetBuildRequest,
    PreprocessedDatasetBuildTables,
)
from phospy.science.datasets.builders.preprocessing import (
    DatasetPreprocessor,
    DatasetProteinAwarePreparationRunner,
)
from phospy.science.datasets.builders.provenance_assembler import (
    DatasetRunProvenanceAssembler,
)
from phospy.science.datasets.builders.report_assembler import (
    DatasetPreprocessingReportAssembler,
)
from phospy.science.datasets.builders.site_sequence_boundary import (
    AnalysisReadySiteSequenceValidator,
)
from phospy.science.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.science.datasets.builders.transformation_state import (
    DatasetTransformationStateResolver,
    ResolvedDatasetTransformationState,
)
from phospy.science.datasets.models import (
    AnalysisReadyPhosphoDataset,
    DatasetPreprocessingReport,
)
from phospy.science.datasets.organism_coherence import (
    normalize_dataset_organism_state,
)
from phospy.science.datasets.preprocessing.models import (
    reject_external_corrected_output_after_downstream_preprocessing,
)
from phospy.science.datasets.preprocessing.protein_aware_preparation import (
    ProteinAwarePreparationResult,
)
from phospy.science.datasets.preprocessing.quantitative_scale_policy import (
    AdditivePreprocessingScaleGuard,
)
from phospy.science.datasets.preprocessing.stage_registry import (
    get_preprocessing_stage_metadata,
)
from phospy.science.transformations.contracts import Transformer
from phospy.science.transformations.quantitative_contracts import (
    QuantitativeContractState,
    initial_quantitative_contract_state,
)
from phospy.science.transformations.transformers import IdentityTransformer


class DatasetBuildExecutor:
    """Construct `AnalysisReadyPhosphoDataset` from interpreted builder input."""

    def __init__(
        self,
        *,
        transformer: Transformer | None = None,
        intensity_scale_resolver: DatasetIntensityScaleResolverContract | None = None,
        preprocessor: DatasetPreprocessorContract | None = None,
        site_sequence_validator: AnalysisReadySiteSequenceValidator | None = None,
        transformation_state_resolver: DatasetTransformationStateResolver | None = None,
        preprocessing_report_assembler: DatasetPreprocessingReportAssembler
        | None = None,
        protein_aware_preparation_runner: DatasetProteinAwarePreparationRunner
        | None = None,
        provenance_assembler: DatasetRunProvenanceAssembler | None = None,
        additive_preprocessing_scale_guard: AdditivePreprocessingScaleGuard
        | None = None,
    ) -> None:
        self._intensity_scale_resolver = (
            intensity_scale_resolver
            or DatasetIntensityScaleResolver(
                transformer=transformer or IdentityTransformer()
            )
        )
        self._preprocessor = preprocessor or DatasetPreprocessor()
        self._site_sequence_validator = (
            site_sequence_validator or AnalysisReadySiteSequenceValidator()
        )
        self._transformation_state_resolver = (
            transformation_state_resolver
            or DatasetTransformationStateResolver(
                intensity_scale_resolver=self._intensity_scale_resolver
            )
        )
        self._preprocessing_report_assembler = (
            preprocessing_report_assembler or DatasetPreprocessingReportAssembler()
        )
        self._protein_aware_preparation_runner = (
            protein_aware_preparation_runner or DatasetProteinAwarePreparationRunner()
        )
        self._provenance_assembler = (
            provenance_assembler or DatasetRunProvenanceAssembler()
        )
        self._additive_preprocessing_scale_guard = (
            additive_preprocessing_scale_guard or AdditivePreprocessingScaleGuard()
        )

    def run(
        self, request: InterpretedDatasetBuildRequest
    ) -> AnalysisReadyPhosphoDataset:
        self._validate_additive_preprocessing_scale_policy(request)
        preprocessed = self._run_preprocessing(request)
        validated_site_metadata = self._validate_analysis_ready_site_sequences(
            preprocessed
        )
        transformed = self._resolve_transformation_state(
            request=request,
            preprocessed=preprocessed,
            validated_site_metadata=validated_site_metadata,
        )
        normalized_organism_state = normalize_dataset_organism_state(
            phospho=transformed.phospho,
            site_metadata=validated_site_metadata,
            error_type=DatasetValidationError,
        )
        transformed = replace(
            transformed,
            phospho=normalized_organism_state.phospho,
        )
        validated_site_metadata = normalized_organism_state.site_metadata
        protein_aware_preparation = self._run_protein_aware_preparation(
            request=request,
            transformed=transformed,
            validated_site_metadata=validated_site_metadata,
        )
        report = self._assemble_preprocessing_report(
            request=request,
            preprocessed=preprocessed,
            transformed=transformed,
            protein_aware_preparation=protein_aware_preparation,
        )
        provenance = self._assemble_run_provenance(
            request=request,
            preprocessed=preprocessed,
            transformed=transformed,
            validated_site_metadata=validated_site_metadata,
            protein_aware_preparation=protein_aware_preparation,
        )
        return self._construct_dataset(
            request=request,
            preprocessed=preprocessed,
            transformed=transformed,
            validated_site_metadata=validated_site_metadata,
            report=report,
            protein_aware_preparation=protein_aware_preparation,
            provenance=provenance,
        )

    def _run_preprocessing(
        self,
        request: InterpretedDatasetBuildRequest,
    ) -> PreprocessedDatasetBuildTables:
        preprocessor_kwargs = {
            "phospho": request.phospho,
            "site_metadata": request.site_metadata,
            "sample_metadata": request.sample_metadata,
            "total": request.total,
            "plan": request.preprocessing_plan,
        }
        if _preprocessor_accepts_quantitative_contract_seed(self._preprocessor):
            preprocessor_kwargs["initial_quantitative_scale_kind"] = (
                request.declared_input_intensity_scale_kind
            )
            preprocessor_kwargs["initial_quantitative_meaning"] = (
                request.quantitative_meaning
            )
        if request.corrected_preprocessing_output is not None:
            preprocessor_kwargs["corrected_preprocessing_output"] = (
                request.corrected_preprocessing_output
            )
        return self._preprocessor.run(**preprocessor_kwargs)

    def _validate_additive_preprocessing_scale_policy(
        self,
        request: InterpretedDatasetBuildRequest,
    ) -> None:
        if request.corrected_preprocessing_output is not None:
            reject_external_corrected_output_after_downstream_preprocessing(
                request.preprocessing_plan.stage_order
            )
        self._additive_preprocessing_scale_guard.run(
            preprocessing_plan=request.preprocessing_plan,
            declared_input_scale_kind=request.declared_input_intensity_scale_kind,
            corrected_preprocessing_output=request.corrected_preprocessing_output,
        )
        _validate_quantitative_operation_contracts_before_preprocessing(
            request,
            preprocessor=self._preprocessor,
        )

    def _validate_analysis_ready_site_sequences(
        self,
        preprocessed: PreprocessedDatasetBuildTables,
    ) -> pd.DataFrame:
        return self._site_sequence_validator.run(
            site_metadata=preprocessed.site_metadata,
            preprocessing_trace=preprocessed.preprocessing_trace,
        )

    def _resolve_transformation_state(
        self,
        *,
        request: InterpretedDatasetBuildRequest,
        preprocessed: PreprocessedDatasetBuildTables,
        validated_site_metadata: pd.DataFrame,
    ) -> ResolvedDatasetTransformationState:
        return self._transformation_state_resolver.run(
            request=request,
            preprocessed=preprocessed,
            validated_site_metadata=validated_site_metadata,
        )

    def _run_protein_aware_preparation(
        self,
        *,
        request: InterpretedDatasetBuildRequest,
        transformed: ResolvedDatasetTransformationState,
        validated_site_metadata: pd.DataFrame,
    ) -> ProteinAwarePreparationResult | None:
        return self._protein_aware_preparation_runner.run(
            phospho=transformed.phospho,
            site_metadata=validated_site_metadata,
            total=transformed.total,
            intensity_scale_state=transformed.intensity_scale_state,
            plan=request.preprocessing_plan,
        )

    def _assemble_preprocessing_report(
        self,
        *,
        request: InterpretedDatasetBuildRequest,
        preprocessed: PreprocessedDatasetBuildTables,
        transformed: ResolvedDatasetTransformationState,
        protein_aware_preparation: ProteinAwarePreparationResult | None,
    ) -> DatasetPreprocessingReport:
        return self._preprocessing_report_assembler.run(
            row_counts=preprocessed.preprocessing_row_counts,
            operations=preprocessed.preprocessing_operations,
            row_audit=preprocessed.row_audit,
            duplicate_site_resolution=preprocessed.duplicate_site_resolution,
            metadata_conflicts=preprocessed.metadata_conflicts,
            comparison_group_stats=preprocessed.comparison_group_stats,
            comparison_pair_stats=preprocessed.comparison_pair_stats,
            preprocessing_trace=preprocessed.preprocessing_trace,
            site_sequence_derivation=request.site_sequence_derivation,
            input_site_count=int(request.site_metadata.shape[0]),
            final_dataset_rows=int(len(transformed.phospho.index)),
            intensity_scale_label=transformed.intensity_scale_state.label,
            intensity_scale_establishment=transformed.intensity_scale_establishment,
            quantitative_meaning_provenance=(
                transformed.quantitative_meaning_provenance
            ),
            declared_input_intensity_scale_kind=(
                None
                if request.declared_input_intensity_scale_kind is None
                else request.declared_input_intensity_scale_kind.value
            ),
            allow_suspicious_declared_input_intensity_scale=(
                request.allow_suspicious_declared_input_intensity_scale
            ),
            effective_declared_input_intensity_scale_diagnostic_policy=(
                _declared_scale_diagnostic_policy_label(request)
            ),
            quantitative_meaning=transformed.quantitative_meaning,
            peptide_evidence_resolution=request.peptide_evidence_resolution,
            preprocessing_plan=request.preprocessing_plan,
            sample_metadata=preprocessed.sample_metadata,
            batch_correction_metadata=preprocessed.batch_correction_metadata,
            batch_correction_report=preprocessed.batch_correction_report,
            protein_aware_preparation_report=(
                None
                if protein_aware_preparation is None
                else protein_aware_preparation.report
            ),
            matrix_shape_before=(
                int(request.phospho.shape[0]),
                int(request.phospho.shape[1]),
            ),
            matrix_shape_after=(
                int(transformed.phospho.shape[0]),
                int(transformed.phospho.shape[1]),
            ),
            intensity_transformation_event=(
                preprocessed.intensity_transformation_event
            ),
        )

    def _assemble_run_provenance(
        self,
        *,
        request: InterpretedDatasetBuildRequest,
        preprocessed: PreprocessedDatasetBuildTables,
        transformed: ResolvedDatasetTransformationState,
        validated_site_metadata: pd.DataFrame,
        protein_aware_preparation: ProteinAwarePreparationResult | None,
    ) -> RunProvenance:
        return self._provenance_assembler.run(
            request=request,
            preprocessed=preprocessed,
            validated_site_metadata=validated_site_metadata,
            resolved_phospho=transformed.phospho,
            resolved_total=transformed.total,
            preprocessing_trace=preprocessed.preprocessing_trace,
            intensity_scale_label=transformed.intensity_scale_state.label,
            intensity_scale_establishment=transformed.intensity_scale_establishment,
            quantitative_meaning_provenance=(
                transformed.quantitative_meaning_provenance
            ),
            quantitative_meaning=transformed.quantitative_meaning,
            allow_opaque_site_values=request.allow_opaque_site_values,
            protein_aware_preparation_report=(
                None
                if protein_aware_preparation is None
                else protein_aware_preparation.report
            ),
        )

    def _construct_dataset(
        self,
        *,
        request: InterpretedDatasetBuildRequest,
        preprocessed: PreprocessedDatasetBuildTables,
        transformed: ResolvedDatasetTransformationState,
        validated_site_metadata: pd.DataFrame,
        report: DatasetPreprocessingReport,
        protein_aware_preparation: ProteinAwarePreparationResult | None,
        provenance: RunProvenance,
    ) -> AnalysisReadyPhosphoDataset:
        return AnalysisReadyPhosphoDataset._from_builder_output(
            phospho=transformed.phospho,
            site_metadata=validated_site_metadata,
            sample_metadata=preprocessed.sample_metadata,
            total=transformed.total,
            comparisons=preprocessed.comparisons,
            imputation_observation_mask=preprocessed.imputation_observation_mask,
            organism=request.organism,
            intensity_scale_state=transformed.intensity_scale_state,
            processing_state=transformed.processing_state,
            preprocessing_report=report,
            protein_aware_preparation=protein_aware_preparation,
            provenance=provenance,
            allow_opaque_site_values=request.allow_opaque_site_values,
        )


def _declared_scale_diagnostic_policy_label(
    request: InterpretedDatasetBuildRequest,
) -> str:
    if request.allow_suspicious_declared_input_intensity_scale:
        return "warn"
    if request.declared_input_intensity_scale_kind is None:
        return "warn"
    if request.declared_input_intensity_scale_kind.value == "log2":
        return "error"
    return "warn"


def _validate_quantitative_operation_contracts_before_preprocessing(
    request: InterpretedDatasetBuildRequest,
    *,
    preprocessor: DatasetPreprocessorContract,
) -> QuantitativeContractState:
    preprocessor_validator = getattr(
        preprocessor,
        "validate_quantitative_contracts",
        None,
    )
    if callable(preprocessor_validator):
        result = preprocessor_validator(
            plan=request.preprocessing_plan,
            initial_quantitative_scale_kind=(
                request.declared_input_intensity_scale_kind
            ),
            initial_quantitative_meaning=request.quantitative_meaning,
        )
        if isinstance(result, QuantitativeContractState):
            return result
    quantitative_state = initial_quantitative_contract_state(
        declared_input_scale_kind=request.declared_input_intensity_scale_kind,
        explicit_quantitative_meaning=request.quantitative_meaning,
    )
    plan = request.preprocessing_plan
    for stage_key in plan.stage_order:
        metadata = get_preprocessing_stage_metadata(stage_key)
        interpreted = metadata.interpret(plan)
        quantitative_state = interpreted.quantitative_contract.validate_and_transition(
            quantitative_state,
            stage=stage_key,
            operation=interpreted.operation,
            evidence=None,
        )
    return quantitative_state


def _preprocessor_accepts_quantitative_contract_seed(
    preprocessor: DatasetPreprocessorContract,
) -> bool:
    try:
        parameters = inspect.signature(preprocessor.run).parameters
    except (TypeError, ValueError):
        return False
    if any(
        parameter.kind is inspect.Parameter.VAR_KEYWORD
        for parameter in parameters.values()
    ):
        return True
    return "initial_quantitative_scale_kind" in parameters

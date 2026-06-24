"""Internal executor for analysis-ready dataset construction."""

from __future__ import annotations

from phospy.science.datasets.builders.contracts import (
    DatasetPreprocessorContract,
    InterpretedDatasetBuildRequest,
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
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.transformations.contracts import Transformer
from phospy.science.transformations.transformers import IdentityTransformer


class DatasetBuildExecutor:
    """Construct `AnalysisReadyPhosphoDataset` from interpreted builder input."""

    def __init__(
        self,
        *,
        transformer: Transformer | None = None,
        intensity_scale_resolver: DatasetIntensityScaleResolver | None = None,
        preprocessor: DatasetPreprocessorContract | None = None,
        site_sequence_validator: AnalysisReadySiteSequenceValidator | None = None,
        transformation_state_resolver: DatasetTransformationStateResolver | None = None,
        preprocessing_report_assembler: DatasetPreprocessingReportAssembler
        | None = None,
        protein_aware_preparation_runner: DatasetProteinAwarePreparationRunner
        | None = None,
        provenance_assembler: DatasetRunProvenanceAssembler | None = None,
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

    def run(
        self, request: InterpretedDatasetBuildRequest
    ) -> AnalysisReadyPhosphoDataset:
        preprocessor_kwargs = {
            "phospho": request.phospho,
            "site_metadata": request.site_metadata,
            "sample_metadata": request.sample_metadata,
            "total": request.total,
            "plan": request.preprocessing_plan,
        }
        if request.corrected_preprocessing_output is not None:
            preprocessor_kwargs["corrected_preprocessing_output"] = (
                request.corrected_preprocessing_output
            )
        preprocessed = self._preprocessor.run(**preprocessor_kwargs)
        validated_site_metadata = self._site_sequence_validator.run(
            site_metadata=preprocessed.site_metadata,
            preprocessing_trace=preprocessed.preprocessing_trace,
        )
        transformed = self._transformation_state_resolver.run(
            request=request,
            preprocessed=preprocessed,
            validated_site_metadata=validated_site_metadata,
        )
        protein_aware_preparation = self._protein_aware_preparation_runner.run(
            phospho=transformed.phospho,
            site_metadata=validated_site_metadata,
            total=transformed.total,
            intensity_scale_state=transformed.intensity_scale_state,
            plan=request.preprocessing_plan,
        )
        report = self._preprocessing_report_assembler.run(
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
            declared_input_intensity_scale_kind=(
                None
                if request.declared_input_intensity_scale_kind is None
                else request.declared_input_intensity_scale_kind.value
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
        )
        provenance = self._provenance_assembler.run(
            request=request,
            preprocessed=preprocessed,
            validated_site_metadata=validated_site_metadata,
            resolved_phospho=transformed.phospho,
            resolved_total=transformed.total,
            preprocessing_trace=preprocessed.preprocessing_trace,
            intensity_scale_label=transformed.intensity_scale_state.label,
            intensity_scale_establishment=transformed.intensity_scale_establishment,
            quantitative_meaning=transformed.quantitative_meaning,
            allow_opaque_site_values=request.allow_opaque_site_values,
            protein_aware_preparation_report=(
                None
                if protein_aware_preparation is None
                else protein_aware_preparation.report
            ),
        )
        return AnalysisReadyPhosphoDataset._from_owned(
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

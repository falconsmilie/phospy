"""Internal executor for analysis-ready dataset construction."""

from __future__ import annotations

import pandas as pd

from phospy.science.datasets.builders.contracts import (
    DatasetPreprocessorContract,
    InterpretedDatasetBuildRequest,
)
from phospy.science.datasets.builders.preprocessing import DatasetPreprocessor
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
        self._provenance_assembler = (
            provenance_assembler or DatasetRunProvenanceAssembler()
        )

    def run(
        self, request: InterpretedDatasetBuildRequest
    ) -> AnalysisReadyPhosphoDataset:
        preprocessed = self._preprocessor.run(
            phospho=request.phospho,
            site_metadata=request.site_metadata,
            sample_metadata=request.sample_metadata,
            total=request.total,
            plan=request.preprocessing_plan,
        )
        validated_site_metadata = self._site_sequence_validator.run(
            site_metadata=preprocessed.site_metadata,
            preprocessing_trace=preprocessed.preprocessing_trace,
        )
        transformed = self._transformation_state_resolver.run(
            request=request,
            preprocessed=preprocessed,
            validated_site_metadata=validated_site_metadata,
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
        )
        public_site_metadata = _project_gene_site_fallback_site_key_to_display_index(
            phospho=transformed.phospho,
            site_metadata=validated_site_metadata,
        )
        return AnalysisReadyPhosphoDataset._from_owned(
            phospho=transformed.phospho,
            site_metadata=public_site_metadata,
            sample_metadata=preprocessed.sample_metadata,
            total=transformed.total,
            comparisons=preprocessed.comparisons,
            organism=request.organism,
            intensity_scale_state=transformed.intensity_scale_state,
            processing_state=transformed.processing_state,
            preprocessing_report=report,
            provenance=provenance,
            allow_opaque_site_values=request.allow_opaque_site_values,
        )


def _project_gene_site_fallback_site_key_to_display_index(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
) -> pd.DataFrame:
    if "site_key" not in site_metadata.columns:
        return site_metadata
    has_explicit_protein_identity = any(
        column in site_metadata.columns
        for column in ("protein_id", "protein_accession", "isoform_id")
    )
    if has_explicit_protein_identity:
        return site_metadata
    site_metadata.loc[:, "site_key"] = phospho.index.astype(str).tolist()
    return site_metadata

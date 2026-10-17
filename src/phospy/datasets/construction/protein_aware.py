"""Compose passive protein-aware output into durable Dataset records."""

from __future__ import annotations

from phospy.datasets.protein_aware import (
    ProteinAwarePreparationReport,
    ProteinAwarePreparationResult,
)
from phospy.errors.validation import DatasetValidationError
from phospy.preprocessing.dataset_construction import (
    ProteinAwarePreparationStageOutput,
)
from phospy.run_lineage.models import RunProvenance


def compose_protein_aware_preparation_result(
    stage_output: ProteinAwarePreparationStageOutput,
    *,
    provenance: RunProvenance | None,
) -> ProteinAwarePreparationResult:
    """Compose a durable Dataset result from a passive preprocessing output."""

    if not isinstance(stage_output, ProteinAwarePreparationStageOutput):
        raise DatasetValidationError(
            "protein-aware preparation stage output must be "
            "ProteinAwarePreparationStageOutput"
        )
    report = ProteinAwarePreparationReport(
        site_eligibility=stage_output.site_eligibility,
        mapping_diagnostics=stage_output.mapping_diagnostics,
        sample_alignment=stage_output.sample_alignment,
        transformation_state=stage_output.transformation_state,
        preparation_policy=stage_output.preparation_policy,
        protein_mapping_policy=stage_output.protein_mapping_policy,
        policy_parameters=stage_output.policy_parameters,
        provenance=provenance,
    )
    matched_pairs = stage_output.matched_pairs
    protein_covariate_matrix = stage_output.protein_covariate_matrix
    return ProteinAwarePreparationResult._from_owned(
        matched_pairs=matched_pairs,
        protein_covariate_matrix=protein_covariate_matrix,
        report=report,
    )


__all__ = ["compose_protein_aware_preparation_result"]

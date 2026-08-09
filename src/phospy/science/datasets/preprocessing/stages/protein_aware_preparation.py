"""Executable protein-aware preparation preprocessing stage."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import DatasetValidationError
from phospy.science.configs.preprocessing.total_protein import (
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_ALLOW_MISSING_WITH_REPORT,
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS,
    DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED,
    DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS,
    DatasetProteinAwarePreparationConfig,
)
from phospy.science.datasets.preprocessing.protein_aware_alignment import (
    PROTEIN_AWARE_REASON_AMBIGUOUS_PROTEIN_MAPPING,
    PROTEIN_AWARE_REASON_MATCHED_PROTEIN_AVAILABLE,
    PROTEIN_AWARE_REASON_MISSING_SITE_PROTEIN_IDENTIFIER,
    PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW,
    ProteinAwareAlignmentConfig,
    ProteinAwareAlignmentEligibilityResolver,
    ProteinAwarePreparationEligibility,
    ProteinAwareSampleAlignmentDiagnostics,
    ProteinAwareSiteEligibilityDiagnostic,
    ProteinAwareTransformationStateDiagnostics,
)
from phospy.science.datasets.preprocessing.protein_aware_models import (
    PROTEIN_AWARE_AMBIGUOUS_MAPPING_DIAGNOSTIC_COLUMNS,
    PROTEIN_AWARE_MATCHED_PAIR_COLUMNS,
    PROTEIN_AWARE_MISSING_PROTEIN_ABUNDANCE_DIAGNOSTIC_COLUMNS,
    ProteinAwareMappingDiagnostics,
    ProteinAwarePreparationReport,
    ProteinAwarePreparationResult,
    ProteinAwareSiteEligibility,
)
from phospy.science.datasets.preprocessing.protein_mapping import (
    ProteinMappingConfig,
    ProteinMappingRecord,
    ProteinMappingResolver,
    ProteinMappingResult,
    ProteinMappingStatus,
)
from phospy.science.transformations.models import IntensityScaleState


class ProteinAwarePreparationStage:
    """Prepare audited protein-aware modelling inputs below the dataset boundary.

    This collaborator resolves phosphosite-to-protein mappings, diagnoses sample
    and transformation compatibility, and builds a total-protein covariate matrix
    for eligible sites. It does not alter phosphosite intensities, subtract total
    protein, fit models, or run differential analysis.
    """

    def __init__(
        self,
        *,
        mapping_resolver: ProteinMappingResolver | None = None,
        eligibility_resolver: ProteinAwareAlignmentEligibilityResolver | None = None,
    ) -> None:
        self._mapping_resolver = mapping_resolver or ProteinMappingResolver()
        self._eligibility_resolver = (
            eligibility_resolver or ProteinAwareAlignmentEligibilityResolver()
        )

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        total: pd.DataFrame | None,
        transformation_state: IntensityScaleState,
        config: DatasetProteinAwarePreparationConfig,
        mapping_config: ProteinMappingConfig | None = None,
    ) -> ProteinAwarePreparationResult | None:
        """Return prepared protein-aware inputs, or ``None`` when disabled."""

        if not isinstance(config, DatasetProteinAwarePreparationConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config."
                "protein_aware_preparation must be "
                "DatasetProteinAwarePreparationConfig"
            )
        if config.policy == DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED:
            return None
        if (
            config.policy
            != DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config."
                "protein_aware_preparation.policy contains an unsupported value"
            )
        if total is None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config."
                "protein_aware_preparation.policy='prepare_model_inputs' requires "
                "total input data"
            )

        mapping_result = self._mapping_resolver.run(
            site_metadata=site_metadata,
            phospho_matrix_index=phospho.index,
            total_protein_matrix_index=total.index,
            config=mapping_config,
        )
        eligibility = self._eligibility_resolver.run(
            phospho=phospho,
            total=total,
            mapping_result=mapping_result,
            intensity_scale_state=transformation_state,
            config=ProteinAwareAlignmentConfig(
                protein_mapping_policy=config.protein_mapping_policy,
                allow_reordered_samples=False,
            ),
        )
        report = _build_preparation_report(
            mapping_result=mapping_result,
            site_eligibility=eligibility.site_eligibility,
            sample_alignment=eligibility.sample_alignment,
            transformation_state=eligibility.transformation_state,
            preparation_policy=config.policy,
            protein_mapping_policy=config.protein_mapping_policy,
        )
        matched_pairs = _matched_pairs_from_eligibility(eligibility.site_eligibility)
        protein_covariates = _protein_covariate_matrix(
            total=total,
            phospho_columns=phospho.columns,
            matched_pairs=matched_pairs,
            sample_order_compatible=(
                eligibility.sample_alignment.sample_order_compatible
            ),
        )
        return ProteinAwarePreparationResult._from_owned(
            matched_pairs=matched_pairs,
            protein_covariate_matrix=protein_covariates,
            report=report,
        )


def _build_preparation_report(
    *,
    mapping_result: ProteinMappingResult,
    site_eligibility: tuple[ProteinAwareSiteEligibilityDiagnostic, ...],
    sample_alignment: ProteinAwareSampleAlignmentDiagnostics,
    transformation_state: ProteinAwareTransformationStateDiagnostics,
    preparation_policy: str,
    protein_mapping_policy: str,
) -> ProteinAwarePreparationReport:
    return ProteinAwarePreparationReport(
        site_eligibility=tuple(
            ProteinAwareSiteEligibility(
                site_key=diagnostic.site_key,
                eligibility=diagnostic.eligibility,
                mapping_status=diagnostic.mapping_status,
                protein_identifier=diagnostic.protein_identifier,
                total_protein_row_key=diagnostic.total_protein_row_key,
                reasons=diagnostic.reasons,
            )
            for diagnostic in site_eligibility
        ),
        mapping_diagnostics=_mapping_diagnostics_from_records(mapping_result.records),
        sample_alignment=sample_alignment,
        transformation_state=transformation_state,
        preparation_policy=preparation_policy,
        protein_mapping_policy=protein_mapping_policy,
        policy_parameters={
            "preparation_mode": "aligned_model_input_preparation_only",
            "allow_reordered_samples": False,
            "missing_mapping_policy": (
                DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_ALLOW_MISSING_WITH_REPORT
                if protein_mapping_policy
                == DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_ALLOW_MISSING_WITH_REPORT
                else DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS
            ),
            "modifies_phospho_matrix": False,
            "performs_total_protein_subtraction": False,
            "performs_normalisation": False,
            "performs_model_adjustment": False,
            "performs_differential_modelling": False,
            "claims_msstatsptm_equivalence": False,
            "limitations": (
                "preparation-only; aligned phosphosite/protein inputs and diagnostics",
                "does not subtract total protein from phosphosite intensities",
                "does not normalise phosphosite intensities",
                "does not run joint PTM/protein differential modelling",
                "does not claim MSstatsPTM-style inference or equivalence",
            ),
        },
    )


def _mapping_diagnostics_from_records(
    records: tuple[ProteinMappingRecord, ...],
) -> ProteinAwareMappingDiagnostics:
    missing_rows: list[dict[str, object]] = []
    ambiguous_rows: list[dict[str, object]] = []
    for record in records:
        reason = _mapping_reason(record.status)
        if record.status in {
            ProteinMappingStatus.MISSING_SITE_PROTEIN_IDENTIFIER,
            ProteinMappingStatus.MISSING_TOTAL_PROTEIN_ROW,
        }:
            missing_rows.append(
                {
                    "site_key": record.site_key,
                    "protein_identifier": record.protein_identifier,
                    "mapping_status": record.status.value,
                    "reason": reason,
                }
            )
            continue
        if record.status in {
            ProteinMappingStatus.AMBIGUOUS_SITE_PROTEIN_MAPPING,
            ProteinMappingStatus.AMBIGUOUS_TOTAL_PROTEIN_MAPPING,
        }:
            ambiguous_rows.append(
                {
                    "site_key": record.site_key,
                    "mapping_status": record.status.value,
                    "protein_identifier": record.protein_identifier,
                    "candidate_protein_identifiers": (
                        record.candidate_protein_identifiers
                    ),
                    "candidate_total_protein_row_keys": (
                        record.candidate_total_protein_row_keys
                    ),
                    "reason": reason,
                }
            )
    return ProteinAwareMappingDiagnostics._from_owned(
        missing_protein_abundance=pd.DataFrame(
            missing_rows,
            columns=list(PROTEIN_AWARE_MISSING_PROTEIN_ABUNDANCE_DIAGNOSTIC_COLUMNS),
        ),
        ambiguous_mapping=pd.DataFrame(
            ambiguous_rows,
            columns=list(PROTEIN_AWARE_AMBIGUOUS_MAPPING_DIAGNOSTIC_COLUMNS),
        ),
    )


def _matched_pairs_from_eligibility(
    site_eligibility: tuple[ProteinAwareSiteEligibilityDiagnostic, ...],
) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for diagnostic in site_eligibility:
        if (
            diagnostic.eligibility
            is not ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION
        ):
            continue
        if diagnostic.protein_identifier is None:
            raise DatasetValidationError(
                "eligible protein-aware site is missing protein_identifier"
            )
        if diagnostic.total_protein_row_key is None:
            raise DatasetValidationError(
                "eligible protein-aware site is missing total_protein_row_key"
            )
        rows.append(
            {
                "site_key": diagnostic.site_key,
                "protein_identifier": diagnostic.protein_identifier,
                "total_protein_row_key": diagnostic.total_protein_row_key,
            }
        )
    return pd.DataFrame(rows, columns=list(PROTEIN_AWARE_MATCHED_PAIR_COLUMNS))


def _protein_covariate_matrix(
    *,
    total: pd.DataFrame,
    phospho_columns: pd.Index,
    matched_pairs: pd.DataFrame,
    sample_order_compatible: bool,
) -> pd.DataFrame:
    output_columns = pd.Index(phospho_columns.tolist())
    if matched_pairs.empty or not sample_order_compatible:
        return pd.DataFrame(
            index=pd.Index([], name=total.index.name),
            columns=output_columns,
            dtype=float,
        )
    row_keys = _dedupe_preserving_order(
        matched_pairs.loc[:, "total_protein_row_key"].astype(str).tolist()
    )
    return total.loc[list(row_keys), output_columns.tolist()].copy(deep=True)


def _mapping_reason(status: ProteinMappingStatus) -> str:
    if status is ProteinMappingStatus.MATCHED:
        return PROTEIN_AWARE_REASON_MATCHED_PROTEIN_AVAILABLE
    if status is ProteinMappingStatus.MISSING_SITE_PROTEIN_IDENTIFIER:
        return PROTEIN_AWARE_REASON_MISSING_SITE_PROTEIN_IDENTIFIER
    if status is ProteinMappingStatus.MISSING_TOTAL_PROTEIN_ROW:
        return PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW
    if status in {
        ProteinMappingStatus.AMBIGUOUS_SITE_PROTEIN_MAPPING,
        ProteinMappingStatus.AMBIGUOUS_TOTAL_PROTEIN_MAPPING,
    }:
        return PROTEIN_AWARE_REASON_AMBIGUOUS_PROTEIN_MAPPING
    return status.value


def _dedupe_preserving_order(values: Sequence[str]) -> tuple[str, ...]:
    seen: set[str] = set()
    deduped: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        deduped.append(value)
    return tuple(deduped)


__all__ = ["ProteinAwarePreparationStage"]

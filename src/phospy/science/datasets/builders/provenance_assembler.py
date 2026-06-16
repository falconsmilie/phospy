"""Run-provenance assembly for dataset builder execution."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_optional_table
from phospy.provenance.models import (
    PREPROCESSING_STAGE_DETERMINISM_PURE,
    JsonValue,
    PreprocessingStageProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.provenance.scientific_policy_models import ScientificPolicyRecord
from phospy.science.datasets.builders.contracts import (
    InterpretedDatasetBuildRequest,
    PreprocessedDatasetBuildTables,
)
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_LOCALISATION,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageExecution,
    PreprocessingStageOrderResolution,
    TotalProteinCorrectionIdentityPolicy,
)
from phospy.science.datasets.preprocessing.protein_aware_preparation import (
    ProteinAwarePreparationReport,
)
from phospy.science.datasets.preprocessing.scientific_policies import (
    PreprocessingStageOrderPolicy,
    build_duplicate_site_resolution_policy,
)

_SUPPORTED_PREPROCESSING_STAGE_ORDER = (
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_LOCALISATION,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
)


class DatasetRunProvenanceAssembler:
    """Assemble `RunProvenance` for the dataset builder workflow."""

    def run(
        self,
        *,
        request: InterpretedDatasetBuildRequest,
        preprocessed: PreprocessedDatasetBuildTables,
        validated_site_metadata: pd.DataFrame,
        resolved_phospho: pd.DataFrame,
        resolved_total: pd.DataFrame | None,
        preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
        intensity_scale_label: str,
        intensity_scale_establishment: Mapping[str, object],
        quantitative_meaning: str,
        allow_opaque_site_values: bool,
        protein_aware_preparation_report: ProteinAwarePreparationReport | None = None,
    ) -> RunProvenance:
        input_tables = _collect_fingerprints(
            (
                ("dataset.phospho", request.phospho),
                ("dataset.site_metadata", request.site_metadata),
                ("dataset.sample_metadata", request.sample_metadata),
                ("dataset.total", request.total),
            )
        )
        output_tables = _collect_fingerprints(
            (
                ("dataset.phospho", resolved_phospho),
                ("dataset.site_metadata", validated_site_metadata),
                ("dataset.sample_metadata", preprocessed.sample_metadata),
                ("dataset.total", resolved_total),
                ("dataset.comparisons", preprocessed.comparisons),
                (
                    "dataset.imputation_observation_mask",
                    preprocessed.imputation_observation_mask,
                ),
            )
        )
        workflow_parameters: dict[str, object] = {
            "preprocessing_plan": _preprocessing_plan_to_payload(
                request.preprocessing_plan
            ),
            "intensity_scale_label": intensity_scale_label,
            "intensity_scale_establishment": dict(intensity_scale_establishment),
            "quantitative_meaning": quantitative_meaning,
            "site_identifier_normalisation": (
                None
                if request.site_identifier_normalisation is None
                else request.site_identifier_normalisation.to_payload()
            ),
            "site_sequence_derivation": request.site_sequence_derivation,
            "site_resolution_mode": request.site_resolution_mode,
            "multi_site_policy": request.multi_site_policy,
            "peptide_evidence_resolution": request.peptide_evidence_resolution,
            "site_token_validation": (
                {"mode": "opaque_opt_in"}
                if allow_opaque_site_values
                else {"mode": "strict_sty_residue_position"}
            ),
        }
        if protein_aware_preparation_report is not None:
            workflow_parameters["protein_aware_preparation"] = (
                _protein_aware_preparation_to_payload(protein_aware_preparation_report)
            )
        return RunProvenance(
            environment=collect_environment_provenance(),
            input_tables=input_tables,
            preprocessing_stages=_stage_trace_to_provenance(
                preprocessing_trace,
                quantitative_meaning=quantitative_meaning,
            ),
            reference=None,
            workflow_name="dataset_builder",
            workflow_parameters=workflow_parameters,
            random_state=None,
            random_seed_policy=None,
            output_tables=output_tables,
            scientific_policies=_dataset_scientific_policies(
                request.preprocessing_plan
            ),
        )


def _collect_fingerprints(
    entries: tuple[tuple[str, pd.DataFrame | None], ...],
) -> tuple[TableFingerprint, ...]:
    fingerprints: list[TableFingerprint] = []
    for name, table in entries:
        fingerprint = fingerprint_optional_table(table, name=name)
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _stage_trace_to_provenance(
    trace: tuple[PreprocessingStageExecution, ...] | None,
    *,
    quantitative_meaning: str,
) -> tuple[PreprocessingStageProvenance, ...]:
    if trace is None:
        return ()
    return tuple(
        PreprocessingStageProvenance(
            stage=item.stage,
            operation=item.operation,
            parameters=dict(item.parameters),
            input_shape=item.input_shape,
            output_shape=item.output_shape,
            input_hash=item.input_hash,
            output_hash=item.output_hash,
            phospho_input_hash=item.phospho_input_hash,
            phospho_output_hash=item.phospho_output_hash,
            dropped_row_ids=item.dropped_row_ids,
            dropped_row_count=int(item.dropped_row_count),
            schema_version=int(item.schema_version),
            consumed_input_tables=tuple(item.consumed_input_tables),
            produced_output_tables=tuple(item.produced_output_tables),
            backend=item.backend,
            random_seed=item.random_seed,
            determinism=(
                str(item.determinism).strip()
                if str(item.determinism).strip()
                else PREPROCESSING_STAGE_DETERMINISM_PURE
            ),
            imputed_cell_count=int(item.imputed_cell_count),
            imputed_row_ids=item.imputed_row_ids,
            notes=item.notes,
            diagnostics=_to_stage_diagnostics_payload(
                item.diagnostics,
                quantitative_meaning=quantitative_meaning,
            ),
        )
        for item in trace
    )


def _to_stage_diagnostics_payload(
    values: Mapping[str, object],
    *,
    quantitative_meaning: str,
) -> dict[str, JsonValue]:
    payload = _to_json_mapping(values)
    transformer_state = payload.get("transformer_state")
    if isinstance(transformer_state, dict):
        state_payload = dict(transformer_state)
        state_payload["quantity"] = quantitative_meaning
        payload["transformer_state"] = state_payload
    return payload


def _to_json_mapping(values: Mapping[str, object]) -> dict[str, JsonValue]:
    return {str(key): _to_json_value(value) for key, value in values.items()}


def _to_json_value(value: object) -> JsonValue:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        return {str(key): _to_json_value(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_to_json_value(item) for item in value]
    if isinstance(value, list):
        return [_to_json_value(item) for item in value]
    return str(value)


def _preprocessing_plan_to_payload(plan: PreprocessingPlan) -> dict[str, object]:
    payload: dict[str, object] = {
        "intensity_transform_policy": plan.intensity_transform_policy.value,
        "intensity_transform_pseudocount": float(plan.intensity_transform_pseudocount),
        "normalisation_policy": plan.normalisation_policy.value,
        "missing_data_policy": plan.missing_data_policy.value,
        "missing_data_min_observed_values": plan.missing_data_min_observed_values,
        "missing_data_q": plan.missing_data_q,
        "missing_data_width": plan.missing_data_width,
        "missing_data_seed": plan.missing_data_seed,
        "missing_data_k": plan.missing_data_k,
        "missing_data_distance": plan.missing_data_distance,
        "missing_data_max_missing_fraction_per_row": (
            plan.missing_data_max_missing_fraction_per_row
        ),
        "localisation_mode": plan.localisation_mode.value,
        "localisation_min_confidence": float(plan.localisation_min_confidence),
        "localisation_confidence_column": plan.localisation_confidence_column,
        "localisation_waiver_reason": plan.localisation_waiver_reason,
        "site_sequence_resolution_enabled": plan.site_sequence_resolution_enabled,
        "site_sequence_resolution_fasta_path": plan.site_sequence_resolution_fasta_path,
        "site_sequence_resolution_mode": plan.site_sequence_resolution_mode.value,
        "site_sequence_resolution_flank_size": int(
            plan.site_sequence_resolution_flank_size
        ),
        "site_sequence_resolution_accession_column": (
            plan.site_sequence_resolution_accession_column
        ),
        "site_sequence_resolution_site_column": (
            plan.site_sequence_resolution_site_column
        ),
        "total_protein_correction_policy": plan.total_protein_correction_policy.value,
        "total_protein_correction_identity_policy": (
            _total_correction_identity_policy_to_payload(
                plan.total_protein_correction_identity_policy
            )
        ),
        "protein_aware_preparation_policy": (plan.protein_aware_preparation_policy),
        "protein_aware_preparation_mapping_policy": (
            plan.protein_aware_preparation_mapping_policy
        ),
        "site_matrix_policy": plan.site_matrix_policy.value,
        "comparison_building_policy": plan.comparison_building_policy.value,
        "site_matrix_duplicate_site_policy": plan.site_matrix_duplicate_site_policy.value,
        "site_matrix_missing_data_policy": plan.site_matrix_missing_data_policy.value,
        "site_matrix_minimum_observed_values": plan.site_matrix_minimum_observed_values,
        "comparison_sample_group_column": plan.comparison_sample_group_column,
        "comparison_pairs": (
            None if plan.comparison_pairs is None else list(plan.comparison_pairs)
        ),
        "ruv_readiness_enabled": bool(plan.ruv_readiness_enabled),
        "ruv_readiness_control_feature_column": (
            plan.ruv_readiness_control_feature_column
        ),
        "ruv_readiness_replicate_group_column": (
            plan.ruv_readiness_replicate_group_column
        ),
        "ruv_readiness_batch_column": plan.ruv_readiness_batch_column,
        "batch_correction_method": plan.batch_correction_method,
        "batch_correction_batch_column": plan.batch_correction_batch_column,
        "batch_correction_condition_column": plan.batch_correction_condition_column,
        "batch_correction_preserve_condition_effects": (
            plan.batch_correction_preserve_condition_effects
        ),
        "stage_order": list(plan.stage_order),
        "resolved_stage_order": _stage_order_resolution_to_payload(
            plan.stage_order_resolution
        ),
    }
    return payload


def _stage_order_resolution_to_payload(
    stage_order_resolution: tuple[PreprocessingStageOrderResolution, ...],
) -> list[dict[str, object]]:
    return [
        {
            "stage": item.stage,
            "order_index": int(item.order_index),
            "rationale": str(item.rationale),
        }
        for item in stage_order_resolution
    ]


def _total_correction_identity_policy_to_payload(
    policy: TotalProteinCorrectionIdentityPolicy,
) -> dict[str, object]:
    return {
        "mode": str(policy.mode),
        "matching_policy": str(policy.matching_policy),
        "phosphosite_key": policy.phosphosite_key,
        "total_protein_key": policy.total_protein_key,
        "mapping_phosphosite_key": policy.mapping_phosphosite_key,
        "mapping_total_protein_key": policy.mapping_total_protein_key,
        "mapping_table_fingerprint": policy.mapping_table_fingerprint,
        "mapping_table_row_count": (
            None if policy.mapping_table is None else int(len(policy.mapping_table))
        ),
        "duplicate_policy": str(policy.duplicate_policy),
        "unmatched_policy": str(policy.unmatched_policy),
    }


def _protein_aware_preparation_to_payload(
    report: ProteinAwarePreparationReport,
) -> dict[str, object]:
    transformation_state = report.transformation_state
    return {
        "status": "prepared",
        "preparation_policy": report.preparation_policy,
        "protein_mapping_policy": report.protein_mapping_policy,
        "eligible_site_count": int(len(report.eligible_site_keys)),
        "fallback_site_count": int(len(report.fallback_site_keys)),
        "excluded_site_count": int(len(report.excluded_site_keys)),
        "sample_order_compatible": report.sample_alignment.sample_order_compatible,
        "exact_sample_order_match": report.sample_alignment.exact_sample_order_match,
        "transformation_state_compatible": (
            None if transformation_state is None else transformation_state.compatible
        ),
        "modifies_phospho_matrix": False,
        "performs_total_protein_subtraction": False,
        "performs_normalisation": False,
        "performs_differential_modelling": False,
        "claims_msstatsptm_equivalence": False,
        "limitations": [
            "preparation-only; aligned phosphosite/protein inputs and diagnostics",
            "does not subtract total protein from phosphosite intensities",
            "does not normalise phosphosite intensities",
            "does not run joint PTM/protein differential modelling",
            "does not claim MSstatsPTM-style inference or equivalence",
        ],
    }


def _dataset_scientific_policies(
    preprocessing_plan: PreprocessingPlan,
) -> tuple[ScientificPolicyRecord, ...]:
    policies = [
        PreprocessingStageOrderPolicy(
            configured_stage_order=tuple(
                str(stage) for stage in preprocessing_plan.stage_order
            ),
            default_stage_order=tuple(
                str(stage) for stage in DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT
            ),
            supported_stage_order=_SUPPORTED_PREPROCESSING_STAGE_ORDER,
        ).record,
    ]
    if DATASET_PREPROCESSING_STAGE_SITE_MATRIX in preprocessing_plan.stage_order:
        policies.append(
            build_duplicate_site_resolution_policy(
                duplicate_site_policy=(
                    preprocessing_plan.site_matrix_duplicate_site_policy.value
                )
            )
        )
    return tuple(policies)

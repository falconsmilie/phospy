"""Resolve validated SPS/RUV-style batch-correction intent into an execution plan."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import cast

import pandas as pd

from phospy.contracts.configs.preprocessing import (
    SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER,
    SUPPORTED_INTERNAL_BATCH_CORRECTION_STAGE_ORDER,
    CorrectionMissingnessPolicy,
    InternalBatchCorrectionRequest,
    InternalBatchCorrectionStageOrder,
    ObservationMask,
    TemporaryImputationPolicy,
)
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteEligibility,
    ControlSiteMapping,
)
from phospy.validation.datasets.batch_correction import (
    ReplicateStructureDiagnostics,
    ResolvedBatchDesignMetadata,
)

REPLICATE_METADATA_ROLE = "provenance_only"
REPLICATE_METADATA_ROLE_DESCRIPTION = (
    "replicate metadata is validated and recorded for provenance and diagnostics "
    "only; it is not used for numerical unwanted-factor estimation and does not "
    "enable RUV-III or replicate-aware RUV-III semantics"
)


@dataclass(frozen=True, slots=True)
class EligibleControlSiteRow:
    """One eligible control row resolved to the dataset site-key axis."""

    site_key: str
    row_position: int
    label: str | None = None
    weight: float | None = None
    group: str | None = None

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible payload for provenance and diagnostics."""

        return {
            "site_key": self.site_key,
            "row_position": self.row_position,
            "label": self.label,
            "weight": self.weight,
            "group": self.group,
        }


@dataclass(frozen=True, slots=True)
class ReplicateStructure:
    """Replicate labels aligned to the batch-correction sample order."""

    replicate_column: str | None
    replicate_by_sample: Mapping[str, str] | None
    replicate_labels: tuple[str, ...] | None
    replicate_groups: Mapping[str, tuple[str, ...]]
    structure_diagnostics: ReplicateStructureDiagnostics | None = None

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible payload."""

        return {
            "replicate_column": self.replicate_column,
            "replicate_metadata_role": REPLICATE_METADATA_ROLE,
            "replicate_metadata_role_description": (
                REPLICATE_METADATA_ROLE_DESCRIPTION
            ),
            "used_for_numerical_factor_estimation": False,
            "ruv_iii_semantics_enabled": False,
            "replicate_by_sample": (
                None
                if self.replicate_by_sample is None
                else dict(self.replicate_by_sample)
            ),
            "replicate_labels": (
                None if self.replicate_labels is None else list(self.replicate_labels)
            ),
            "replicate_groups": {
                replicate: list(samples)
                for replicate, samples in self.replicate_groups.items()
            },
            "structure_diagnostics": (
                None
                if self.structure_diagnostics is None
                else self.structure_diagnostics.to_payload()
            ),
        }


@dataclass(frozen=True, slots=True)
class BatchCorrectionDiagnosticRequirements:
    """Diagnostic payloads a future executor must preserve or emit."""

    diagnostics_enabled: bool
    required_payloads: tuple[str, ...]

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible payload."""

        return {
            "diagnostics_enabled": self.diagnostics_enabled,
            "required_payloads": list(self.required_payloads),
        }


@dataclass(frozen=True, slots=True)
class ResolvedBatchCorrectionPlan:
    """Execution plan for SPS/RUV-style correction without numerical correction."""

    method: str
    resolved_design_matrix: pd.DataFrame
    batch_terms: tuple[str, ...]
    condition_terms_to_preserve: tuple[str, ...]
    replicate_structure: ReplicateStructure
    eligible_control_site_rows: tuple[EligibleControlSiteRow, ...]
    observation_mask: ObservationMask
    temporary_imputation_policy: TemporaryImputationPolicy
    n_unwanted_factors: int | None
    stage_order: tuple[str, ...]
    stage_order_policy: str
    diagnostic_requirements: BatchCorrectionDiagnosticRequirements
    provenance_seed_data: Mapping[str, object]

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible plan payload."""

        return {
            "method": self.method,
            "resolved_design_matrix": _frame_payload(self.resolved_design_matrix),
            "batch_terms": list(self.batch_terms),
            "condition_terms_to_preserve": list(self.condition_terms_to_preserve),
            "replicate_structure": self.replicate_structure.to_payload(),
            "eligible_control_site_rows": [
                row.to_payload() for row in self.eligible_control_site_rows
            ],
            "observation_mask": self.observation_mask.to_payload(),
            "temporary_imputation_policy": _temporary_imputation_payload(
                self.temporary_imputation_policy
            ),
            "n_unwanted_factors": self.n_unwanted_factors,
            "executed_stage_order": list(self.stage_order),
            "stage_order": list(self.stage_order),
            "requested_stage_order": self.stage_order_policy,
            "stage_order_policy": self.stage_order_policy,
            "diagnostic_requirements": self.diagnostic_requirements.to_payload(),
            "provenance_seed_data": dict(self.provenance_seed_data),
        }


class BatchCorrectionPlanInterpreter:
    """Resolve validated correction intent into a deterministic correction plan."""

    def run(
        self,
        *,
        config: InternalBatchCorrectionRequest,
        dataset_metadata: ResolvedBatchDesignMetadata,
        control_site_mapping: ControlSiteMapping,
        missingness_policy: CorrectionMissingnessPolicy,
    ) -> ResolvedBatchCorrectionPlan:
        sample_order = tuple(dataset_metadata.sample_order)
        batch_labels = tuple(dataset_metadata.batch_by_sample[s] for s in sample_order)
        condition_labels = tuple(
            dataset_metadata.condition_by_sample[s] for s in sample_order
        )
        condition_design = _treatment_coded_frame(
            labels=condition_labels,
            sample_order=sample_order,
            term_prefix="condition",
            include_intercept=True,
        )
        batch_design = _treatment_coded_frame(
            labels=batch_labels,
            sample_order=sample_order,
            term_prefix="batch",
            include_intercept=False,
        )
        resolved_design = pd.concat((condition_design, batch_design), axis=1)
        eligible_controls = _eligible_control_site_rows(control_site_mapping)
        observation_mask = _resolve_observation_mask(
            missingness_policy=missingness_policy,
            control_site_mapping=control_site_mapping,
            sample_order=sample_order,
        )
        replicate_structure = _replicate_structure(
            replicate_column=config.replicate_column,
            replicate_by_sample=dataset_metadata.replicate_by_sample,
            structure_diagnostics=(dataset_metadata.replicate_structure_diagnostics),
            sample_order=sample_order,
        )
        diagnostic_requirements = _diagnostic_requirements(
            diagnostics_enabled=bool(config.diagnostics_enabled)
        )
        stage_order = _stage_order(config.stage_order)
        plan = ResolvedBatchCorrectionPlan(
            method=config.method.value,
            resolved_design_matrix=cast(pd.DataFrame, resolved_design.copy(deep=True)),
            batch_terms=tuple(str(column) for column in batch_design.columns),
            condition_terms_to_preserve=tuple(
                str(column) for column in condition_design.columns
            ),
            replicate_structure=replicate_structure,
            eligible_control_site_rows=eligible_controls,
            observation_mask=observation_mask,
            temporary_imputation_policy=missingness_policy.temporary_imputation,
            n_unwanted_factors=config.n_unwanted_factors,
            stage_order=stage_order,
            stage_order_policy=config.stage_order.value,
            diagnostic_requirements=diagnostic_requirements,
            provenance_seed_data={},
        )
        provenance_seed_data = _provenance_seed_data(
            config=config,
            dataset_metadata=dataset_metadata,
            plan=plan,
            control_site_mapping=control_site_mapping,
            missingness_policy=missingness_policy,
        )
        return ResolvedBatchCorrectionPlan(
            method=plan.method,
            resolved_design_matrix=plan.resolved_design_matrix,
            batch_terms=plan.batch_terms,
            condition_terms_to_preserve=plan.condition_terms_to_preserve,
            replicate_structure=plan.replicate_structure,
            eligible_control_site_rows=plan.eligible_control_site_rows,
            observation_mask=plan.observation_mask,
            temporary_imputation_policy=plan.temporary_imputation_policy,
            n_unwanted_factors=plan.n_unwanted_factors,
            stage_order=plan.stage_order,
            stage_order_policy=plan.stage_order_policy,
            diagnostic_requirements=plan.diagnostic_requirements,
            provenance_seed_data=provenance_seed_data,
        )


def _treatment_coded_frame(
    *,
    labels: Sequence[str],
    sample_order: tuple[str, ...],
    term_prefix: str,
    include_intercept: bool,
) -> pd.DataFrame:
    levels = _levels_in_order(labels)
    columns: list[str] = []
    if include_intercept:
        columns.append("intercept")
    columns.extend(f"{term_prefix}[{level}]" for level in levels[1:])

    rows: list[list[float]] = []
    for label in labels:
        row: list[float] = []
        if include_intercept:
            row.append(1.0)
        row.extend(1.0 if label == level else 0.0 for level in levels[1:])
        rows.append(row)
    return pd.DataFrame(
        rows,
        index=pd.Index(sample_order, name="sample"),
        columns=pd.Index(columns, name="term"),
        dtype="float64",
    )


def _eligible_control_site_rows(
    control_site_mapping: ControlSiteMapping,
) -> tuple[EligibleControlSiteRow, ...]:
    rows: list[EligibleControlSiteRow] = []
    for row in control_site_mapping.row_eligibility:
        if not row.is_control:
            continue
        site_key = _require_site_key(row)
        row_position = 0 if row.row_position is None else int(row.row_position)
        rows.append(
            EligibleControlSiteRow(
                site_key=site_key,
                row_position=row_position,
                label=row.label,
                weight=row.weight,
                group=row.group,
            )
        )
    return tuple(rows)


def _resolve_observation_mask(
    *,
    missingness_policy: CorrectionMissingnessPolicy,
    control_site_mapping: ControlSiteMapping,
    sample_order: tuple[str, ...],
) -> ObservationMask:
    if missingness_policy.observation_mask is not None:
        return missingness_policy.observation_mask
    return ObservationMask(
        feature_ids=tuple(
            _require_site_key(row) for row in control_site_mapping.row_eligibility
        ),
        sample_ids=sample_order,
        originally_missing_cells=(),
    )


def _replicate_structure(
    *,
    replicate_column: str | None,
    replicate_by_sample: Mapping[str, str] | None,
    structure_diagnostics: ReplicateStructureDiagnostics | None,
    sample_order: tuple[str, ...],
) -> ReplicateStructure:
    if replicate_by_sample is None:
        return ReplicateStructure(
            replicate_column=replicate_column,
            replicate_by_sample=None,
            replicate_labels=None,
            replicate_groups={},
            structure_diagnostics=None,
        )
    labels = tuple(replicate_by_sample[sample] for sample in sample_order)
    grouped: dict[str, list[str]] = {}
    for sample, label in zip(sample_order, labels, strict=True):
        grouped.setdefault(label, []).append(sample)
    return ReplicateStructure(
        replicate_column=replicate_column,
        replicate_by_sample={
            sample: replicate_by_sample[sample] for sample in sample_order
        },
        replicate_labels=labels,
        replicate_groups={label: tuple(samples) for label, samples in grouped.items()},
        structure_diagnostics=structure_diagnostics,
    )


def _diagnostic_requirements(
    *, diagnostics_enabled: bool
) -> BatchCorrectionDiagnosticRequirements:
    required = (
        "resolved_design_matrix",
        "eligible_control_site_rows",
        "observation_mask",
        "temporary_imputation_policy",
        "stage_order",
        "provenance_seed_data",
    )
    if diagnostics_enabled:
        required = (*required, "diagnostic_tables")
    return BatchCorrectionDiagnosticRequirements(
        diagnostics_enabled=diagnostics_enabled,
        required_payloads=required,
    )


def _stage_order(stage_order: InternalBatchCorrectionStageOrder) -> tuple[str, ...]:
    if stage_order is SUPPORTED_INTERNAL_BATCH_CORRECTION_STAGE_ORDER:
        return SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER
    raise PhosPyInputError(
        "batch-correction workflow stage_order="
        f"{stage_order.value!r} is unsupported by the current dataset preprocessing "
        "pipeline; requested stage order implies "
        f"{' -> '.join(_implied_stage_order(stage_order))}; supported stage "
        "order is "
        f"{' -> '.join(SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER)}; "
        "provenance must match the actual executed pipeline"
    )


def _implied_stage_order(
    stage_order: InternalBatchCorrectionStageOrder,
) -> tuple[str, ...]:
    if stage_order is (
        InternalBatchCorrectionStageOrder.AFTER_INTENSITY_TRANSFORM_BEFORE_MISSING_DATA
    ):
        return ("intensity_transform", "batch_correction", "missing_data")
    if stage_order is (
        InternalBatchCorrectionStageOrder.AFTER_TOTAL_PROTEIN_CORRECTION_BEFORE_DOWNSTREAM
    ):
        return (
            "total_protein_correction",
            "batch_correction",
            "downstream_workflows",
        )
    return SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER


def _provenance_seed_data(
    *,
    config: InternalBatchCorrectionRequest,
    dataset_metadata: ResolvedBatchDesignMetadata,
    plan: ResolvedBatchCorrectionPlan,
    control_site_mapping: ControlSiteMapping,
    missingness_policy: CorrectionMissingnessPolicy,
) -> dict[str, object]:
    return {
        "method": config.method.value,
        "batch_column": config.batch_column,
        "condition_columns": list(config.condition_columns),
        "replicate_column": config.replicate_column,
        "replicate_metadata_role": REPLICATE_METADATA_ROLE,
        "replicate_metadata_role_description": REPLICATE_METADATA_ROLE_DESCRIPTION,
        "replicate_metadata_used_for_numerical_factor_estimation": False,
        "replicate_metadata_enables_ruv_iii_semantics": False,
        "control_site_source": config.control_site_source.value,
        "control_site_mode": config.control_site_mode.value,
        "missing_value_policy": config.missing_value_policy.value,
        "imputation_policy": config.imputation_policy.value,
        "n_unwanted_factors": config.n_unwanted_factors,
        "executed_stage_order": list(plan.stage_order),
        "requested_stage_order": config.stage_order.value,
        "stage_order_policy": config.stage_order.value,
        "sample_order": list(dataset_metadata.sample_order),
        "batch_by_sample": dict(dataset_metadata.batch_by_sample),
        "condition_by_sample": dict(dataset_metadata.condition_by_sample),
        "replicate_by_sample": (
            None
            if dataset_metadata.replicate_by_sample is None
            else dict(dataset_metadata.replicate_by_sample)
        ),
        "replicate_structure_diagnostics": (
            None
            if dataset_metadata.replicate_structure_diagnostics is None
            else dataset_metadata.replicate_structure_diagnostics.to_payload()
        ),
        "resolved_design_matrix": _frame_payload(plan.resolved_design_matrix),
        "batch_terms": list(plan.batch_terms),
        "condition_terms_to_preserve": list(plan.condition_terms_to_preserve),
        "eligible_control_site_rows": [
            row.to_payload() for row in plan.eligible_control_site_rows
        ],
        "control_site_mapping": control_site_mapping.to_payload(),
        "observation_mask": plan.observation_mask.to_payload(),
        "missingness_policy": missingness_policy.to_payload(),
        "diagnostic_requirements": plan.diagnostic_requirements.to_payload(),
    }


def _frame_payload(frame: pd.DataFrame) -> dict[str, object]:
    return {
        "index": [str(value) for value in frame.index.tolist()],
        "columns": [str(value) for value in frame.columns.tolist()],
        "data": [
            [float(value) for value in row]
            for row in frame.to_numpy(dtype="float64").tolist()
        ],
    }


def _temporary_imputation_payload(
    value: TemporaryImputationPolicy,
) -> dict[str, object]:
    return value.to_payload()


def _levels_in_order(labels: Sequence[str]) -> tuple[str, ...]:
    levels: list[str] = []
    seen: set[str] = set()
    for label in labels:
        if label in seen:
            continue
        seen.add(label)
        levels.append(label)
    return tuple(levels)


def _require_site_key(row: ControlSiteEligibility) -> str:
    site_key = row.site_key
    if site_key is None:
        raise ValueError("validated control-site mapping contains a missing site_key")
    return str(site_key)


__all__ = [
    "BatchCorrectionDiagnosticRequirements",
    "BatchCorrectionPlanInterpreter",
    "EligibleControlSiteRow",
    "REPLICATE_METADATA_ROLE",
    "REPLICATE_METADATA_ROLE_DESCRIPTION",
    "ReplicateStructure",
    "ResolvedBatchCorrectionPlan",
]

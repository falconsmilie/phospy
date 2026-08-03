"""Typed evidence validation for preprocessing quantitative-operation contracts."""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, NoReturn, cast

from phospy.errors.build import DatasetBuildError
from phospy.provenance.models import BatchCorrectionProvenance, TableFingerprint
from phospy.science.transformations.models import (
    IntensityScaleState,
    IntensityTransformationEvent,
    QuantitativeMeaningEvidenceMode,
)
from phospy.science.transformations.quantitative_contracts import (
    QuantitativeContractState,
    QuantitativeEvidenceRequirement,
    QuantitativeOperationContract,
    QuantitativeTransitionEvidence,
)

if TYPE_CHECKING:
    from phospy.science.datasets.preprocessing.stage_contract import (
        InterpretedPreprocessingStageContract,
    )

_DATASET_IMPUTATION_OBSERVATION_MASK = "dataset.imputation_observation_mask"
_DATASET_SAMPLE_METADATA = "dataset.sample_metadata"
_REPORT_ROW_AUDIT = "report.row_audit"
_SAMPLE_DESIGN_OUTPUT_TABLES = frozenset(
    {
        "dataset.comparisons",
        "report.comparison_group_stats",
        "report.comparison_pair_stats",
    }
)


@dataclass(frozen=True, slots=True)
class ObservationMaskEvidence:
    """Typed or fingerprinted observation/missingness-mask evidence."""

    mask_fingerprint: TableFingerprint

    def __post_init__(self) -> None:
        if not isinstance(self.mask_fingerprint, TableFingerprint):
            raise DatasetBuildError(
                "observation-mask evidence requires a TableFingerprint"
            )


@dataclass(frozen=True, slots=True)
class RowAuditEvidence:
    """Typed row-audit evidence for a quantitative preprocessing operation."""

    record_count: int
    row_audit_fingerprint: TableFingerprint | None = None

    def __post_init__(self) -> None:
        record_count = _non_negative_int(
            self.record_count,
            field_name="row_audit_evidence.record_count",
        )
        if self.row_audit_fingerprint is not None and not isinstance(
            self.row_audit_fingerprint,
            TableFingerprint,
        ):
            raise DatasetBuildError(
                "row-audit evidence row_audit_fingerprint must be a TableFingerprint"
            )
        object.__setattr__(self, "record_count", record_count)


@dataclass(frozen=True, slots=True)
class TotalProteinRowMappingEvidence:
    """Typed phosphosite-to-total row-mapping evidence."""

    input_phosphosite_row_count: int
    corrected_phosphosite_row_ids: tuple[str, ...]
    uncorrected_phosphosite_row_ids: tuple[str, ...]
    corrected_phosphosite_to_total_row_ids: tuple[tuple[str, str], ...]
    total_protein_row_count: int | None = None

    def __post_init__(self) -> None:
        input_row_count = _non_negative_int(
            self.input_phosphosite_row_count,
            field_name="total_protein_row_mapping.input_phosphosite_row_count",
        )
        corrected = _normalize_text_tuple(
            self.corrected_phosphosite_row_ids,
            field_name="total_protein_row_mapping.corrected_phosphosite_row_ids",
            unique=False,
        )
        uncorrected = _normalize_text_tuple(
            self.uncorrected_phosphosite_row_ids,
            field_name="total_protein_row_mapping.uncorrected_phosphosite_row_ids",
            unique=False,
        )
        mapping = _normalize_text_pair_tuple(
            self.corrected_phosphosite_to_total_row_ids,
            field_name=(
                "total_protein_row_mapping.corrected_phosphosite_to_total_row_ids"
            ),
            unique_left=False,
        )
        total_row_count = (
            None
            if self.total_protein_row_count is None
            else _non_negative_int(
                self.total_protein_row_count,
                field_name="total_protein_row_mapping.total_protein_row_count",
            )
        )
        corrected_set = set(corrected)
        uncorrected_set = set(uncorrected)
        if corrected_set & uncorrected_set:
            raise DatasetBuildError(
                "total-protein row-mapping evidence contains overlapping corrected "
                "and uncorrected phosphosite rows"
            )
        mapped_phosphosites = Counter(left for left, _right in mapping)
        if mapped_phosphosites != Counter(corrected):
            raise DatasetBuildError(
                "total-protein row-mapping evidence must map exactly every "
                "corrected phosphosite row to a total-protein row"
            )
        if len(corrected) + len(uncorrected) != input_row_count:
            raise DatasetBuildError(
                "total-protein row-mapping evidence corrected and uncorrected row "
                "counts must equal the input phosphosite row count"
            )
        object.__setattr__(self, "input_phosphosite_row_count", input_row_count)
        object.__setattr__(self, "corrected_phosphosite_row_ids", corrected)
        object.__setattr__(self, "uncorrected_phosphosite_row_ids", uncorrected)
        object.__setattr__(
            self,
            "corrected_phosphosite_to_total_row_ids",
            mapping,
        )
        object.__setattr__(self, "total_protein_row_count", total_row_count)

    @property
    def corrected_row_count(self) -> int:
        return len(self.corrected_phosphosite_row_ids)

    @property
    def uncorrected_row_count(self) -> int:
        return len(self.uncorrected_phosphosite_row_ids)


@dataclass(frozen=True, slots=True)
class SampleMetadataDesignEvidence:
    """Typed sample-metadata design evidence."""

    sample_metadata_fingerprint: TableFingerprint
    resolved_design_fingerprints: tuple[TableFingerprint, ...] = ()
    resolved_sample_count: int | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.sample_metadata_fingerprint, TableFingerprint):
            raise DatasetBuildError(
                "sample-metadata design evidence requires a sample metadata "
                "TableFingerprint"
            )
        design_fingerprints = _normalize_fingerprint_tuple(
            self.resolved_design_fingerprints,
            field_name="sample_metadata_design.resolved_design_fingerprints",
        )
        sample_count = (
            None
            if self.resolved_sample_count is None
            else _non_negative_int(
                self.resolved_sample_count,
                field_name="sample_metadata_design.resolved_sample_count",
            )
        )
        object.__setattr__(
            self,
            "resolved_design_fingerprints",
            design_fingerprints,
        )
        object.__setattr__(self, "resolved_sample_count", sample_count)


@dataclass(frozen=True, slots=True)
class ControlSiteSetEvidence:
    """Typed resolved control-site-set evidence."""

    selected_site_key_rows: tuple[str, ...]
    control_site_fingerprint: TableFingerprint | None = None

    def __post_init__(self) -> None:
        selected = _normalize_text_tuple(
            self.selected_site_key_rows,
            field_name="control_site_set.selected_site_key_rows",
        )
        if not selected:
            raise DatasetBuildError(
                "control-site-set evidence requires selected site-key rows"
            )
        if self.control_site_fingerprint is not None and not isinstance(
            self.control_site_fingerprint,
            TableFingerprint,
        ):
            raise DatasetBuildError(
                "control-site-set evidence control_site_fingerprint must be a "
                "TableFingerprint"
            )
        object.__setattr__(self, "selected_site_key_rows", selected)


@dataclass(frozen=True, slots=True)
class QuantitativeOperationEvidence:
    """Stage-supplied typed evidence sidecars for quantitative requirements."""

    observation_mask: ObservationMaskEvidence | None = None
    row_audit: RowAuditEvidence | None = None
    total_protein_row_mapping: TotalProteinRowMappingEvidence | None = None
    sample_metadata_design: SampleMetadataDesignEvidence | None = None
    control_site_set: ControlSiteSetEvidence | None = None


@dataclass(frozen=True, slots=True)
class QuantitativeOperationEvidenceContext:
    """Authoritative execution context used to validate evidence requirements."""

    stage: str
    operation: str
    quantitative_contract: QuantitativeOperationContract
    trace_record: object
    input_quantitative_state: QuantitativeContractState | None = None
    input_intensity_scale_state: IntensityScaleState | None = None
    input_quantitative_meaning_evidence_mode: QuantitativeMeaningEvidenceMode | None = (
        None
    )
    interpreted_contract: InterpretedPreprocessingStageContract | None = None
    previous_preprocessing_state: object | None = None
    current_preprocessing_state: object | None = None


class QuantitativeOperationEvidenceValidator:
    """Validate required quantitative-operation evidence by requirement enum."""

    @classmethod
    def supported_requirements(cls) -> frozenset[QuantitativeEvidenceRequirement]:
        return frozenset(_VALIDATOR_BY_REQUIREMENT)

    @classmethod
    def require_all_enum_members_handled(cls) -> None:
        observed = cls.supported_requirements()
        expected = frozenset(QuantitativeEvidenceRequirement)
        missing = expected - observed
        extra = observed - expected
        if missing or extra:
            details: list[str] = []
            if missing:
                details.append(
                    "missing validators for "
                    + ", ".join(sorted(item.value for item in missing))
                )
            if extra:
                details.append(
                    "unknown validators for "
                    + ", ".join(sorted(item.value for item in extra))
                )
            raise DatasetBuildError(
                "dataset preprocessing quantitative evidence validator is not "
                "exhaustive: " + "; ".join(details)
            )

    @classmethod
    def require_supported_contract(
        cls,
        contract: QuantitativeOperationContract,
        *,
        stage: str,
        operation: str,
    ) -> None:
        cls.require_all_enum_members_handled()
        requirements = contract.required_evidence
        unsupported = tuple(
            requirement
            for requirement in requirements
            if requirement not in _VALIDATOR_BY_REQUIREMENT
        )
        if unsupported:
            raise DatasetBuildError(
                "dataset preprocessing quantitative evidence contract declares "
                "unsupported evidence requirements: "
                + ", ".join(repr(item.value) for item in unsupported)
                + f"; stage={stage!r}, operation={operation!r}"
            )

    def validate(self, context: QuantitativeOperationEvidenceContext) -> None:
        self.require_supported_contract(
            context.quantitative_contract,
            stage=context.stage,
            operation=context.operation,
        )
        requirements = context.quantitative_contract.required_evidence
        if QuantitativeEvidenceRequirement.NONE in requirements:
            if len(requirements) != 1:
                _raise_missing(
                    context,
                    QuantitativeEvidenceRequirement.NONE,
                    expected=(
                        "QuantitativeEvidenceRequirement.NONE must be the only "
                        "declared evidence requirement"
                    ),
                    corrective_action=(
                        "remove NONE when the operation requires execution evidence"
                    ),
                )
            return

        for requirement in QuantitativeEvidenceRequirement:
            if requirement is QuantitativeEvidenceRequirement.NONE:
                continue
            if requirement not in requirements:
                continue
            _VALIDATOR_BY_REQUIREMENT[requirement](context)


def resolve_quantitative_operation_evidence(
    *,
    provided: QuantitativeOperationEvidence | None,
    consumed_input_tables: tuple[TableFingerprint, ...],
    produced_output_tables: tuple[TableFingerprint, ...],
    batch_correction_provenance: BatchCorrectionProvenance | None,
) -> QuantitativeOperationEvidence | None:
    """Merge stage-provided typed evidence with generic fingerprint evidence."""

    base = provided or QuantitativeOperationEvidence()
    observation_mask = base.observation_mask or _resolve_observation_mask_evidence(
        consumed_input_tables=consumed_input_tables,
        produced_output_tables=produced_output_tables,
        batch_correction_provenance=batch_correction_provenance,
    )
    row_audit = base.row_audit or _resolve_row_audit_evidence(
        produced_output_tables=produced_output_tables
    )
    sample_design = (
        base.sample_metadata_design
        or _resolve_sample_metadata_design_evidence(
            consumed_input_tables=consumed_input_tables,
            produced_output_tables=produced_output_tables,
            batch_correction_provenance=batch_correction_provenance,
        )
    )
    control_site_set = base.control_site_set or _resolve_control_site_set_evidence(
        batch_correction_provenance=batch_correction_provenance
    )
    resolved = QuantitativeOperationEvidence(
        observation_mask=observation_mask,
        row_audit=row_audit,
        total_protein_row_mapping=base.total_protein_row_mapping,
        sample_metadata_design=sample_design,
        control_site_set=control_site_set,
    )
    if (
        resolved.observation_mask is None
        and resolved.row_audit is None
        and resolved.total_protein_row_mapping is None
        and resolved.sample_metadata_design is None
        and resolved.control_site_set is None
    ):
        return None
    return resolved


def _require_none(_context: QuantitativeOperationEvidenceContext) -> None:
    return None


def _require_established_input_scale(
    context: QuantitativeOperationEvidenceContext,
) -> None:
    input_scale_state = context.input_intensity_scale_state
    if (
        isinstance(input_scale_state, IntensityScaleState)
        and input_scale_state.is_established
        and input_scale_state.establishment_provenance is not None
    ):
        return
    if isinstance(context.input_quantitative_state, QuantitativeContractState):
        return
    _raise_missing(
        context,
        QuantitativeEvidenceRequirement.ESTABLISHED_INPUT_SCALE,
        expected="typed established input-scale state",
        corrective_action=(
            "pass an established IntensityScaleState or typed quantitative "
            "contract state into preprocessing evidence validation"
        ),
    )


def _require_declared_or_inferred_input_meaning(
    context: QuantitativeOperationEvidenceContext,
) -> None:
    input_state = context.input_quantitative_state
    mode = context.input_quantitative_meaning_evidence_mode
    if (
        isinstance(input_state, QuantitativeContractState)
        and isinstance(mode, QuantitativeMeaningEvidenceMode)
        and mode is not QuantitativeMeaningEvidenceMode.LEGACY_UNVERIFIED
    ):
        return
    input_scale_state = context.input_intensity_scale_state
    provenance = (
        None
        if input_scale_state is None
        else input_scale_state.quantitative_meaning_provenance
    )
    if (
        isinstance(input_scale_state, IntensityScaleState)
        and input_scale_state.quantity is not None
        and provenance is not None
        and provenance.evidence_mode
        is not QuantitativeMeaningEvidenceMode.LEGACY_UNVERIFIED
    ):
        return
    _raise_missing(
        context,
        QuantitativeEvidenceRequirement.DECLARED_OR_INFERRED_INPUT_MEANING,
        expected="typed quantitative meaning with recognised evidence mode",
        corrective_action=(
            "establish quantitative meaning through caller declaration, scale "
            "contract inference, or a typed PhosPy operation transition before "
            "executing this stage"
        ),
    )


def _require_typed_intensity_transformation_event(
    context: QuantitativeOperationEvidenceContext,
) -> None:
    event = _trace_attr(context.trace_record, "intensity_transformation_event")
    if not isinstance(event, IntensityTransformationEvent):
        _raise_missing(
            context,
            QuantitativeEvidenceRequirement.TYPED_INTENSITY_TRANSFORMATION_EVENT,
            expected="valid IntensityTransformationEvent for this execution",
            corrective_action=(
                "return a typed IntensityTransformationEvent from the stage result"
            ),
        )
    input_hash = _trace_attr(context.trace_record, "phospho_input_hash")
    output_hash = _trace_attr(context.trace_record, "phospho_output_hash")
    if (
        isinstance(input_hash, str)
        and input_hash
        and isinstance(output_hash, str)
        and output_hash
        and event.input_fingerprint == input_hash
        and event.output_fingerprint == output_hash
    ):
        return
    _raise_missing(
        context,
        QuantitativeEvidenceRequirement.TYPED_INTENSITY_TRANSFORMATION_EVENT,
        expected=(
            "IntensityTransformationEvent with input/output fingerprints matching "
            "this execution trace"
        ),
        corrective_action=(
            "attach the observed transformation event after input and output table "
            "fingerprints are resolved"
        ),
    )


def _require_table_fingerprints(context: QuantitativeOperationEvidenceContext) -> None:
    contract = context.quantitative_contract
    if not contract.provenance_input_tables or contract.provenance_output_table is None:
        _raise_missing(
            context,
            QuantitativeEvidenceRequirement.TABLE_FINGERPRINTS,
            expected=(
                "contract-declared provenance input tables and provenance output table"
            ),
            corrective_action=(
                "declare provenance_input_tables and provenance_output_table on "
                "the quantitative operation contract"
            ),
        )
    consumed = _fingerprint_names(_consumed_input_tables(context.trace_record))
    produced = _fingerprint_names(_produced_output_tables(context.trace_record))
    missing_inputs = tuple(
        table for table in contract.provenance_input_tables if table not in consumed
    )
    missing_output = contract.provenance_output_table not in produced
    if not missing_inputs and not missing_output:
        return
    details: list[str] = []
    if missing_inputs:
        details.append(
            "input fingerprints missing for "
            + ", ".join(repr(item) for item in missing_inputs)
        )
    if missing_output:
        details.append(
            f"output fingerprint missing for {contract.provenance_output_table!r}"
        )
    _raise_missing(
        context,
        QuantitativeEvidenceRequirement.TABLE_FINGERPRINTS,
        expected="; ".join(details),
        corrective_action=(
            "include the required tables in the stage metadata and preserve their "
            "TableFingerprint records in the execution trace"
        ),
    )


def _require_total_protein_row_mapping(
    context: QuantitativeOperationEvidenceContext,
) -> None:
    evidence = _operation_evidence(context.trace_record)
    mapping = None if evidence is None else evidence.total_protein_row_mapping
    transition_evidence = _trace_attr(
        context.trace_record,
        "quantitative_transition_evidence",
    )
    if not isinstance(mapping, TotalProteinRowMappingEvidence):
        _raise_missing(
            context,
            QuantitativeEvidenceRequirement.TOTAL_PROTEIN_ROW_MAPPING,
            expected="typed corrected/uncorrected total-protein row-mapping evidence",
            corrective_action=(
                "record TotalProteinRowMappingEvidence with corrected rows, "
                "uncorrected rows, and phosphosite-to-total row links"
            ),
        )
    input_rows = _trace_attr(context.trace_record, "input_rows")
    output_rows = _trace_attr(context.trace_record, "output_rows")
    if (
        isinstance(input_rows, int)
        and mapping.input_phosphosite_row_count != input_rows
    ):
        _raise_missing(
            context,
            QuantitativeEvidenceRequirement.TOTAL_PROTEIN_ROW_MAPPING,
            expected="row-mapping input row count matching the execution input rows",
            corrective_action=(
                "record total-protein row-mapping evidence from the executed "
                "input phosphosite rows"
            ),
        )
    if (
        isinstance(output_rows, int)
        and mapping.input_phosphosite_row_count != output_rows
    ):
        _raise_missing(
            context,
            QuantitativeEvidenceRequirement.TOTAL_PROTEIN_ROW_MAPPING,
            expected="row-mapping row count coherent with the execution output rows",
            corrective_action=(
                "record total-protein row mapping for the exact rows retained by "
                "the correction execution"
            ),
        )
    if isinstance(transition_evidence, QuantitativeTransitionEvidence):
        if (
            transition_evidence.total_protein_corrected_row_count
            != mapping.corrected_row_count
            or transition_evidence.total_protein_uncorrected_row_count
            != mapping.uncorrected_row_count
        ):
            _raise_missing(
                context,
                QuantitativeEvidenceRequirement.TOTAL_PROTEIN_ROW_MAPPING,
                expected=(
                    "typed row-mapping counts matching QuantitativeTransitionEvidence"
                ),
                corrective_action=(
                    "record corrected and uncorrected row counts from the same "
                    "executed row mapping used for the semantic transition"
                ),
            )
        return
    _raise_missing(
        context,
        QuantitativeEvidenceRequirement.TOTAL_PROTEIN_ROW_MAPPING,
        expected="QuantitativeTransitionEvidence paired with row-mapping evidence",
        corrective_action=(
            "return QuantitativeTransitionEvidence and TotalProteinRowMappingEvidence "
            "from total-protein correction execution"
        ),
    )


def _require_missingness_mask(context: QuantitativeOperationEvidenceContext) -> None:
    evidence = _operation_evidence(context.trace_record)
    if evidence is not None and isinstance(
        evidence.observation_mask, ObservationMaskEvidence
    ):
        return
    if _find_fingerprint(
        (
            *_consumed_input_tables(context.trace_record),
            *_produced_output_tables(context.trace_record),
        ),
        name=_DATASET_IMPUTATION_OBSERVATION_MASK,
    ):
        return
    provenance = _batch_correction_provenance(context.trace_record)
    if provenance is not None and provenance.observation_masks:
        return
    _raise_missing(
        context,
        QuantitativeEvidenceRequirement.MISSINGNESS_MASK,
        expected=(
            "typed or fingerprinted observation/missingness mask associated with "
            "this execution"
        ),
        corrective_action=(
            "retain the observation mask as a typed evidence object, a table "
            "fingerprint, or BatchCorrectionProvenance observation-mask fingerprint"
        ),
    )


def _require_sample_metadata_design(
    context: QuantitativeOperationEvidenceContext,
) -> None:
    evidence = _operation_evidence(context.trace_record)
    if evidence is not None and isinstance(
        evidence.sample_metadata_design, SampleMetadataDesignEvidence
    ):
        return
    sample_metadata_fingerprint = _find_fingerprint(
        _consumed_input_tables(context.trace_record),
        name=_DATASET_SAMPLE_METADATA,
    )
    if sample_metadata_fingerprint is None:
        _raise_missing(
            context,
            QuantitativeEvidenceRequirement.SAMPLE_METADATA_DESIGN,
            expected="fingerprinted sample metadata consumed by this execution",
            corrective_action=(
                "include dataset.sample_metadata in consumed input tables and "
                "preserve its TableFingerprint"
            ),
        )
    produced_design = tuple(
        fingerprint
        for fingerprint in _produced_output_tables(context.trace_record)
        if fingerprint.name in _SAMPLE_DESIGN_OUTPUT_TABLES
    )
    provenance = _batch_correction_provenance(context.trace_record)
    if produced_design:
        return
    if (
        provenance is not None
        and bool(provenance.design_metadata)
        and bool(provenance.batch_metadata)
    ):
        return
    _raise_missing(
        context,
        QuantitativeEvidenceRequirement.SAMPLE_METADATA_DESIGN,
        expected=("fingerprinted sample metadata plus typed resolved design evidence"),
        corrective_action=(
            "record resolved design outputs or typed BatchCorrectionProvenance "
            "derived from the consumed sample metadata"
        ),
    )


def _require_control_site_set(context: QuantitativeOperationEvidenceContext) -> None:
    evidence = _operation_evidence(context.trace_record)
    if evidence is not None and isinstance(
        evidence.control_site_set, ControlSiteSetEvidence
    ):
        return
    provenance = _batch_correction_provenance(context.trace_record)
    if (
        provenance is not None
        and bool(provenance.selected_site_key_rows)
        and bool(provenance.control_site_source)
    ):
        return
    _raise_missing(
        context,
        QuantitativeEvidenceRequirement.CONTROL_SITE_SET,
        expected="typed or fingerprinted resolved control-site-set evidence",
        corrective_action=(
            "record ControlSiteSetEvidence or BatchCorrectionProvenance with "
            "selected control site rows and source metadata"
        ),
    )


def _require_row_audit(context: QuantitativeOperationEvidenceContext) -> None:
    evidence = _operation_evidence(context.trace_record)
    if evidence is not None and isinstance(evidence.row_audit, RowAuditEvidence):
        if evidence.row_audit.row_audit_fingerprint is not None:
            return
        if (
            evidence.row_audit.record_count == 0
            and _trace_attr(context.trace_record, "dropped_row_count") == 0
            and _trace_attr(context.trace_record, "input_rows")
            == _trace_attr(context.trace_record, "output_rows")
        ):
            return
    if _find_fingerprint(
        _produced_output_tables(context.trace_record),
        name=_REPORT_ROW_AUDIT,
    ):
        return
    _raise_missing(
        context,
        QuantitativeEvidenceRequirement.ROW_AUDIT,
        expected="typed row-audit output or row-audit table fingerprint",
        corrective_action=(
            "emit row-audit report rows, preserve the row-audit table fingerprint, "
            "or record typed zero-row audit evidence when no row action occurred"
        ),
    )


def _require_random_seed(context: QuantitativeOperationEvidenceContext) -> None:
    random_seed = _trace_attr(context.trace_record, "random_seed")
    if isinstance(random_seed, int) and not isinstance(random_seed, bool):
        return
    _raise_missing(
        context,
        QuantitativeEvidenceRequirement.RANDOM_SEED,
        expected="concrete recorded execution seed",
        corrective_action=(
            "record the resolved random seed in typed execution diagnostics so it "
            "is copied to the accepted trace"
        ),
    )


_RequirementValidator = Callable[[QuantitativeOperationEvidenceContext], None]

_VALIDATOR_BY_REQUIREMENT: Mapping[
    QuantitativeEvidenceRequirement,
    _RequirementValidator,
] = {
    QuantitativeEvidenceRequirement.NONE: _require_none,
    QuantitativeEvidenceRequirement.ESTABLISHED_INPUT_SCALE: (
        _require_established_input_scale
    ),
    QuantitativeEvidenceRequirement.DECLARED_OR_INFERRED_INPUT_MEANING: (
        _require_declared_or_inferred_input_meaning
    ),
    QuantitativeEvidenceRequirement.TYPED_INTENSITY_TRANSFORMATION_EVENT: (
        _require_typed_intensity_transformation_event
    ),
    QuantitativeEvidenceRequirement.TABLE_FINGERPRINTS: _require_table_fingerprints,
    QuantitativeEvidenceRequirement.TOTAL_PROTEIN_ROW_MAPPING: (
        _require_total_protein_row_mapping
    ),
    QuantitativeEvidenceRequirement.MISSINGNESS_MASK: _require_missingness_mask,
    QuantitativeEvidenceRequirement.SAMPLE_METADATA_DESIGN: (
        _require_sample_metadata_design
    ),
    QuantitativeEvidenceRequirement.CONTROL_SITE_SET: _require_control_site_set,
    QuantitativeEvidenceRequirement.ROW_AUDIT: _require_row_audit,
    QuantitativeEvidenceRequirement.RANDOM_SEED: _require_random_seed,
}


def _resolve_observation_mask_evidence(
    *,
    consumed_input_tables: tuple[TableFingerprint, ...],
    produced_output_tables: tuple[TableFingerprint, ...],
    batch_correction_provenance: BatchCorrectionProvenance | None,
) -> ObservationMaskEvidence | None:
    fingerprint = _find_fingerprint(
        (*produced_output_tables, *consumed_input_tables),
        name=_DATASET_IMPUTATION_OBSERVATION_MASK,
    )
    if fingerprint is not None:
        return ObservationMaskEvidence(mask_fingerprint=fingerprint)
    if batch_correction_provenance is not None:
        for fingerprint in batch_correction_provenance.observation_masks:
            return ObservationMaskEvidence(mask_fingerprint=fingerprint)
    return None


def _resolve_row_audit_evidence(
    *,
    produced_output_tables: tuple[TableFingerprint, ...],
) -> RowAuditEvidence | None:
    fingerprint = _find_fingerprint(
        produced_output_tables,
        name=_REPORT_ROW_AUDIT,
    )
    if fingerprint is None:
        return None
    return RowAuditEvidence(
        record_count=int(fingerprint.rows),
        row_audit_fingerprint=fingerprint,
    )


def _resolve_sample_metadata_design_evidence(
    *,
    consumed_input_tables: tuple[TableFingerprint, ...],
    produced_output_tables: tuple[TableFingerprint, ...],
    batch_correction_provenance: BatchCorrectionProvenance | None,
) -> SampleMetadataDesignEvidence | None:
    sample_metadata_fingerprint = _find_fingerprint(
        consumed_input_tables,
        name=_DATASET_SAMPLE_METADATA,
    )
    if sample_metadata_fingerprint is None:
        return None
    design_fingerprints = tuple(
        fingerprint
        for fingerprint in produced_output_tables
        if fingerprint.name in _SAMPLE_DESIGN_OUTPUT_TABLES
    )
    if design_fingerprints:
        return SampleMetadataDesignEvidence(
            sample_metadata_fingerprint=sample_metadata_fingerprint,
            resolved_design_fingerprints=design_fingerprints,
            resolved_sample_count=int(sample_metadata_fingerprint.rows),
        )
    if batch_correction_provenance is not None and (
        batch_correction_provenance.design_metadata
        or batch_correction_provenance.batch_metadata
    ):
        return SampleMetadataDesignEvidence(
            sample_metadata_fingerprint=sample_metadata_fingerprint,
            resolved_design_fingerprints=(),
            resolved_sample_count=int(sample_metadata_fingerprint.rows),
        )
    return None


def _resolve_control_site_set_evidence(
    *,
    batch_correction_provenance: BatchCorrectionProvenance | None,
) -> ControlSiteSetEvidence | None:
    if batch_correction_provenance is None:
        return None
    selected = tuple(
        str(item).strip()
        for item in batch_correction_provenance.selected_site_key_rows
        if str(item).strip()
    )
    if not selected:
        return None
    return ControlSiteSetEvidence(selected_site_key_rows=selected)


def _raise_missing(
    context: QuantitativeOperationEvidenceContext,
    requirement: QuantitativeEvidenceRequirement,
    *,
    expected: str,
    corrective_action: str,
) -> NoReturn:
    raise DatasetBuildError(
        "quantitative operation evidence contract failed: "
        f"stage={context.stage!r}, operation={context.operation!r}, "
        f"missing_requirement={requirement.value!r}; "
        f"{expected}; corrective_action={corrective_action}"
    )


def _operation_evidence(trace_record: object) -> QuantitativeOperationEvidence | None:
    evidence = _trace_attr(trace_record, "quantitative_evidence")
    if isinstance(evidence, QuantitativeOperationEvidence):
        return evidence
    return None


def _batch_correction_provenance(
    trace_record: object,
) -> BatchCorrectionProvenance | None:
    provenance = _trace_attr(trace_record, "batch_correction_provenance")
    if isinstance(provenance, BatchCorrectionProvenance):
        return provenance
    return None


def _consumed_input_tables(trace_record: object) -> tuple[TableFingerprint, ...]:
    return _fingerprint_tuple(_trace_attr(trace_record, "consumed_input_tables"))


def _produced_output_tables(trace_record: object) -> tuple[TableFingerprint, ...]:
    return _fingerprint_tuple(_trace_attr(trace_record, "produced_output_tables"))


def _fingerprint_tuple(value: object) -> tuple[TableFingerprint, ...]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        return ()
    return tuple(item for item in value if isinstance(item, TableFingerprint))


def _find_fingerprint(
    fingerprints: Sequence[TableFingerprint],
    *,
    name: str,
) -> TableFingerprint | None:
    for fingerprint in fingerprints:
        if fingerprint.name == name:
            return fingerprint
    return None


def _fingerprint_names(fingerprints: Sequence[TableFingerprint]) -> frozenset[str]:
    return frozenset(fingerprint.name for fingerprint in fingerprints)


def _normalize_fingerprint_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[TableFingerprint, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise DatasetBuildError(f"{field_name} must be a sequence")
    fingerprints = tuple(values)
    for fingerprint in fingerprints:
        if not isinstance(fingerprint, TableFingerprint):
            raise DatasetBuildError(
                f"{field_name} must contain only TableFingerprint values"
            )
    return fingerprints


def _normalize_text_tuple(
    values: object,
    *,
    field_name: str,
    unique: bool = True,
) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise DatasetBuildError(f"{field_name} must be a sequence")
    normalized: list[str] = []
    seen: set[str] = set()
    for item in values:
        text = str(item).strip()
        if not text:
            raise DatasetBuildError(f"{field_name} entries must be non-empty strings")
        if unique and text in seen:
            raise DatasetBuildError(f"{field_name} entries must be unique")
        seen.add(text)
        normalized.append(text)
    return tuple(normalized)


def _normalize_text_pair_tuple(
    values: object,
    *,
    field_name: str,
    unique_left: bool = True,
) -> tuple[tuple[str, str], ...]:
    if isinstance(values, Mapping):
        raw_pairs = tuple(cast(Mapping[object, object], values).items())
    elif isinstance(values, Sequence) and not isinstance(
        values,
        (str, bytes, bytearray),
    ):
        raw_pairs = tuple(values)
    else:
        raise DatasetBuildError(f"{field_name} must be a mapping or sequence")
    pairs: list[tuple[str, str]] = []
    seen_left: set[str] = set()
    for pair in raw_pairs:
        if (
            isinstance(pair, (str, bytes, bytearray))
            or not isinstance(pair, Sequence)
            or len(pair) != 2
        ):
            raise DatasetBuildError(f"{field_name} entries must be pairs")
        left = str(pair[0]).strip()
        right = str(pair[1]).strip()
        if not left or not right:
            raise DatasetBuildError(f"{field_name} entries must be non-empty strings")
        if unique_left and left in seen_left:
            raise DatasetBuildError(
                f"{field_name} may contain each phosphosite row at most once"
            )
        seen_left.add(left)
        pairs.append((left, right))
    return tuple(sorted(pairs))


def _non_negative_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise DatasetBuildError(f"{field_name} must be an int")
    if value < 0:
        raise DatasetBuildError(f"{field_name} must be >= 0")
    return int(value)


def _trace_attr(trace_record: object, name: str) -> object:
    return getattr(trace_record, name, None)


__all__ = [
    "ControlSiteSetEvidence",
    "ObservationMaskEvidence",
    "QuantitativeOperationEvidence",
    "QuantitativeOperationEvidenceContext",
    "QuantitativeOperationEvidenceValidator",
    "RowAuditEvidence",
    "SampleMetadataDesignEvidence",
    "TotalProteinRowMappingEvidence",
    "resolve_quantitative_operation_evidence",
]

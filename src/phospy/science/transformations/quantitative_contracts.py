"""Typed quantitative-operation contracts for preprocessing semantics."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from enum import Enum
from typing import Final, cast

from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import InvalidTransformationStateError
from phospy.science.transformations.models import (
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
    QuantitativeMeaningEvidenceMode,
    default_quantitative_meaning_for_scale_kind,
)


class QuantitativeScaleTransitionKind(str, Enum):
    """Supported scale-kind transitions for quantitative preprocessing."""

    PRESERVE = "preserve"
    LINEAR_TO_LOG2 = "linear_to_log2"


class QuantitativeMeaningTransitionKind(str, Enum):
    """Supported quantitative-meaning transitions."""

    PRESERVE = "preserve"
    STATIC_MAPPING = "static_mapping"
    EVIDENCE_RESOLVED = "evidence_resolved"


class NegativeDomainPolicy(str, Enum):
    """How an operation treats negative or non-positive numeric domains."""

    NOT_APPLICABLE = "not_applicable"
    PRESERVES_INPUT_DOMAIN = "preserves_input_domain"
    REQUIRES_POSITIVE_INPUT_AFTER_PSEUDOCOUNT = (
        "requires_positive_input_after_pseudocount"
    )
    REQUIRES_LOG2_DOMAIN = "requires_log2_domain"
    MAY_INTRODUCE_NEGATIVE_VALUES = "may_introduce_negative_values"
    ALLOWS_NEGATIVE_OUTPUT = "allows_negative_output"


class QuantitativeEvidenceRequirement(str, Enum):
    """Typed evidence a quantitative operation must provide or consume."""

    NONE = "none"
    ESTABLISHED_INPUT_SCALE = "established_input_scale"
    DECLARED_OR_INFERRED_INPUT_MEANING = "declared_or_inferred_input_meaning"
    TYPED_INTENSITY_TRANSFORMATION_EVENT = "typed_intensity_transformation_event"
    TABLE_FINGERPRINTS = "table_fingerprints"
    TOTAL_PROTEIN_ROW_MAPPING = "total_protein_row_mapping"
    MISSINGNESS_MASK = "missingness_mask"
    SAMPLE_METADATA_DESIGN = "sample_metadata_design"
    CONTROL_SITE_SET = "control_site_set"
    ROW_AUDIT = "row_audit"
    RANDOM_SEED = "random_seed"


class QuantitativeReversibilityKind(str, Enum):
    """Whether an operation can be reversed from retained state."""

    REVERSIBLE = "reversible"
    CONDITIONALLY_REVERSIBLE = "conditionally_reversible"
    IRREVERSIBLE = "irreversible"
    NOT_APPLICABLE = "not_applicable"


class QuantitativeInformationLossKind(str, Enum):
    """Information-loss category introduced by a quantitative operation."""

    NONE = "none"
    ROW_FILTERING = "row_filtering"
    IMPUTATION = "imputation"
    SCALE_TRANSFORMATION = "scale_transformation"
    DISTRIBUTION_RESHAPING = "distribution_reshaping"
    ADDITIVE_RESIDUALIZATION = "additive_residualization"
    RATIO_TRANSFORMATION = "ratio_transformation"
    COMPARISON_DERIVATION = "comparison_derivation"


ALL_INTENSITY_SCALE_KINDS: Final[frozenset[IntensityScaleKind]] = frozenset(
    IntensityScaleKind
)
ALL_QUANTITATIVE_MEANINGS: Final[frozenset[QuantitativeMeaning]] = frozenset(
    QuantitativeMeaning
)

_EvidenceMeaningResolver = Callable[
    ["QuantitativeTransitionEvidence | None"],
    QuantitativeMeaning,
]
_CaveatResolver = Callable[
    [QuantitativeMeaning, "QuantitativeTransitionEvidence | None"],
    tuple[str, ...],
]


def _empty_payload() -> dict[str, object]:
    return {}


def _no_caveats(
    _target: QuantitativeMeaning,
    _evidence: QuantitativeTransitionEvidence | None,
) -> tuple[str, ...]:
    return ()


@dataclass(frozen=True, slots=True)
class QuantitativeTransitionEvidence:
    """Typed operation evidence used by semantic contract transitions."""

    total_protein_corrected_row_count: int | None = None
    total_protein_uncorrected_row_count: int | None = None
    metadata: Mapping[str, object] = field(default_factory=_empty_payload)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "total_protein_corrected_row_count",
            _optional_non_negative_int(
                self.total_protein_corrected_row_count,
                field_name=(
                    "quantitative_transition_evidence.total_protein_corrected_row_count"
                ),
            ),
        )
        object.__setattr__(
            self,
            "total_protein_uncorrected_row_count",
            _optional_non_negative_int(
                self.total_protein_uncorrected_row_count,
                field_name=(
                    "quantitative_transition_evidence."
                    "total_protein_uncorrected_row_count"
                ),
            ),
        )
        object.__setattr__(
            self,
            "metadata",
            _normalize_evidence_metadata(self.metadata),
        )

    def to_payload(self) -> dict[str, object]:
        """Return JSON-safe evidence payload for provenance parameters."""

        payload: dict[str, object] = {}
        if self.total_protein_corrected_row_count is not None:
            payload["total_protein_corrected_row_count"] = int(
                self.total_protein_corrected_row_count
            )
        if self.total_protein_uncorrected_row_count is not None:
            payload["total_protein_uncorrected_row_count"] = int(
                self.total_protein_uncorrected_row_count
            )
        payload.update(dict(self.metadata))
        return payload


@dataclass(frozen=True, slots=True)
class QuantitativeContractState:
    """Scale and meaning snapshot used while folding operation contracts."""

    scale_kind: IntensityScaleKind
    meaning: QuantitativeMeaning

    def __post_init__(self) -> None:
        scale_kind = _normalize_scale_kind(self.scale_kind)
        meaning = _normalize_meaning(self.meaning)
        _validate_meaning_scale_coherence(
            scale_kind=scale_kind,
            meaning=meaning,
        )
        object.__setattr__(self, "scale_kind", scale_kind)
        object.__setattr__(self, "meaning", meaning)


@dataclass(frozen=True, slots=True)
class QuantitativeScaleTransition:
    """Typed output scale transition."""

    kind: QuantitativeScaleTransitionKind
    output_scale_by_input: Mapping[IntensityScaleKind, IntensityScaleKind]
    output_scale_label: str

    def __post_init__(self) -> None:
        kind = _normalize_scale_transition_kind(self.kind)
        output_scale_by_input = {
            _normalize_scale_kind(source): _normalize_scale_kind(target)
            for source, target in self.output_scale_by_input.items()
        }
        if not output_scale_by_input:
            raise InvalidTransformationStateError(
                "quantitative scale transition requires at least one input scale"
            )
        label = str(self.output_scale_label).strip()
        if not label:
            raise InvalidTransformationStateError(
                "quantitative scale transition requires output_scale_label"
            )
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "output_scale_by_input", output_scale_by_input)
        object.__setattr__(self, "output_scale_label", label)

    def apply(
        self,
        input_scale_kind: IntensityScaleKind,
        *,
        stage: str,
        operation: str,
    ) -> IntensityScaleKind:
        """Return output scale kind or raise an unsupported-transition error."""

        scale_kind = _normalize_scale_kind(input_scale_kind)
        output_kind = self.output_scale_by_input.get(scale_kind)
        if output_kind is not None:
            return output_kind
        supported = ", ".join(sorted(item.value for item in self.output_scale_by_input))
        raise PhosPyInputError(
            "dataset preprocessing quantitative contract rejected unsupported "
            f"scale transition before numerical execution: stage={stage!r}, "
            f"operation={operation!r}, input_scale={scale_kind.value!r}, "
            f"accepted_input_scales=[{supported}]"
        )


@dataclass(frozen=True, slots=True)
class QuantitativeMeaningTransition:
    """Typed output quantitative-meaning transition."""

    kind: QuantitativeMeaningTransitionKind
    output_meaning_by_input: Mapping[QuantitativeMeaning, QuantitativeMeaning]
    evidence_resolver: _EvidenceMeaningResolver | None = None
    default_output_meaning: QuantitativeMeaning | None = None

    def __post_init__(self) -> None:
        kind = _normalize_meaning_transition_kind(self.kind)
        output_meaning_by_input = {
            _normalize_meaning(source): _normalize_meaning(target)
            for source, target in self.output_meaning_by_input.items()
        }
        if not output_meaning_by_input:
            raise InvalidTransformationStateError(
                "quantitative meaning transition requires at least one input meaning"
            )
        default_output_meaning = (
            None
            if self.default_output_meaning is None
            else _normalize_meaning(self.default_output_meaning)
        )
        if kind is QuantitativeMeaningTransitionKind.EVIDENCE_RESOLVED:
            if self.evidence_resolver is None:
                raise InvalidTransformationStateError(
                    "evidence-resolved quantitative meaning transitions require "
                    "an evidence resolver"
                )
            if default_output_meaning is None:
                raise InvalidTransformationStateError(
                    "evidence-resolved quantitative meaning transitions require "
                    "a default output meaning for pre-execution validation"
                )
        elif self.evidence_resolver is not None:
            raise InvalidTransformationStateError(
                "only evidence-resolved quantitative meaning transitions may "
                "declare an evidence resolver"
            )
        object.__setattr__(self, "kind", kind)
        object.__setattr__(
            self,
            "output_meaning_by_input",
            output_meaning_by_input,
        )
        object.__setattr__(self, "default_output_meaning", default_output_meaning)

    def apply(
        self,
        input_meaning: QuantitativeMeaning,
        *,
        evidence: QuantitativeTransitionEvidence | None,
        stage: str,
        operation: str,
    ) -> QuantitativeMeaning:
        """Return output quantitative meaning or reject unsupported input."""

        meaning = _normalize_meaning(input_meaning)
        mapped = self.output_meaning_by_input.get(meaning)
        if mapped is None:
            supported = ", ".join(
                sorted(item.value for item in self.output_meaning_by_input)
            )
            raise PhosPyInputError(
                "dataset preprocessing quantitative contract rejected unsupported "
                "quantitative-meaning transition before numerical execution: "
                f"stage={stage!r}, operation={operation!r}, "
                f"input_meaning={meaning.value!r}, "
                f"accepted_input_meanings=[{supported}]"
            )
        if self.kind is not QuantitativeMeaningTransitionKind.EVIDENCE_RESOLVED:
            return mapped
        resolver = self.evidence_resolver
        if resolver is None:  # pragma: no cover - guarded by __post_init__.
            raise InvalidTransformationStateError(
                "evidence-resolved quantitative meaning transition missing resolver"
            )
        if evidence is None:
            default = self.default_output_meaning
            if default is None:  # pragma: no cover - guarded by __post_init__.
                raise InvalidTransformationStateError(
                    "evidence-resolved quantitative meaning transition missing "
                    "default output meaning"
                )
            return default
        return _normalize_meaning(resolver(evidence))


@dataclass(frozen=True, slots=True)
class QuantitativeOperationContract:
    """Operation-level quantitative scale and meaning contract."""

    accepted_input_scale_kinds: frozenset[IntensityScaleKind]
    accepted_quantitative_meanings: frozenset[QuantitativeMeaning]
    output_scale_transition: QuantitativeScaleTransition
    output_meaning_transition: QuantitativeMeaningTransition
    preserves_abundance: bool
    negative_domain_policy: NegativeDomainPolicy
    required_evidence: frozenset[QuantitativeEvidenceRequirement]
    reversibility: QuantitativeReversibilityKind
    information_loss: QuantitativeInformationLossKind
    operation_id: str | None = None
    producer_id: str | None = None
    evidence_mode: QuantitativeMeaningEvidenceMode = (
        QuantitativeMeaningEvidenceMode.DERIVED_BY_PHOSPY_OPERATION
    )
    provenance_input_tables: tuple[str, ...] = ()
    provenance_output_table: str | None = None
    provenance_diagnostic_fields: tuple[str, ...] = ()
    emits_state_transition_event: bool = True
    caveat_resolver: _CaveatResolver = _no_caveats

    def __post_init__(self) -> None:
        accepted_scales = frozenset(
            _normalize_scale_kind(item) for item in self.accepted_input_scale_kinds
        )
        accepted_meanings = frozenset(
            _normalize_meaning(item) for item in self.accepted_quantitative_meanings
        )
        if not accepted_scales:
            raise InvalidTransformationStateError(
                "quantitative operation contract requires accepted input scales"
            )
        if not accepted_meanings:
            raise InvalidTransformationStateError(
                "quantitative operation contract requires accepted meanings"
            )
        _require_mapping_covers(
            observed=frozenset(self.output_scale_transition.output_scale_by_input),
            required=accepted_scales,
            field_name="output_scale_transition",
        )
        _require_mapping_covers(
            observed=frozenset(self.output_meaning_transition.output_meaning_by_input),
            required=accepted_meanings,
            field_name="output_meaning_transition",
        )
        required_evidence = frozenset(
            _normalize_evidence_requirement(item) for item in self.required_evidence
        )
        if not required_evidence:
            raise InvalidTransformationStateError(
                "quantitative operation contract requires required_evidence; use "
                "QuantitativeEvidenceRequirement.NONE for no external evidence"
            )
        negative_domain_policy = _normalize_negative_domain_policy(
            self.negative_domain_policy
        )
        reversibility = _normalize_reversibility(self.reversibility)
        information_loss = _normalize_information_loss(self.information_loss)
        evidence_mode = _normalize_evidence_mode(self.evidence_mode)
        operation_id = _normalize_optional_identifier(
            self.operation_id,
            field_name="quantitative_operation_contract.operation_id",
        )
        producer_id = _normalize_optional_identifier(
            self.producer_id,
            field_name="quantitative_operation_contract.producer_id",
        )
        provenance_output_table = (
            None
            if self.provenance_output_table is None
            else _normalize_required_text(
                self.provenance_output_table,
                field_name="quantitative_operation_contract.provenance_output_table",
            )
        )
        provenance_input_tables = tuple(
            _normalize_required_text(
                item,
                field_name="quantitative_operation_contract.provenance_input_tables",
            )
            for item in self.provenance_input_tables
        )
        provenance_diagnostic_fields = tuple(
            _normalize_required_text(
                item,
                field_name=(
                    "quantitative_operation_contract.provenance_diagnostic_fields"
                ),
            )
            for item in self.provenance_diagnostic_fields
        )
        if self.emits_quantitative_meaning_state_event and (
            operation_id is None
            or producer_id is None
            or provenance_output_table is None
        ):
            raise InvalidTransformationStateError(
                "quantitative operation contracts that can change meaning require "
                "operation_id, producer_id, and provenance_output_table"
            )
        object.__setattr__(self, "accepted_input_scale_kinds", accepted_scales)
        object.__setattr__(self, "accepted_quantitative_meanings", accepted_meanings)
        object.__setattr__(self, "required_evidence", required_evidence)
        object.__setattr__(self, "negative_domain_policy", negative_domain_policy)
        object.__setattr__(self, "reversibility", reversibility)
        object.__setattr__(self, "information_loss", information_loss)
        object.__setattr__(self, "evidence_mode", evidence_mode)
        object.__setattr__(self, "operation_id", operation_id)
        object.__setattr__(self, "producer_id", producer_id)
        object.__setattr__(self, "provenance_output_table", provenance_output_table)
        object.__setattr__(self, "provenance_input_tables", provenance_input_tables)
        object.__setattr__(
            self,
            "provenance_diagnostic_fields",
            provenance_diagnostic_fields,
        )

    @property
    def emits_meaning_transition(self) -> bool:
        """Return whether this contract can produce a new meaning value."""

        for (
            source,
            target,
        ) in self.output_meaning_transition.output_meaning_by_input.items():
            if source is not target:
                return True
        return (
            self.output_meaning_transition.kind
            is QuantitativeMeaningTransitionKind.EVIDENCE_RESOLVED
        )

    @property
    def emits_quantitative_meaning_state_event(self) -> bool:
        """Return whether state-builder should mint a semantic transition event."""

        return bool(self.emits_state_transition_event and self.emits_meaning_transition)

    def validate_and_transition(
        self,
        state: QuantitativeContractState,
        *,
        stage: str,
        operation: str,
        evidence: QuantitativeTransitionEvidence | None = None,
    ) -> QuantitativeContractState:
        """Validate the input state and return the contract output state."""

        if state.scale_kind not in self.accepted_input_scale_kinds:
            supported = ", ".join(
                sorted(item.value for item in self.accepted_input_scale_kinds)
            )
            raise PhosPyInputError(
                "dataset preprocessing quantitative contract rejected unsupported "
                f"input scale before numerical execution: stage={stage!r}, "
                f"operation={operation!r}, input_scale={state.scale_kind.value!r}, "
                f"accepted_input_scales=[{supported}]"
            )
        if state.meaning not in self.accepted_quantitative_meanings:
            supported = ", ".join(
                sorted(item.value for item in self.accepted_quantitative_meanings)
            )
            raise PhosPyInputError(
                "dataset preprocessing quantitative contract rejected unsupported "
                "input quantitative meaning before numerical execution: "
                f"stage={stage!r}, operation={operation!r}, "
                f"input_meaning={state.meaning.value!r}, "
                f"accepted_input_meanings=[{supported}]"
            )
        output_scale_kind = self.output_scale_transition.apply(
            state.scale_kind,
            stage=stage,
            operation=operation,
        )
        output_meaning = self.output_meaning_transition.apply(
            state.meaning,
            evidence=evidence,
            stage=stage,
            operation=operation,
        )
        return QuantitativeContractState(
            scale_kind=output_scale_kind,
            meaning=output_meaning,
        )

    def caveat_codes(
        self,
        *,
        target: QuantitativeMeaning,
        evidence: QuantitativeTransitionEvidence | None,
    ) -> tuple[str, ...]:
        """Return normalized caveat codes for the resolved transition target."""

        return tuple(
            str(item).strip()
            for item in self.caveat_resolver(target, evidence)
            if str(item).strip()
        )


def initial_quantitative_contract_state(
    *,
    declared_input_scale_kind: IntensityScaleKind | None,
    explicit_quantitative_meaning: QuantitativeMeaning | None,
) -> QuantitativeContractState:
    """Return the initial semantic state used for pre-execution contract folding."""

    scale_kind = (
        IntensityScaleKind.LINEAR
        if declared_input_scale_kind is None
        else _normalize_scale_kind(declared_input_scale_kind)
    )
    meaning = (
        default_quantitative_meaning_for_scale_kind(scale_kind)
        if explicit_quantitative_meaning is None
        else _normalize_meaning(explicit_quantitative_meaning)
    )
    try:
        return QuantitativeContractState(scale_kind=scale_kind, meaning=meaning)
    except InvalidTransformationStateError as exc:
        raise PhosPyInputError(
            "dataset preprocessing quantitative contract rejected invalid initial "
            f"quantitative state: {exc}"
        ) from exc


def preserve_scale_transition(
    scale_kinds: frozenset[IntensityScaleKind] = ALL_INTENSITY_SCALE_KINDS,
    *,
    output_scale_label: str = "preserved",
) -> QuantitativeScaleTransition:
    """Return a scale transition preserving the input scale kind."""

    return QuantitativeScaleTransition(
        kind=QuantitativeScaleTransitionKind.PRESERVE,
        output_scale_by_input={item: item for item in scale_kinds},
        output_scale_label=output_scale_label,
    )


def linear_to_log2_scale_transition() -> QuantitativeScaleTransition:
    """Return the supported linear-to-log2 scale transition."""

    return QuantitativeScaleTransition(
        kind=QuantitativeScaleTransitionKind.LINEAR_TO_LOG2,
        output_scale_by_input={IntensityScaleKind.LINEAR: IntensityScaleKind.LOG2},
        output_scale_label="log2",
    )


def preserve_meaning_transition(
    meanings: frozenset[QuantitativeMeaning] = ALL_QUANTITATIVE_MEANINGS,
) -> QuantitativeMeaningTransition:
    """Return a transition preserving quantitative meaning."""

    return QuantitativeMeaningTransition(
        kind=QuantitativeMeaningTransitionKind.PRESERVE,
        output_meaning_by_input={item: item for item in meanings},
    )


def static_meaning_transition(
    mapping: Mapping[QuantitativeMeaning, QuantitativeMeaning],
) -> QuantitativeMeaningTransition:
    """Return a static quantitative-meaning mapping transition."""

    return QuantitativeMeaningTransition(
        kind=QuantitativeMeaningTransitionKind.STATIC_MAPPING,
        output_meaning_by_input=dict(mapping),
    )


def total_protein_evidence_meaning_resolver(
    evidence: QuantitativeTransitionEvidence | None,
) -> QuantitativeMeaning:
    """Resolve total-protein correction output meaning from typed row evidence."""

    if evidence is None:
        return QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO
    uncorrected = evidence.total_protein_uncorrected_row_count
    if uncorrected is not None and uncorrected > 0:
        return QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE
    return QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO


def total_protein_evidence_caveats(
    target: QuantitativeMeaning,
    _evidence: QuantitativeTransitionEvidence | None,
) -> tuple[str, ...]:
    """Return caveats for total-protein semantic outputs."""

    if (
        target
        is QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE
    ):
        return ("quantitative_meaning_mixed_total_protein_correction",)
    return ()


def evidence_resolved_meaning_transition(
    *,
    accepted_input_meanings: frozenset[QuantitativeMeaning],
    resolver: _EvidenceMeaningResolver,
    default_output_meaning: QuantitativeMeaning,
) -> QuantitativeMeaningTransition:
    """Return an evidence-resolved output meaning transition."""

    return QuantitativeMeaningTransition(
        kind=QuantitativeMeaningTransitionKind.EVIDENCE_RESOLVED,
        output_meaning_by_input={
            item: default_output_meaning for item in accepted_input_meanings
        },
        evidence_resolver=resolver,
        default_output_meaning=default_output_meaning,
    )


def preserve_quantitative_contract(
    *,
    information_loss: QuantitativeInformationLossKind = (
        QuantitativeInformationLossKind.NONE
    ),
    preserves_abundance: bool = True,
    required_evidence: frozenset[QuantitativeEvidenceRequirement] = frozenset(
        {QuantitativeEvidenceRequirement.NONE}
    ),
    negative_domain_policy: NegativeDomainPolicy = (
        NegativeDomainPolicy.PRESERVES_INPUT_DOMAIN
    ),
    reversibility: QuantitativeReversibilityKind = (
        QuantitativeReversibilityKind.REVERSIBLE
    ),
) -> QuantitativeOperationContract:
    """Return a contract that preserves quantitative scale and meaning."""

    return QuantitativeOperationContract(
        accepted_input_scale_kinds=ALL_INTENSITY_SCALE_KINDS,
        accepted_quantitative_meanings=ALL_QUANTITATIVE_MEANINGS,
        output_scale_transition=preserve_scale_transition(),
        output_meaning_transition=preserve_meaning_transition(),
        preserves_abundance=preserves_abundance,
        negative_domain_policy=negative_domain_policy,
        required_evidence=required_evidence,
        reversibility=reversibility,
        information_loss=information_loss,
    )


def _optional_non_negative_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        raise InvalidTransformationStateError(f"{field_name} must be an int or None")
    if value < 0:
        raise InvalidTransformationStateError(f"{field_name} must be >= 0")
    return int(value)


def _normalize_evidence_metadata(value: object) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise InvalidTransformationStateError(
            "quantitative_transition_evidence.metadata must be a mapping"
        )
    parsed: dict[str, object] = {}
    for key, item in cast(Mapping[object, object], value).items():
        if not isinstance(key, str) or not key.strip():
            raise InvalidTransformationStateError(
                "quantitative_transition_evidence.metadata keys must be "
                "non-empty strings"
            )
        if item is not None and not isinstance(item, str | int | float | bool):
            raise InvalidTransformationStateError(
                "quantitative_transition_evidence.metadata values must be JSON "
                "scalar values"
            )
        parsed[key.strip()] = item
    return parsed


def _normalize_scale_kind(value: IntensityScaleKind | str) -> IntensityScaleKind:
    if isinstance(value, IntensityScaleKind):
        return value
    try:
        return IntensityScaleKind(str(value))
    except ValueError as exc:
        supported = ", ".join(item.value for item in IntensityScaleKind)
        raise InvalidTransformationStateError(
            f"unsupported quantitative scale kind {value!r}; supported: {supported}"
        ) from exc


def _normalize_meaning(value: QuantitativeMeaning | str) -> QuantitativeMeaning:
    if isinstance(value, QuantitativeMeaning):
        return value
    try:
        return QuantitativeMeaning(str(value))
    except ValueError as exc:
        supported = ", ".join(item.value for item in QuantitativeMeaning)
        raise InvalidTransformationStateError(
            f"unsupported quantitative meaning {value!r}; supported: {supported}"
        ) from exc


def _normalize_scale_transition_kind(
    value: QuantitativeScaleTransitionKind | str,
) -> QuantitativeScaleTransitionKind:
    if isinstance(value, QuantitativeScaleTransitionKind):
        return value
    try:
        return QuantitativeScaleTransitionKind(str(value))
    except ValueError as exc:
        supported = ", ".join(item.value for item in QuantitativeScaleTransitionKind)
        raise InvalidTransformationStateError(
            f"unsupported quantitative scale transition kind {value!r}; "
            f"supported: {supported}"
        ) from exc


def _normalize_meaning_transition_kind(
    value: QuantitativeMeaningTransitionKind | str,
) -> QuantitativeMeaningTransitionKind:
    if isinstance(value, QuantitativeMeaningTransitionKind):
        return value
    try:
        return QuantitativeMeaningTransitionKind(str(value))
    except ValueError as exc:
        supported = ", ".join(item.value for item in QuantitativeMeaningTransitionKind)
        raise InvalidTransformationStateError(
            f"unsupported quantitative meaning transition kind {value!r}; "
            f"supported: {supported}"
        ) from exc


def _normalize_negative_domain_policy(
    value: NegativeDomainPolicy | str,
) -> NegativeDomainPolicy:
    if isinstance(value, NegativeDomainPolicy):
        return value
    try:
        return NegativeDomainPolicy(str(value))
    except ValueError as exc:
        supported = ", ".join(item.value for item in NegativeDomainPolicy)
        raise InvalidTransformationStateError(
            f"unsupported negative-domain policy {value!r}; supported: {supported}"
        ) from exc


def _normalize_evidence_requirement(
    value: QuantitativeEvidenceRequirement | str,
) -> QuantitativeEvidenceRequirement:
    if isinstance(value, QuantitativeEvidenceRequirement):
        return value
    try:
        return QuantitativeEvidenceRequirement(str(value))
    except ValueError as exc:
        supported = ", ".join(item.value for item in QuantitativeEvidenceRequirement)
        raise InvalidTransformationStateError(
            f"unsupported quantitative evidence requirement {value!r}; "
            f"supported: {supported}"
        ) from exc


def _normalize_reversibility(
    value: QuantitativeReversibilityKind | str,
) -> QuantitativeReversibilityKind:
    if isinstance(value, QuantitativeReversibilityKind):
        return value
    try:
        return QuantitativeReversibilityKind(str(value))
    except ValueError as exc:
        supported = ", ".join(item.value for item in QuantitativeReversibilityKind)
        raise InvalidTransformationStateError(
            f"unsupported quantitative reversibility {value!r}; supported: {supported}"
        ) from exc


def _normalize_information_loss(
    value: QuantitativeInformationLossKind | str,
) -> QuantitativeInformationLossKind:
    if isinstance(value, QuantitativeInformationLossKind):
        return value
    try:
        return QuantitativeInformationLossKind(str(value))
    except ValueError as exc:
        supported = ", ".join(item.value for item in QuantitativeInformationLossKind)
        raise InvalidTransformationStateError(
            f"unsupported quantitative information-loss kind {value!r}; "
            f"supported: {supported}"
        ) from exc


def _normalize_evidence_mode(
    value: QuantitativeMeaningEvidenceMode | str,
) -> QuantitativeMeaningEvidenceMode:
    if isinstance(value, QuantitativeMeaningEvidenceMode):
        return value
    try:
        return QuantitativeMeaningEvidenceMode(str(value))
    except ValueError as exc:
        supported = ", ".join(item.value for item in QuantitativeMeaningEvidenceMode)
        raise InvalidTransformationStateError(
            f"unsupported quantitative evidence mode {value!r}; supported: {supported}"
        ) from exc


def _normalize_optional_identifier(value: str | None, *, field_name: str) -> str | None:
    if value is None:
        return None
    return _normalize_required_text(value, field_name=field_name)


def _normalize_required_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise InvalidTransformationStateError(f"{field_name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise InvalidTransformationStateError(
            f"{field_name} must be a non-empty string"
        )
    return normalized


def _require_mapping_covers(
    *,
    observed: frozenset[object],
    required: frozenset[object],
    field_name: str,
) -> None:
    missing = required - observed
    if missing:
        values = ", ".join(sorted(str(item) for item in missing))
        raise InvalidTransformationStateError(
            f"quantitative operation contract {field_name} does not cover "
            f"accepted inputs: {values}"
        )


def _validate_meaning_scale_coherence(
    *,
    scale_kind: IntensityScaleKind,
    meaning: QuantitativeMeaning,
) -> None:
    state = IntensityScaleState(
        phospho=_matrix_state_for_scale(scale_kind),
        total=None,
        quantity=meaning,
    )
    if state.quantity is not meaning:  # pragma: no cover - constructor normalizes.
        raise InvalidTransformationStateError(
            "quantitative meaning normalization error"
        )


def _matrix_state_for_scale(
    scale_kind: IntensityScaleKind,
) -> MatrixIntensityScaleState:
    if scale_kind is IntensityScaleKind.LOG2:
        return MatrixIntensityScaleState.log2(
            established_by="phospy.science.transformations.quantitative_contracts"
        )
    return MatrixIntensityScaleState.linear(
        established_by="phospy.science.transformations.quantitative_contracts"
    )


__all__ = [
    "ALL_INTENSITY_SCALE_KINDS",
    "ALL_QUANTITATIVE_MEANINGS",
    "NegativeDomainPolicy",
    "QuantitativeContractState",
    "QuantitativeEvidenceRequirement",
    "QuantitativeInformationLossKind",
    "QuantitativeMeaningTransition",
    "QuantitativeMeaningTransitionKind",
    "QuantitativeOperationContract",
    "QuantitativeReversibilityKind",
    "QuantitativeScaleTransition",
    "QuantitativeScaleTransitionKind",
    "QuantitativeTransitionEvidence",
    "evidence_resolved_meaning_transition",
    "initial_quantitative_contract_state",
    "linear_to_log2_scale_transition",
    "preserve_meaning_transition",
    "preserve_quantitative_contract",
    "preserve_scale_transition",
    "static_meaning_transition",
    "total_protein_evidence_caveats",
    "total_protein_evidence_meaning_resolver",
]

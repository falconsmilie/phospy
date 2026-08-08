"""Activity method identity metadata and computability counters."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from phospy._deprecations import warn_deprecated


@dataclass(frozen=True, slots=True)
class ActivityMethodMetadata:
    """Stable scientific identity metadata for an activity-like scoring method."""

    activity_method_id: str
    activity_method_family: str
    activity_method_label: str
    is_ksea: bool
    is_phosr_kinase_activity_equivalent: bool

    def to_payload(self) -> dict[str, object]:
        """Return a scalar metadata snapshot, not an export/report payload."""

        return {
            "activity_method_id": self.activity_method_id,
            "activity_method_family": self.activity_method_family,
            "activity_method_label": self.activity_method_label,
            "is_ksea": bool(self.is_ksea),
            "is_phosr_kinase_activity_equivalent": bool(
                self.is_phosr_kinase_activity_equivalent
            ),
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityMethodMetadata:
        method_id = str(payload.get("activity_method_id", "")).strip()
        method_family = str(payload.get("activity_method_family", "")).strip()
        method_label = str(payload.get("activity_method_label", "")).strip()
        is_ksea = payload.get("is_ksea")
        is_phosr_equivalent = payload.get("is_phosr_kinase_activity_equivalent")
        if not method_id:
            raise ValueError("activity_method_id must be a non-empty string")
        if not method_family:
            raise ValueError("activity_method_family must be a non-empty string")
        if not method_label:
            raise ValueError("activity_method_label must be a non-empty string")
        if not isinstance(is_ksea, bool):
            raise ValueError("is_ksea must be a bool")
        if not isinstance(is_phosr_equivalent, bool):
            raise ValueError("is_phosr_kinase_activity_equivalent must be a bool")
        return cls(
            activity_method_id=method_id,
            activity_method_family=method_family,
            activity_method_label=method_label,
            is_ksea=is_ksea,
            is_phosr_kinase_activity_equivalent=is_phosr_equivalent,
        )


SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD = ActivityMethodMetadata(
    activity_method_id="simplified_weighted_substrate_activity_v1",
    activity_method_family="heuristic_weighted_substrate_score",
    activity_method_label="simplified weighted substrate activity-like score",
    is_ksea=False,
    is_phosr_kinase_activity_equivalent=False,
)

KSEA_ZSCORE_ACTIVITY_METHOD = ActivityMethodMetadata(
    activity_method_id="ksea_zscore_v1",
    activity_method_family="substrate_set_enrichment",
    activity_method_label="KSEA-style z-score kinase activity score",
    is_ksea=True,
    is_phosr_kinase_activity_equivalent=False,
)

SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD = ActivityMethodMetadata(
    activity_method_id="ssgsea_substrate_enrichment_activity_v1",
    activity_method_family="substrate_set_enrichment",
    activity_method_label="ssGSEA substrate enrichment activity-like score",
    is_ksea=False,
    is_phosr_kinase_activity_equivalent=False,
)


@dataclass(frozen=True, slots=True, init=False)
class ActivityMethodSummary:
    """Method-level score computability counters."""

    kinases_evaluated: int
    kinase_profile_pairs_evaluated: int
    kinase_profile_pairs_computed: int
    kinase_profile_pairs_insufficient_substrates: int
    kinase_profile_pairs_invalid_background_variance: int
    kinase_profile_pairs_no_finite_background_values: int
    kinase_profile_pairs_no_finite_substrate_values: int

    def __init__(
        self,
        *,
        kinases_evaluated: int,
        kinase_profile_pairs_evaluated: int | None = None,
        kinase_profile_pairs_computed: int | None = None,
        kinase_profile_pairs_insufficient_substrates: int | None = None,
        kinase_profile_pairs_invalid_background_variance: int | None = None,
        kinase_profile_pairs_no_finite_background_values: int | None = None,
        kinase_profile_pairs_no_finite_substrate_values: int | None = None,
        kinase_condition_pairs_evaluated: int | None = None,
        kinase_condition_pairs_computed: int | None = None,
        kinase_condition_pairs_insufficient_substrates: int | None = None,
        kinase_condition_pairs_invalid_background_variance: int | None = None,
        kinase_condition_pairs_no_finite_background_values: int | None = None,
        kinase_condition_pairs_no_finite_substrate_values: int | None = None,
    ) -> None:
        object.__setattr__(
            self,
            "kinases_evaluated",
            _coerce_summary_counter_value(
                kinases_evaluated,
                field_name="kinases_evaluated",
            ),
        )
        object.__setattr__(
            self,
            "kinase_profile_pairs_evaluated",
            _resolve_summary_counter_alias(
                field_name="kinase_profile_pairs_evaluated",
                field_value=kinase_profile_pairs_evaluated,
                fallback_field_name="kinase_condition_pairs_evaluated",
                fallback_value=kinase_condition_pairs_evaluated,
            ),
        )
        object.__setattr__(
            self,
            "kinase_profile_pairs_computed",
            _resolve_summary_counter_alias(
                field_name="kinase_profile_pairs_computed",
                field_value=kinase_profile_pairs_computed,
                fallback_field_name="kinase_condition_pairs_computed",
                fallback_value=kinase_condition_pairs_computed,
            ),
        )
        object.__setattr__(
            self,
            "kinase_profile_pairs_insufficient_substrates",
            _resolve_summary_counter_alias(
                field_name="kinase_profile_pairs_insufficient_substrates",
                field_value=kinase_profile_pairs_insufficient_substrates,
                fallback_field_name="kinase_condition_pairs_insufficient_substrates",
                fallback_value=kinase_condition_pairs_insufficient_substrates,
            ),
        )
        object.__setattr__(
            self,
            "kinase_profile_pairs_invalid_background_variance",
            _resolve_summary_counter_alias(
                field_name="kinase_profile_pairs_invalid_background_variance",
                field_value=kinase_profile_pairs_invalid_background_variance,
                fallback_field_name=(
                    "kinase_condition_pairs_invalid_background_variance"
                ),
                fallback_value=kinase_condition_pairs_invalid_background_variance,
            ),
        )
        object.__setattr__(
            self,
            "kinase_profile_pairs_no_finite_background_values",
            _resolve_summary_counter_alias(
                field_name="kinase_profile_pairs_no_finite_background_values",
                field_value=kinase_profile_pairs_no_finite_background_values,
                fallback_field_name=(
                    "kinase_condition_pairs_no_finite_background_values"
                ),
                fallback_value=kinase_condition_pairs_no_finite_background_values,
            ),
        )
        object.__setattr__(
            self,
            "kinase_profile_pairs_no_finite_substrate_values",
            _resolve_summary_counter_alias(
                field_name="kinase_profile_pairs_no_finite_substrate_values",
                field_value=kinase_profile_pairs_no_finite_substrate_values,
                fallback_field_name=(
                    "kinase_condition_pairs_no_finite_substrate_values"
                ),
                fallback_value=kinase_condition_pairs_no_finite_substrate_values,
            ),
        )

    def to_payload(self) -> dict[str, int]:
        """Return computability counters as a plain defensive snapshot."""

        return {
            "kinases_evaluated": int(self.kinases_evaluated),
            "kinase_profile_pairs_evaluated": int(self.kinase_profile_pairs_evaluated),
            "kinase_profile_pairs_computed": int(self.kinase_profile_pairs_computed),
            "kinase_profile_pairs_insufficient_substrates": int(
                self.kinase_profile_pairs_insufficient_substrates
            ),
            "kinase_profile_pairs_invalid_background_variance": int(
                self.kinase_profile_pairs_invalid_background_variance
            ),
            "kinase_profile_pairs_no_finite_background_values": int(
                self.kinase_profile_pairs_no_finite_background_values
            ),
            "kinase_profile_pairs_no_finite_substrate_values": int(
                self.kinase_profile_pairs_no_finite_substrate_values
            ),
        }

    @property
    def kinase_condition_pairs_evaluated(self) -> int:
        _warn_legacy_activity_summary_alias(
            legacy_field_name="kinase_condition_pairs_evaluated",
            field_name="kinase_profile_pairs_evaluated",
        )
        return self.kinase_profile_pairs_evaluated

    @property
    def kinase_condition_pairs_computed(self) -> int:
        _warn_legacy_activity_summary_alias(
            legacy_field_name="kinase_condition_pairs_computed",
            field_name="kinase_profile_pairs_computed",
        )
        return self.kinase_profile_pairs_computed

    @property
    def kinase_condition_pairs_insufficient_substrates(self) -> int:
        _warn_legacy_activity_summary_alias(
            legacy_field_name="kinase_condition_pairs_insufficient_substrates",
            field_name="kinase_profile_pairs_insufficient_substrates",
        )
        return self.kinase_profile_pairs_insufficient_substrates

    @property
    def kinase_condition_pairs_invalid_background_variance(self) -> int:
        _warn_legacy_activity_summary_alias(
            legacy_field_name="kinase_condition_pairs_invalid_background_variance",
            field_name="kinase_profile_pairs_invalid_background_variance",
        )
        return self.kinase_profile_pairs_invalid_background_variance

    @property
    def kinase_condition_pairs_no_finite_background_values(self) -> int:
        _warn_legacy_activity_summary_alias(
            legacy_field_name="kinase_condition_pairs_no_finite_background_values",
            field_name="kinase_profile_pairs_no_finite_background_values",
        )
        return self.kinase_profile_pairs_no_finite_background_values

    @property
    def kinase_condition_pairs_no_finite_substrate_values(self) -> int:
        _warn_legacy_activity_summary_alias(
            legacy_field_name="kinase_condition_pairs_no_finite_substrate_values",
            field_name="kinase_profile_pairs_no_finite_substrate_values",
        )
        return self.kinase_profile_pairs_no_finite_substrate_values

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityMethodSummary:
        return cls(
            kinases_evaluated=_coerce_payload_int(
                payload=payload,
                field_name="kinases_evaluated",
            ),
            kinase_profile_pairs_evaluated=_coerce_payload_int_with_fallback(
                payload=payload,
                field_name="kinase_profile_pairs_evaluated",
                fallback_field_name="kinase_condition_pairs_evaluated",
            ),
            kinase_profile_pairs_computed=_coerce_payload_int_with_fallback(
                payload=payload,
                field_name="kinase_profile_pairs_computed",
                fallback_field_name="kinase_condition_pairs_computed",
            ),
            kinase_profile_pairs_insufficient_substrates=(
                _coerce_payload_int_with_fallback(
                    payload=payload,
                    field_name="kinase_profile_pairs_insufficient_substrates",
                    fallback_field_name=(
                        "kinase_condition_pairs_insufficient_substrates"
                    ),
                )
            ),
            kinase_profile_pairs_invalid_background_variance=(
                _coerce_payload_int_with_fallback(
                    payload=payload,
                    field_name="kinase_profile_pairs_invalid_background_variance",
                    fallback_field_name=(
                        "kinase_condition_pairs_invalid_background_variance"
                    ),
                )
            ),
            kinase_profile_pairs_no_finite_background_values=(
                _coerce_payload_int_with_fallback(
                    payload=payload,
                    field_name="kinase_profile_pairs_no_finite_background_values",
                    fallback_field_name=(
                        "kinase_condition_pairs_no_finite_background_values"
                    ),
                )
            ),
            kinase_profile_pairs_no_finite_substrate_values=(
                _coerce_payload_int_with_fallback(
                    payload=payload,
                    field_name="kinase_profile_pairs_no_finite_substrate_values",
                    fallback_field_name=(
                        "kinase_condition_pairs_no_finite_substrate_values"
                    ),
                )
            ),
        )


def _coerce_payload_int_with_fallback(
    *,
    payload: Mapping[str, object],
    field_name: str,
    fallback_field_name: str,
) -> int:
    if field_name in payload:
        value = _coerce_payload_int(payload=payload, field_name=field_name)
        if fallback_field_name in payload:
            fallback_value = _coerce_payload_int(
                payload=payload,
                field_name=fallback_field_name,
            )
            if value != fallback_value:
                raise ValueError(
                    f"{field_name} conflicts with legacy alias {fallback_field_name}"
                )
        return value
    return _coerce_payload_int(payload=payload, field_name=fallback_field_name)


def _resolve_summary_counter_alias(
    *,
    field_name: str,
    field_value: int | None,
    fallback_field_name: str,
    fallback_value: int | None,
) -> int:
    if field_value is not None:
        value = _coerce_summary_counter_value(field_value, field_name=field_name)
        if fallback_value is not None:
            legacy_value = _coerce_summary_counter_value(
                fallback_value,
                field_name=fallback_field_name,
            )
            if value != legacy_value:
                raise ValueError(
                    f"{field_name} conflicts with legacy alias {fallback_field_name}"
                )
        return value
    if fallback_value is not None:
        _warn_legacy_activity_summary_alias(
            legacy_field_name=fallback_field_name,
            field_name=field_name,
        )
        return _coerce_summary_counter_value(
            fallback_value,
            field_name=fallback_field_name,
        )
    raise ValueError(f"{field_name} must be provided")


def _coerce_summary_counter_value(value: object, *, field_name: str) -> int:
    if not isinstance(value, bool | int | float | str | bytes | bytearray):
        raise ValueError(f"{field_name} must be int-compatible")
    return int(value)


def _warn_legacy_activity_summary_alias(
    *,
    legacy_field_name: str,
    field_name: str,
) -> None:
    warn_deprecated(
        f"activities.method_summary.{legacy_field_name}",
        stacklevel=3,
    )


def _coerce_payload_int(*, payload: Mapping[str, object], field_name: str) -> int:
    value = payload.get(field_name, 0)
    if isinstance(value, bool | int | float | str | bytes | bytearray):
        return int(value)
    raise ValueError(f"{field_name} must be int-compatible")


__all__ = [
    "ActivityMethodMetadata",
    "ActivityMethodSummary",
    "KSEA_ZSCORE_ACTIVITY_METHOD",
    "SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD",
    "SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD",
]

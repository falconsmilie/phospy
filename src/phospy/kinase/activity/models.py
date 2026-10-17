"""Activity method identity metadata and computability counters."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ActivityMethodMetadata:
    """Stable scientific identity metadata for an activity-like scoring method."""

    activity_method_id: str
    activity_method_family: str
    activity_method_label: str
    is_ksea: bool
    is_phosr_kinase_activity_equivalent: bool

    def to_payload(self) -> dict[str, object]:
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
    def from_payload(cls, payload: Mapping[str, object]) -> ActivityMethodMetadata:
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
        kinase_profile_pairs_evaluated: int,
        kinase_profile_pairs_computed: int,
        kinase_profile_pairs_insufficient_substrates: int,
        kinase_profile_pairs_invalid_background_variance: int,
        kinase_profile_pairs_no_finite_background_values: int,
        kinase_profile_pairs_no_finite_substrate_values: int,
    ) -> None:
        for field_name, field_value in (
            ("kinases_evaluated", kinases_evaluated),
            ("kinase_profile_pairs_evaluated", kinase_profile_pairs_evaluated),
            ("kinase_profile_pairs_computed", kinase_profile_pairs_computed),
            (
                "kinase_profile_pairs_insufficient_substrates",
                kinase_profile_pairs_insufficient_substrates,
            ),
            (
                "kinase_profile_pairs_invalid_background_variance",
                kinase_profile_pairs_invalid_background_variance,
            ),
            (
                "kinase_profile_pairs_no_finite_background_values",
                kinase_profile_pairs_no_finite_background_values,
            ),
            (
                "kinase_profile_pairs_no_finite_substrate_values",
                kinase_profile_pairs_no_finite_substrate_values,
            ),
        ):
            object.__setattr__(
                self,
                field_name,
                _coerce_summary_counter_value(field_value, field_name=field_name),
            )

    def to_payload(self) -> dict[str, int]:
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

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> ActivityMethodSummary:
        legacy_fields = sorted(
            str(name) for name in payload if str(name).startswith("kinase_condition_")
        )
        if legacy_fields:
            raise ValueError(
                "unsupported activity summary fields: " + ", ".join(legacy_fields)
            )
        return cls(
            kinases_evaluated=_coerce_payload_int(payload, "kinases_evaluated"),
            kinase_profile_pairs_evaluated=_coerce_payload_int(
                payload, "kinase_profile_pairs_evaluated"
            ),
            kinase_profile_pairs_computed=_coerce_payload_int(
                payload, "kinase_profile_pairs_computed"
            ),
            kinase_profile_pairs_insufficient_substrates=_coerce_payload_int(
                payload, "kinase_profile_pairs_insufficient_substrates"
            ),
            kinase_profile_pairs_invalid_background_variance=_coerce_payload_int(
                payload, "kinase_profile_pairs_invalid_background_variance"
            ),
            kinase_profile_pairs_no_finite_background_values=_coerce_payload_int(
                payload, "kinase_profile_pairs_no_finite_background_values"
            ),
            kinase_profile_pairs_no_finite_substrate_values=_coerce_payload_int(
                payload, "kinase_profile_pairs_no_finite_substrate_values"
            ),
        )


def _coerce_summary_counter_value(value: object, *, field_name: str) -> int:
    if not isinstance(value, bool | int | float | str | bytes | bytearray):
        raise ValueError(f"{field_name} must be int-compatible")
    return int(value)


def _coerce_payload_int(payload: Mapping[str, object], field_name: str) -> int:
    return _coerce_summary_counter_value(
        payload.get(field_name, 0), field_name=field_name
    )


__all__ = [
    "ActivityMethodMetadata",
    "ActivityMethodSummary",
    "KSEA_ZSCORE_ACTIVITY_METHOD",
    "SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD",
    "SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD",
]

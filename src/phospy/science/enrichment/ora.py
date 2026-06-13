"""Offline over-representation analysis engine."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace
from typing import Literal, cast

from scipy.stats import hypergeom

from phospy.science.enrichment.models import (
    ENRICHMENT_METHOD_OVER_REPRESENTATION,
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    SUPPORTED_MULTIPLE_TESTING_CORRECTIONS,
    EnrichmentCollectionKind,
    EnrichmentIdentifierKind,
    EnrichmentSetCollection,
    MultipleTestingCorrection,
)
from phospy.science.statistics.multiple_testing import (
    run as run_multiple_testing_correction,
)

ORA_STATISTICAL_TEST_HYPERGEOMETRIC = "hypergeometric"
OraStatisticalTest = Literal["hypergeometric"]

ORA_OUTSIDE_BACKGROUND_POLICY_ERROR = "error"
ORA_OUTSIDE_BACKGROUND_POLICY_DROP = "drop"
OraOutsideBackgroundPolicy = Literal["error", "drop"]

SUPPORTED_ORA_STATISTICAL_TESTS: tuple[OraStatisticalTest, ...] = (
    ORA_STATISTICAL_TEST_HYPERGEOMETRIC,
)
SUPPORTED_ORA_OUTSIDE_BACKGROUND_POLICIES: tuple[
    OraOutsideBackgroundPolicy,
    ...,
] = (
    ORA_OUTSIDE_BACKGROUND_POLICY_ERROR,
    ORA_OUTSIDE_BACKGROUND_POLICY_DROP,
)


@dataclass(frozen=True, slots=True)
class OraConfig:
    """Method configuration for offline over-representation analysis."""

    statistical_test: OraStatisticalTest = ORA_STATISTICAL_TEST_HYPERGEOMETRIC
    selected_outside_background_policy: OraOutsideBackgroundPolicy = (
        ORA_OUTSIDE_BACKGROUND_POLICY_ERROR
    )
    set_outside_background_policy: OraOutsideBackgroundPolicy = (
        ORA_OUTSIDE_BACKGROUND_POLICY_DROP
    )
    multiple_testing_correction: MultipleTestingCorrection = (
        MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "statistical_test",
            cast(
                OraStatisticalTest,
                _require_supported_value(
                    self.statistical_test,
                    field_name="ora_config.statistical_test",
                    supported=SUPPORTED_ORA_STATISTICAL_TESTS,
                ),
            ),
        )
        object.__setattr__(
            self,
            "selected_outside_background_policy",
            cast(
                OraOutsideBackgroundPolicy,
                _require_supported_value(
                    self.selected_outside_background_policy,
                    field_name="ora_config.selected_outside_background_policy",
                    supported=SUPPORTED_ORA_OUTSIDE_BACKGROUND_POLICIES,
                ),
            ),
        )
        object.__setattr__(
            self,
            "set_outside_background_policy",
            cast(
                OraOutsideBackgroundPolicy,
                _require_supported_value(
                    self.set_outside_background_policy,
                    field_name="ora_config.set_outside_background_policy",
                    supported=SUPPORTED_ORA_OUTSIDE_BACKGROUND_POLICIES,
                ),
            ),
        )
        object.__setattr__(
            self,
            "multiple_testing_correction",
            cast(
                MultipleTestingCorrection,
                _require_supported_value(
                    self.multiple_testing_correction,
                    field_name="ora_config.multiple_testing_correction",
                    supported=SUPPORTED_MULTIPLE_TESTING_CORRECTIONS,
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class OraResultRecord:
    """One ORA row for one enrichment set."""

    set_id: str
    name: str
    collection_kind: EnrichmentCollectionKind
    identifier_kind: EnrichmentIdentifierKind
    background_size: int
    selected_size: int
    raw_set_size: int
    set_size: int
    overlap_size: int
    overlap_identifiers: tuple[str, ...]
    p_value: float
    enrichment_ratio: float | None
    set_identifiers_outside_background_count: int
    adjusted_p_value: float | None = None
    correction_method: MultipleTestingCorrection | None = None


@dataclass(frozen=True, slots=True)
class OraResult:
    """Offline ORA result with deterministic rows and background diagnostics."""

    method: str
    config: OraConfig
    background_size: int
    selected_size: int
    selected_identifiers: tuple[str, ...]
    dropped_selected_identifiers: tuple[str, ...]
    records: tuple[OraResultRecord, ...]


@dataclass(frozen=True, slots=True)
class OraEngine:
    """Pure offline over-representation analysis engine."""

    def run(
        self,
        *,
        selected_identifiers: Sequence[str],
        background_universe: Sequence[str],
        enrichment_sets: EnrichmentSetCollection,
        config: OraConfig | None = None,
    ) -> OraResult:
        return run(
            selected_identifiers=selected_identifiers,
            background_universe=background_universe,
            enrichment_sets=enrichment_sets,
            config=config,
        )


def run(
    *,
    selected_identifiers: Sequence[str],
    background_universe: Sequence[str],
    enrichment_sets: EnrichmentSetCollection,
    config: OraConfig | None = None,
) -> OraResult:
    """Run offline ORA using a hypergeometric survival-function p-value."""

    resolved_config = OraConfig() if config is None else _require_ora_config(config)
    if resolved_config.statistical_test != ORA_STATISTICAL_TEST_HYPERGEOMETRIC:
        raise ValueError("ora statistical_test must be 'hypergeometric'")
    if not isinstance(enrichment_sets, EnrichmentSetCollection):
        raise ValueError("enrichment_sets must be an EnrichmentSetCollection")

    background = _normalise_identifiers(
        background_universe,
        field_name="background_universe",
        allow_empty=False,
    )
    background_set = frozenset(background)
    selected_raw = _normalise_identifiers(
        selected_identifiers,
        field_name="selected_identifiers",
        allow_empty=True,
    )
    selected, dropped_selected = _apply_selected_background_policy(
        selected=selected_raw,
        background=background_set,
        policy=resolved_config.selected_outside_background_policy,
    )
    selected_set = frozenset(selected)

    records: list[OraResultRecord] = []
    for enrichment_set in enrichment_sets.enrichment_sets:
        set_members_raw = tuple(enrichment_set.identifiers)
        set_members_outside_background = tuple(
            identifier
            for identifier in set_members_raw
            if identifier not in background_set
        )
        if (
            set_members_outside_background
            and resolved_config.set_outside_background_policy
            == ORA_OUTSIDE_BACKGROUND_POLICY_ERROR
        ):
            formatted = ", ".join(repr(item) for item in set_members_outside_background)
            raise ValueError(
                "enrichment set identifiers must be within background_universe; "
                f"set_id={enrichment_set.set_id!r}, outside={formatted}"
            )

        set_members = tuple(
            identifier for identifier in set_members_raw if identifier in background_set
        )
        set_member_set = frozenset(set_members)
        overlap_identifiers = tuple(sorted(selected_set.intersection(set_member_set)))
        p_value = _hypergeometric_survival_p_value(
            background_size=len(background_set),
            set_size=len(set_member_set),
            selected_size=len(selected_set),
            overlap_size=len(overlap_identifiers),
        )
        enrichment_ratio = _enrichment_ratio(
            background_size=len(background_set),
            set_size=len(set_member_set),
            selected_size=len(selected_set),
            overlap_size=len(overlap_identifiers),
        )
        records.append(
            OraResultRecord(
                set_id=enrichment_set.set_id,
                name=enrichment_set.name,
                collection_kind=enrichment_sets.collection_kind,
                identifier_kind=enrichment_sets.identifier_kind,
                background_size=len(background_set),
                selected_size=len(selected_set),
                raw_set_size=len(set_members_raw),
                set_size=len(set_member_set),
                overlap_size=len(overlap_identifiers),
                overlap_identifiers=overlap_identifiers,
                p_value=p_value,
                enrichment_ratio=enrichment_ratio,
                set_identifiers_outside_background_count=len(
                    set_members_outside_background
                ),
            )
        )

    corrected_records = _apply_multiple_testing_correction(
        records=tuple(records),
        method=resolved_config.multiple_testing_correction,
    )

    return OraResult(
        method=ENRICHMENT_METHOD_OVER_REPRESENTATION,
        config=resolved_config,
        background_size=len(background_set),
        selected_size=len(selected_set),
        selected_identifiers=selected,
        dropped_selected_identifiers=dropped_selected,
        records=tuple(
            sorted(
                corrected_records, key=lambda record: (record.p_value, record.set_id)
            )
        ),
    )


def _apply_multiple_testing_correction(
    *,
    records: tuple[OraResultRecord, ...],
    method: MultipleTestingCorrection,
) -> tuple[OraResultRecord, ...]:
    if not records:
        return ()
    adjusted_p_values = run_multiple_testing_correction(
        tuple(record.p_value for record in records),
        method=method,
    )
    return tuple(
        replace(
            record,
            adjusted_p_value=adjusted_p_value,
            correction_method=method,
        )
        for record, adjusted_p_value in zip(
            records,
            adjusted_p_values,
            strict=True,
        )
    )


def _normalise_identifiers(
    value: object,
    *,
    field_name: str,
    allow_empty: bool,
) -> tuple[str, ...]:
    if isinstance(value, str | bytes | bytearray) or not isinstance(value, Sequence):
        raise ValueError(f"{field_name} must be a sequence of strings")
    seen: set[str] = set()
    normalised: list[str] = []
    for raw_item in value:
        if not isinstance(raw_item, str):
            raise ValueError(f"{field_name} values must be strings")
        item = raw_item.strip()
        if not item:
            raise ValueError(f"{field_name} values must be non-empty strings")
        if item in seen:
            continue
        seen.add(item)
        normalised.append(item)
    if not allow_empty and not normalised:
        raise ValueError(f"{field_name} must not be empty")
    return tuple(normalised)


def _apply_selected_background_policy(
    *,
    selected: tuple[str, ...],
    background: frozenset[str],
    policy: OraOutsideBackgroundPolicy,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    outside_background = tuple(
        identifier for identifier in selected if identifier not in background
    )
    if outside_background and policy == ORA_OUTSIDE_BACKGROUND_POLICY_ERROR:
        formatted = ", ".join(repr(item) for item in outside_background)
        raise ValueError(
            "selected_identifiers must be within background_universe; "
            f"outside={formatted}"
        )
    if policy == ORA_OUTSIDE_BACKGROUND_POLICY_DROP:
        return (
            tuple(identifier for identifier in selected if identifier in background),
            tuple(sorted(outside_background)),
        )
    return selected, ()


def _hypergeometric_survival_p_value(
    *,
    background_size: int,
    set_size: int,
    selected_size: int,
    overlap_size: int,
) -> float:
    if overlap_size <= 0 or selected_size <= 0 or set_size <= 0:
        return 1.0
    p_value = float(
        hypergeom.sf(
            overlap_size - 1,
            background_size,
            set_size,
            selected_size,
        )
    )
    return min(max(p_value, 0.0), 1.0)


def _enrichment_ratio(
    *,
    background_size: int,
    set_size: int,
    selected_size: int,
    overlap_size: int,
) -> float | None:
    if background_size <= 0 or set_size <= 0 or selected_size <= 0:
        return None
    expected_overlap = float(selected_size) * float(set_size) / float(background_size)
    if expected_overlap <= 0.0:
        return None
    return float(overlap_size) / expected_overlap


def _require_ora_config(value: object) -> OraConfig:
    if not isinstance(value, OraConfig):
        raise ValueError("config must be OraConfig or None")
    return value


def _require_supported_value(
    value: object,
    *,
    field_name: str,
    supported: tuple[str, ...],
) -> str:
    if not isinstance(value, str) or value not in supported:
        allowed = ", ".join(repr(item) for item in supported)
        raise ValueError(f"{field_name} must be one of: {allowed}")
    return value


__all__ = [
    "ORA_OUTSIDE_BACKGROUND_POLICY_DROP",
    "ORA_OUTSIDE_BACKGROUND_POLICY_ERROR",
    "ORA_STATISTICAL_TEST_HYPERGEOMETRIC",
    "OraConfig",
    "OraEngine",
    "OraOutsideBackgroundPolicy",
    "OraResult",
    "OraResultRecord",
    "OraStatisticalTest",
    "SUPPORTED_ORA_OUTSIDE_BACKGROUND_POLICIES",
    "SUPPORTED_ORA_STATISTICAL_TESTS",
    "run",
]

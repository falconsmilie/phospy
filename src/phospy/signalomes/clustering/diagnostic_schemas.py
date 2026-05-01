"""Typed schemas and validators for stable clustering diagnostics payloads."""

from __future__ import annotations

from collections.abc import Mapping
from math import isfinite
from numbers import Integral, Real
from typing import Literal, TypedDict, cast

from phospy.signalomes.clustering.policies import (
    SIGNALOME_CANDIDATE_SCORING_APPLIES_TO,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY,
)

SIGNALOME_TREE_GENERATION_MODE_FULL_EXACT = "full_exact_tree_construction"
SIGNALOME_TREE_GENERATION_SCOPE_MODULE_SELECTION_AND_ASSIGNMENT = (
    "module_count_selection_and_final_assignment"
)
SignalomeTreeGenerationMode = Literal["full_exact_tree_construction"]
SignalomeTreeGenerationScope = Literal["module_count_selection_and_final_assignment"]
SignalomeCandidateScoringScope = Literal["candidate_module_count_evaluation_only"]
SignalomeCandidateScoringSamplingMethod = Literal[
    "deterministic_uniform_without_replacement"
]
SignalomeCandidateScoringSamplingSeedPolicy = Literal[
    "order_invariant_seed_from_row_hashes_and_sample_size"
]


class SignalomeCandidateScoringSampleCountSummary(TypedDict):
    """Stable sampled-count summary fields."""

    min: int
    max: int
    mean: float
    total: int


class SignalomeCandidateScoringSamplingDiagnostics(TypedDict):
    """Stable sampled candidate-scoring diagnostics payload."""

    sampling_cap: int
    sampling_method: SignalomeCandidateScoringSamplingMethod
    deterministic_seed_policy: SignalomeCandidateScoringSamplingSeedPolicy
    actual_sampled_pair_count: int
    per_cluster_sample_count_summary: SignalomeCandidateScoringSampleCountSummary


class SignalomeTreeEngineDiagnostics(TypedDict):
    """Stable tree-engine diagnostics emitted by one backend implementation."""

    uses_scipy: bool
    linkage_method: str
    distance_metric: str


class SignalomeBackendDiagnostics(TypedDict):
    """Stable clustering-backend diagnostics payload."""

    backend_name: str
    backend_version: str
    tree_implementation: str
    tree_implementation_version: str
    # Legacy keys retained for historical/internal metadata compatibility.
    tree_engine: str
    tree_engine_version: str
    uses_scipy: bool
    linkage_method: str
    distance_metric: str
    selected_module_count: int
    input_site_count: int
    exact_tree_path_used: bool
    tree_generation_mode: SignalomeTreeGenerationMode
    tree_generation_is_approximate: bool
    tree_generation_scope: SignalomeTreeGenerationScope
    candidate_scoring_scope: SignalomeCandidateScoringScope


class SignalomeClusteringThresholdMetadata(TypedDict):
    """Stable threshold-metadata payload."""

    primary_threshold: float
    fallback_threshold: float


class SignalomeClusteringLimitMetadata(TypedDict):
    """Stable limit-metadata payload."""

    max_exact_tree_sites: int | None
    max_full_candidate_scoring_sites: int
    max_clusters: int


def build_tree_engine_diagnostics(
    *,
    uses_scipy: bool,
    linkage_method: str,
    distance_metric: str,
) -> SignalomeTreeEngineDiagnostics:
    """Build and validate a stable tree-engine diagnostics record."""

    payload: SignalomeTreeEngineDiagnostics = {
        "uses_scipy": _require_bool(uses_scipy, field_name="uses_scipy"),
        "linkage_method": _require_non_empty_str(
            linkage_method,
            field_name="linkage_method",
        ),
        "distance_metric": _require_non_empty_str(
            distance_metric,
            field_name="distance_metric",
        ),
    }
    return validate_tree_engine_diagnostics(
        payload,
        field_name="tree_engine_diagnostics",
    )


def build_backend_diagnostics(
    *,
    backend_name: str,
    backend_version: str,
    tree_engine: str,
    tree_engine_version: str,
    tree_engine_diagnostics: SignalomeTreeEngineDiagnostics,
    selected_module_count: int,
    input_site_count: int,
    exact_tree_path_used: bool,
) -> SignalomeBackendDiagnostics:
    """Build and validate backend diagnostics payload."""

    tree_payload = validate_tree_engine_diagnostics(
        tree_engine_diagnostics,
        field_name="tree_engine_diagnostics",
    )
    payload: SignalomeBackendDiagnostics = {
        "backend_name": _require_non_empty_str(
            backend_name,
            field_name="backend_name",
        ),
        "backend_version": _require_non_empty_str(
            backend_version,
            field_name="backend_version",
        ),
        "tree_implementation": _require_non_empty_str(
            tree_engine,
            field_name="tree_implementation",
        ),
        "tree_implementation_version": _require_non_empty_str(
            tree_engine_version,
            field_name="tree_implementation_version",
        ),
        "tree_engine": _require_non_empty_str(
            tree_engine,
            field_name="tree_engine",
        ),
        "tree_engine_version": _require_non_empty_str(
            tree_engine_version,
            field_name="tree_engine_version",
        ),
        "uses_scipy": tree_payload["uses_scipy"],
        "linkage_method": tree_payload["linkage_method"],
        "distance_metric": tree_payload["distance_metric"],
        "selected_module_count": _require_int(
            selected_module_count,
            field_name="selected_module_count",
        ),
        "input_site_count": _require_int(
            input_site_count,
            field_name="input_site_count",
        ),
        "exact_tree_path_used": _require_bool(
            exact_tree_path_used,
            field_name="exact_tree_path_used",
        ),
        "tree_generation_mode": SIGNALOME_TREE_GENERATION_MODE_FULL_EXACT,
        "tree_generation_is_approximate": False,
        "tree_generation_scope": (
            SIGNALOME_TREE_GENERATION_SCOPE_MODULE_SELECTION_AND_ASSIGNMENT
        ),
        "candidate_scoring_scope": SIGNALOME_CANDIDATE_SCORING_APPLIES_TO,
    }
    return validate_backend_diagnostics(
        payload,
        field_name="backend_diagnostics",
    )


def validate_candidate_scoring_sampling_diagnostics(
    payload: Mapping[str, object],
    *,
    field_name: str,
) -> SignalomeCandidateScoringSamplingDiagnostics:
    mapping = _require_mapping(payload, field_name=field_name)
    _require_exact_keys(
        mapping,
        expected={
            "sampling_cap",
            "sampling_method",
            "deterministic_seed_policy",
            "actual_sampled_pair_count",
            "per_cluster_sample_count_summary",
        },
        field_name=field_name,
    )
    sampling_method = _require_non_empty_str(
        mapping["sampling_method"],
        field_name=f"{field_name}.sampling_method",
    )
    if sampling_method != SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD:
        raise ValueError(
            f"{field_name}.sampling_method={sampling_method!r} is unsupported; "
            f"expected {SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD!r}"
        )
    seed_policy = _require_non_empty_str(
        mapping["deterministic_seed_policy"],
        field_name=f"{field_name}.deterministic_seed_policy",
    )
    if seed_policy != SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY:
        raise ValueError(
            f"{field_name}.deterministic_seed_policy={seed_policy!r} is unsupported; "
            f"expected {SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY!r}"
        )
    summary = validate_candidate_scoring_sample_count_summary(
        _require_mapping(
            mapping["per_cluster_sample_count_summary"],
            field_name=f"{field_name}.per_cluster_sample_count_summary",
        ),
        field_name=f"{field_name}.per_cluster_sample_count_summary",
    )
    return {
        "sampling_cap": _require_int(
            mapping["sampling_cap"],
            field_name=f"{field_name}.sampling_cap",
        ),
        "sampling_method": sampling_method,
        "deterministic_seed_policy": seed_policy,
        "actual_sampled_pair_count": _require_int(
            mapping["actual_sampled_pair_count"],
            field_name=f"{field_name}.actual_sampled_pair_count",
        ),
        "per_cluster_sample_count_summary": summary,
    }


def validate_candidate_scoring_sample_count_summary(
    payload: Mapping[str, object],
    *,
    field_name: str,
) -> SignalomeCandidateScoringSampleCountSummary:
    mapping = _require_mapping(payload, field_name=field_name)
    _require_exact_keys(
        mapping,
        expected={"min", "max", "mean", "total"},
        field_name=field_name,
    )
    return {
        "min": _require_int(mapping["min"], field_name=f"{field_name}.min"),
        "max": _require_int(mapping["max"], field_name=f"{field_name}.max"),
        "mean": _require_float(mapping["mean"], field_name=f"{field_name}.mean"),
        "total": _require_int(mapping["total"], field_name=f"{field_name}.total"),
    }


def validate_tree_engine_diagnostics(
    payload: Mapping[str, object],
    *,
    field_name: str,
) -> SignalomeTreeEngineDiagnostics:
    mapping = _require_mapping(payload, field_name=field_name)
    _require_exact_keys(
        mapping,
        expected={"uses_scipy", "linkage_method", "distance_metric"},
        field_name=field_name,
    )
    return {
        "uses_scipy": _require_bool(
            mapping["uses_scipy"],
            field_name=f"{field_name}.uses_scipy",
        ),
        "linkage_method": _require_non_empty_str(
            mapping["linkage_method"],
            field_name=f"{field_name}.linkage_method",
        ),
        "distance_metric": _require_non_empty_str(
            mapping["distance_metric"],
            field_name=f"{field_name}.distance_metric",
        ),
    }


def validate_backend_diagnostics(
    payload: Mapping[str, object],
    *,
    field_name: str,
) -> SignalomeBackendDiagnostics:
    mapping = _require_mapping(payload, field_name=field_name)
    required_keys = {
        "backend_name",
        "uses_scipy",
        "linkage_method",
        "distance_metric",
        "selected_module_count",
        "input_site_count",
        "exact_tree_path_used",
        "tree_generation_mode",
        "tree_generation_is_approximate",
        "tree_generation_scope",
        "candidate_scoring_scope",
    }
    missing = sorted(key for key in required_keys if key not in mapping)
    if missing:
        raise ValueError(f"{field_name} schema mismatch (missing keys: {missing})")

    legacy_tree_engine = _require_non_empty_str(
        mapping.get("tree_engine"),
        field_name=f"{field_name}.tree_engine",
    )
    legacy_tree_engine_version = _require_non_empty_str(
        mapping.get("tree_engine_version"),
        field_name=f"{field_name}.tree_engine_version",
    )
    tree_implementation = _require_non_empty_str(
        mapping.get("tree_implementation", legacy_tree_engine),
        field_name=f"{field_name}.tree_implementation",
    )
    tree_implementation_version = _require_non_empty_str(
        mapping.get("tree_implementation_version", legacy_tree_engine_version),
        field_name=f"{field_name}.tree_implementation_version",
    )
    backend_version = _require_non_empty_str(
        mapping.get("backend_version", "1"),
        field_name=f"{field_name}.backend_version",
    )
    tree_generation_mode = _require_non_empty_str(
        mapping["tree_generation_mode"],
        field_name=f"{field_name}.tree_generation_mode",
    )
    if tree_generation_mode != SIGNALOME_TREE_GENERATION_MODE_FULL_EXACT:
        raise ValueError(
            f"{field_name}.tree_generation_mode={tree_generation_mode!r} is "
            f"unsupported; expected {SIGNALOME_TREE_GENERATION_MODE_FULL_EXACT!r}"
        )
    tree_generation_scope = _require_non_empty_str(
        mapping["tree_generation_scope"],
        field_name=f"{field_name}.tree_generation_scope",
    )
    if (
        tree_generation_scope
        != SIGNALOME_TREE_GENERATION_SCOPE_MODULE_SELECTION_AND_ASSIGNMENT
    ):
        raise ValueError(
            f"{field_name}.tree_generation_scope={tree_generation_scope!r} is "
            "unsupported; expected "
            f"{SIGNALOME_TREE_GENERATION_SCOPE_MODULE_SELECTION_AND_ASSIGNMENT!r}"
        )
    candidate_scoring_scope = _require_non_empty_str(
        mapping["candidate_scoring_scope"],
        field_name=f"{field_name}.candidate_scoring_scope",
    )
    if candidate_scoring_scope != SIGNALOME_CANDIDATE_SCORING_APPLIES_TO:
        raise ValueError(
            f"{field_name}.candidate_scoring_scope={candidate_scoring_scope!r} is "
            f"unsupported; expected {SIGNALOME_CANDIDATE_SCORING_APPLIES_TO!r}"
        )
    return {
        "backend_name": _require_non_empty_str(
            mapping["backend_name"],
            field_name=f"{field_name}.backend_name",
        ),
        "backend_version": backend_version,
        "tree_implementation": tree_implementation,
        "tree_implementation_version": tree_implementation_version,
        "tree_engine": legacy_tree_engine,
        "tree_engine_version": legacy_tree_engine_version,
        "uses_scipy": _require_bool(
            mapping["uses_scipy"],
            field_name=f"{field_name}.uses_scipy",
        ),
        "linkage_method": _require_non_empty_str(
            mapping["linkage_method"],
            field_name=f"{field_name}.linkage_method",
        ),
        "distance_metric": _require_non_empty_str(
            mapping["distance_metric"],
            field_name=f"{field_name}.distance_metric",
        ),
        "selected_module_count": _require_int(
            mapping["selected_module_count"],
            field_name=f"{field_name}.selected_module_count",
        ),
        "input_site_count": _require_int(
            mapping["input_site_count"],
            field_name=f"{field_name}.input_site_count",
        ),
        "exact_tree_path_used": _require_bool(
            mapping["exact_tree_path_used"],
            field_name=f"{field_name}.exact_tree_path_used",
        ),
        "tree_generation_mode": tree_generation_mode,
        "tree_generation_is_approximate": _require_bool(
            mapping["tree_generation_is_approximate"],
            field_name=f"{field_name}.tree_generation_is_approximate",
        ),
        "tree_generation_scope": tree_generation_scope,
        "candidate_scoring_scope": candidate_scoring_scope,
    }


def validate_threshold_metadata(
    payload: Mapping[str, object],
    *,
    field_name: str,
) -> SignalomeClusteringThresholdMetadata:
    mapping = _require_mapping(payload, field_name=field_name)
    _require_exact_keys(
        mapping,
        expected={"primary_threshold", "fallback_threshold"},
        field_name=field_name,
    )
    return {
        "primary_threshold": _require_float(
            mapping["primary_threshold"],
            field_name=f"{field_name}.primary_threshold",
        ),
        "fallback_threshold": _require_float(
            mapping["fallback_threshold"],
            field_name=f"{field_name}.fallback_threshold",
        ),
    }


def validate_limit_metadata(
    payload: Mapping[str, object],
    *,
    field_name: str,
) -> SignalomeClusteringLimitMetadata:
    mapping = _require_mapping(payload, field_name=field_name)
    _require_exact_keys(
        mapping,
        expected={
            "max_exact_tree_sites",
            "max_full_candidate_scoring_sites",
            "max_clusters",
        },
        field_name=field_name,
    )
    max_exact_tree_sites_raw = mapping["max_exact_tree_sites"]
    max_exact_tree_sites = (
        None
        if max_exact_tree_sites_raw is None
        else _require_int(
            max_exact_tree_sites_raw,
            field_name=f"{field_name}.max_exact_tree_sites",
        )
    )
    return {
        "max_exact_tree_sites": max_exact_tree_sites,
        "max_full_candidate_scoring_sites": _require_int(
            mapping["max_full_candidate_scoring_sites"],
            field_name=f"{field_name}.max_full_candidate_scoring_sites",
        ),
        "max_clusters": _require_int(
            mapping["max_clusters"],
            field_name=f"{field_name}.max_clusters",
        ),
    }


def candidate_scoring_sampling_diagnostics_to_payload(
    payload: SignalomeCandidateScoringSamplingDiagnostics,
) -> dict[str, object]:
    normalized = validate_candidate_scoring_sampling_diagnostics(
        payload,
        field_name="candidate_scoring_sampling",
    )
    summary = normalized["per_cluster_sample_count_summary"]
    return {
        "sampling_cap": int(normalized["sampling_cap"]),
        "sampling_method": str(normalized["sampling_method"]),
        "deterministic_seed_policy": str(normalized["deterministic_seed_policy"]),
        "actual_sampled_pair_count": int(normalized["actual_sampled_pair_count"]),
        "per_cluster_sample_count_summary": {
            "min": int(summary["min"]),
            "max": int(summary["max"]),
            "mean": float(summary["mean"]),
            "total": int(summary["total"]),
        },
    }


def backend_diagnostics_to_payload(
    payload: SignalomeBackendDiagnostics,
) -> dict[str, object]:
    normalized = validate_backend_diagnostics(
        payload,
        field_name="backend_diagnostics",
    )
    return {
        "backend_name": str(normalized["backend_name"]),
        "backend_version": str(normalized["backend_version"]),
        "tree_implementation": str(normalized["tree_implementation"]),
        "tree_implementation_version": str(normalized["tree_implementation_version"]),
        "tree_engine": str(normalized["tree_engine"]),
        "tree_engine_version": str(normalized["tree_engine_version"]),
        "uses_scipy": bool(normalized["uses_scipy"]),
        "linkage_method": str(normalized["linkage_method"]),
        "distance_metric": str(normalized["distance_metric"]),
        "selected_module_count": int(normalized["selected_module_count"]),
        "input_site_count": int(normalized["input_site_count"]),
        "exact_tree_path_used": bool(normalized["exact_tree_path_used"]),
        "tree_generation_mode": str(normalized["tree_generation_mode"]),
        "tree_generation_is_approximate": bool(
            normalized["tree_generation_is_approximate"]
        ),
        "tree_generation_scope": str(normalized["tree_generation_scope"]),
        "candidate_scoring_scope": str(normalized["candidate_scoring_scope"]),
    }


def _require_mapping(value: object, *, field_name: str) -> Mapping[str, object]:
    if isinstance(value, Mapping):
        mapping = cast(Mapping[object, object], value)
        return {str(key): item for key, item in mapping.items()}
    raise ValueError(f"{field_name} must be a mapping")


def _require_exact_keys(
    payload: Mapping[str, object],
    *,
    expected: set[str],
    field_name: str,
) -> None:
    actual = {str(key) for key in payload}
    missing = sorted(expected - actual)
    unknown = sorted(actual - expected)
    if missing or unknown:
        details: list[str] = []
        if missing:
            details.append(f"missing keys: {missing}")
        if unknown:
            details.append(f"unexpected keys: {unknown}")
        raise ValueError(f"{field_name} schema mismatch ({'; '.join(details)})")


def _require_non_empty_str(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{field_name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise ValueError(f"{field_name} must be a non-empty string")
    return normalized


def _require_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{field_name} must be an int")
    return int(value)


def _require_float(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{field_name} must be a float")
    resolved = float(value)
    if not isfinite(resolved):
        raise ValueError(f"{field_name} must be finite")
    return resolved


def _require_bool(value: object, *, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{field_name} must be a bool")
    return bool(value)


__all__ = [
    "SIGNALOME_TREE_GENERATION_MODE_FULL_EXACT",
    "SIGNALOME_TREE_GENERATION_SCOPE_MODULE_SELECTION_AND_ASSIGNMENT",
    "SignalomeBackendDiagnostics",
    "SignalomeCandidateScoringSampleCountSummary",
    "SignalomeCandidateScoringSamplingDiagnostics",
    "SignalomeClusteringLimitMetadata",
    "SignalomeClusteringThresholdMetadata",
    "SignalomeTreeEngineDiagnostics",
    "backend_diagnostics_to_payload",
    "build_backend_diagnostics",
    "build_tree_engine_diagnostics",
    "candidate_scoring_sampling_diagnostics_to_payload",
    "validate_backend_diagnostics",
    "validate_candidate_scoring_sample_count_summary",
    "validate_candidate_scoring_sampling_diagnostics",
    "validate_limit_metadata",
    "validate_threshold_metadata",
    "validate_tree_engine_diagnostics",
]

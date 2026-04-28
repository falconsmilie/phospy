"""Signalome clustering request validation helpers."""

from __future__ import annotations

from phospy.errors.workflows import SignalomeModuleCountValidationError


def validate_requested_module_count(
    *,
    requested_module_count: int | None,
    available_clustering_site_count: int,
    field_name: str,
    require_explicit: bool = False,
) -> int | None:
    """Validate module-count requests against available clustering sites.

    This validator enforces strict module-count policy for signalome clustering:
    request values are never clamped, coerced, or silently rewritten.
    """

    available = int(available_clustering_site_count)
    if available < 0:
        raise ValueError("available_clustering_site_count must be >= 0")

    if requested_module_count is None:
        if require_explicit:
            raise SignalomeModuleCountValidationError(
                _module_count_error_message(
                    field_name=field_name,
                    requested_module_count=None,
                    available_clustering_site_count=available,
                    reason=(
                        "module count is required for this workflow and cannot be "
                        "omitted"
                    ),
                )
            )
        return None

    requested = int(requested_module_count)
    if requested < 1:
        raise SignalomeModuleCountValidationError(
            _module_count_error_message(
                field_name=field_name,
                requested_module_count=requested,
                available_clustering_site_count=available,
                reason="module count must be greater than or equal to 1",
            )
        )
    if requested > available:
        raise SignalomeModuleCountValidationError(
            _module_count_error_message(
                field_name=field_name,
                requested_module_count=requested,
                available_clustering_site_count=available,
                reason=(
                    "requested module count exceeds available clustering sites and "
                    "cannot produce non-empty site modules"
                ),
            )
        )
    return requested


def validate_cluster_count_for_site_count(
    *,
    cluster_count: int,
    available_clustering_site_count: int,
    field_name: str,
) -> int:
    """Validate a resolved cluster-count request for strict non-clamping behaviour."""

    available = int(available_clustering_site_count)
    if available < 0:
        raise ValueError("available_clustering_site_count must be >= 0")
    resolved_count = int(cluster_count)
    if resolved_count < 1:
        raise ValueError(f"{field_name} must be greater than or equal to 1")
    if resolved_count > available:
        raise ValueError(
            f"{field_name}={resolved_count} exceeds available clustering "
            f"site count ({available})"
        )
    return resolved_count


def _module_count_error_message(
    *,
    field_name: str,
    requested_module_count: int | None,
    available_clustering_site_count: int,
    reason: str,
) -> str:
    requested_label = (
        "missing"
        if requested_module_count is None
        else str(int(requested_module_count))
    )
    return (
        "invalid signalome module-count request: "
        f"field={field_name}; "
        f"requested_module_count={requested_label}; "
        f"available_clustering_site_count={int(available_clustering_site_count)}; "
        f"{reason}; "
        "choose a module count between 1 and the number of available clustering sites."
    )


__all__ = [
    "validate_cluster_count_for_site_count",
    "validate_requested_module_count",
]

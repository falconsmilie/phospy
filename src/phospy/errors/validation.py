from __future__ import annotations

from collections.abc import Sequence

from pydantic import ValidationError

from .base import PhospyError


class PhospyValidationError(PhospyError, ValueError):
    """Base class for phospy validation failures."""


class RequestValidationError(PhospyValidationError):
    """Raised when a validated request object cannot be created."""

    @classmethod
    def from_pydantic(
        cls,
        *,
        context: str,
        error: ValidationError,
    ) -> RequestValidationError:
        details = _format_pydantic_errors(error)
        return cls(f"{context}: {details}")


class TableSchemaError(PhospyValidationError):
    """Raised when a tabular input fails schema validation."""


class InputCompatibilityError(PhospyValidationError):
    """Raised when otherwise valid inputs are incompatible with each other."""


class NoCandidateKinasesError(InputCompatibilityError):
    """Raised when prediction thresholds leave no kinase candidates to score."""


class PredictionConfigurationError(PhospyValidationError):
    """Raised when predictor configuration cannot produce a valid training run."""


class TraceError(PhospyValidationError):
    """Raised when prediction tracing is misconfigured or cannot be loaded."""


def format_overlap_failure_message(
    *,
    pred_context: str,
    matrix_context: str,
    overlap_count: int,
    pred_mat_rows: int,
    matrix_rows: int,
    min_overlap: int | None = None,
    min_fraction: float | None = None,
    seam: str | None = None,
) -> str:
    """Build a seam-focused overlap diagnostic message."""

    seam_name = seam or f"{pred_context}/{matrix_context} overlap seam"
    matrix_percent = overlap_count / max(matrix_rows, 1)
    pred_percent = overlap_count / max(pred_mat_rows, 1)
    shared_counts = (
        f"shared={overlap_count}, {matrix_context} rows={matrix_rows}, "
        f"{pred_context} rows={pred_mat_rows}"
    )

    if overlap_count == 0:
        detail = (
            f"{pred_context} and {matrix_context} have no overlapping phosphosite IDs "
            f"at the {seam_name} ({shared_counts})."
        )
    else:
        detail = (
            f"{pred_context} and {matrix_context} have insufficient overlapping "
            "phosphosite IDs "
            f"at the {seam_name}: {shared_counts}, "
            f"{matrix_context} coverage={matrix_percent * 100.0:.1f}%, "
            f"{pred_context} coverage={pred_percent * 100.0:.1f}%."
        )

    if min_overlap is not None and min_fraction is not None:
        detail = (
            f"{detail} Required thresholds: min_overlap={min_overlap}, "
            f"min_fraction={min_fraction:.3f} ({min_fraction * 100.0:.1f}%) "
            f"for {matrix_context}."
        )

    return (
        f"{detail} Next step: confirm both inputs use the same phosphosite ID space "
        "(for example 'ENTITY;SITE;') or regenerate inputs from one reference seam."
    )


def format_no_candidate_kinases_message(
    *,
    source_name: str,
    top: int,
    score_threshold: float,
    inclusion: int,
    kinase_count: int | None = None,
    site_count: int | None = None,
    effective_top: int | None = None,
    qualifying_kinases: int | None = None,
    max_qualifying_sites: int | None = None,
    near_miss_kinases: Sequence[tuple[str, int]] = (),
) -> str:
    """Build a strict-threshold candidate shortfall diagnostic message."""

    message = (
        "No candidate kinases qualified for prediction from "
        f"{source_name} using top={top}, score_threshold={score_threshold}, "
        f"and inclusion={inclusion}."
    )

    if kinase_count is not None and site_count is not None:
        message = (
            f"{message} Evaluated {kinase_count} kinase column(s) across "
            f"{site_count} phosphosite row(s)."
        )
    if effective_top is not None:
        message = f"{message} Effective top window per kinase={effective_top}."
    if qualifying_kinases is not None:
        message = (
            f"{message} Kinases with at least one site above score_threshold: "
            f"{qualifying_kinases}."
        )
    if max_qualifying_sites is not None:
        message = (
            f"{message} Best-support kinase had {max_qualifying_sites} qualifying "
            f"site(s), below inclusion={inclusion}."
        )
    if near_miss_kinases:
        preview = ", ".join(
            f"{kinase} ({count})" for kinase, count in near_miss_kinases
        )
        message = (
            f"{message} Near-miss kinases below inclusion: {preview} "
            "(qualifying site counts)."
        )

    return f"{message} Lower score_threshold or inclusion, or increase top."


def format_empty_prediction_matrix_message(
    *,
    context: str,
    phosphosite_rows: int,
    source_hint: str = "prediction output",
) -> str:
    """Build diagnostics for predMat-like objects with zero kinase columns."""

    return (
        f"{context} does not contain any kinase columns because no candidate kinases "
        f"qualified for prediction ({context} columns=0, phosphosite rows="
        f"{phosphosite_rows}, source={source_hint}). "
        "Regenerate predMat with less restrictive top, score_threshold, or inclusion "
        "settings, and verify you loaded the expected prediction artifact."
    )


def _format_pydantic_errors(error: ValidationError) -> str:
    parts: list[str] = []
    for item in error.errors(include_url=False):
        location = _format_error_location(item.get("loc", ()))
        message = str(item.get("msg", "Invalid value"))
        parts.append(f"{location}: {message}" if location else message)
    return "; ".join(parts) if parts else str(error)


def _format_error_location(location: Sequence[object]) -> str:
    return ".".join(str(part) for part in location if part not in {None, ""})


__all__ = [
    "format_empty_prediction_matrix_message",
    "format_no_candidate_kinases_message",
    "format_overlap_failure_message",
    "InputCompatibilityError",
    "NoCandidateKinasesError",
    "PhospyValidationError",
    "PredictionConfigurationError",
    "RequestValidationError",
    "TableSchemaError",
    "TraceError",
]

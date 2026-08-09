"""Public facade for common structured workflow-result caveats."""

from __future__ import annotations

from collections.abc import Iterable

from phospy.errors.validation import ContractValidationError
from phospy.science.result_caveats import (
    ResultCaveat,
    ResultCaveatSeverity,
    coerce_result_caveats,
    result_caveats_from_payloads,
)


def validate_result_caveats(
    caveats: Iterable[ResultCaveat],
    *,
    field_name: str,
    error_type: type[Exception] = ContractValidationError,
) -> tuple[ResultCaveat, ...]:
    """Return a tuple of validated common result caveats."""

    return coerce_result_caveats(
        caveats,
        field_name=field_name,
        error_type=error_type,
    )


__all__ = [
    "ResultCaveat",
    "ResultCaveatSeverity",
    "result_caveats_from_payloads",
    "validate_result_caveats",
]

"""Common structured caveats for workflow results."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Literal, TypeAlias, cast

from phospy.errors.validation import ContractValidationError

ResultCaveatSeverity: TypeAlias = Literal["info", "warning", "error"]

_RESULT_CAVEAT_SEVERITIES = frozenset({"info", "warning", "error"})


class _ImmutableDetailsMapping(Mapping[str, object]):
    """Read-only mapping that deep-copies as a plain dict for dataclass helpers."""

    _data: Mapping[str, object]
    __slots__ = ("_data",)

    def __init__(self, data: Mapping[str, object]) -> None:
        object.__setattr__(self, "_data", MappingProxyType(dict(data)))

    def __setattr__(self, name: str, value: object) -> None:
        raise AttributeError("result caveat details are read-only")

    def __delattr__(self, name: str) -> None:
        raise AttributeError("result caveat details are read-only")

    def __getitem__(self, key: str) -> object:
        return self._data[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)

    def __repr__(self) -> str:
        return repr(dict(self._data))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Mapping):
            return False
        return dict(self._data) == dict(other)

    def __deepcopy__(self, memo: dict[int, object]) -> dict[str, object]:
        return dict(self._data)


@dataclass(frozen=True, slots=True)
class ResultCaveat:
    """Machine-readable scientific caveat attached to a workflow result."""

    code: str
    severity: ResultCaveatSeverity
    message: str
    details: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "code",
            _require_non_empty_text(
                self.code,
                field_name="result_caveat.code",
            ),
        )
        object.__setattr__(
            self,
            "severity",
            _require_result_caveat_severity(self.severity),
        )
        object.__setattr__(
            self,
            "message",
            _require_non_empty_text(
                self.message,
                field_name="result_caveat.message",
            ),
        )
        object.__setattr__(
            self,
            "details",
            _ImmutableDetailsMapping(
                _require_details_mapping(
                    self.details,
                    field_name="result_caveat.details",
                )
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible caveat payload."""

        return {
            "code": self.code,
            "severity": self.severity,
            "message": self.message,
            "details": dict(self.details),
        }


def validate_result_caveats(
    caveats: Iterable[ResultCaveat],
    *,
    field_name: str,
    error_type: type[Exception] = ContractValidationError,
) -> tuple[ResultCaveat, ...]:
    """Return a tuple of validated common result caveats."""

    owned = tuple(caveats)
    for caveat in owned:
        if isinstance(caveat, ResultCaveat):
            continue
        raise error_type(f"{field_name} must contain ResultCaveat values")
    return owned


def result_caveats_from_payloads(
    raw_payloads: object,
) -> tuple[ResultCaveat, ...]:
    """Parse persisted caveat payloads into validated ResultCaveat objects."""

    if not isinstance(raw_payloads, list):
        return ()
    caveats: list[ResultCaveat] = []
    for raw_payload in raw_payloads:
        if not isinstance(raw_payload, Mapping):
            continue
        code = raw_payload.get("code")
        severity = raw_payload.get("severity")
        message = raw_payload.get("message")
        details = raw_payload.get("details", {})
        if not isinstance(code, str) or code.strip() == "":
            continue
        if not isinstance(severity, str) or severity not in _RESULT_CAVEAT_SEVERITIES:
            continue
        if not isinstance(message, str) or message.strip() == "":
            continue
        if not isinstance(details, Mapping):
            continue
        caveats.append(
            ResultCaveat(
                code=code,
                severity=cast(ResultCaveatSeverity, severity),
                message=message,
                details=dict(details),
            )
        )
    return tuple(caveats)


def _require_result_caveat_severity(value: object) -> ResultCaveatSeverity:
    if isinstance(value, str) and value in _RESULT_CAVEAT_SEVERITIES:
        return cast(ResultCaveatSeverity, value)
    allowed = ", ".join(sorted(_RESULT_CAVEAT_SEVERITIES))
    raise ContractValidationError(f"result_caveat.severity must be one of: {allowed}")


def _require_non_empty_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise ContractValidationError(f"{field_name} must be a non-empty string")
    return value.strip()


def _require_details_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, object]:
    if not isinstance(value, Mapping):
        raise ContractValidationError(f"{field_name} must be a mapping")
    details: dict[str, object] = {}
    for key, item in value.items():
        if not isinstance(key, str) or key.strip() == "":
            raise ContractValidationError(
                f"{field_name} keys must be non-empty strings"
            )
        details[key] = item
    return details


__all__ = [
    "ResultCaveat",
    "ResultCaveatSeverity",
    "result_caveats_from_payloads",
    "validate_result_caveats",
]

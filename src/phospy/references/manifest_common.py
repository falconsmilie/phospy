"""Reference manifest scalar coercion helpers."""

from __future__ import annotations

from datetime import date


def _optional_string(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text if text else None


def _required_string(value: object, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"reference manifest {field_name} must be non-empty")
    return value.strip()


def _coerce_date(value: date | str) -> date:
    if isinstance(value, date):
        return value
    if not isinstance(value, str) or not value.strip():
        raise ValueError("reference manifest retrieved_at must be YYYY-MM-DD")
    return date.fromisoformat(value.strip())


__all__ = [
    "_coerce_date",
    "_optional_string",
    "_required_string",
]

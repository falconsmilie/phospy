"""Deterministic phosphosite identity message formatting."""

from __future__ import annotations


def identity_incoherence_message(
    *,
    context_label: str,
    field_name: str,
    row_position: int,
    row_label: object,
    detail: str,
) -> str:
    """Format row-scoped site-key incoherence messages."""

    return (
        f"{context_label} is inconsistent with site_key "
        f"at row {row_position} ({row_label!r}) in {field_name}: {detail}"
    )


__all__ = ["identity_incoherence_message"]

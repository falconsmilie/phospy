"""Stable provenance identifiers for quantitative transformers."""

from __future__ import annotations

from collections.abc import Mapping


def resolve_transformer_identifier(
    transformer: object,
    *,
    built_in_identifiers: Mapping[type[object], str],
) -> str:
    """Return a built-in or explicit ID, otherwise the transformer's class path."""

    transformer_type = type(transformer)
    built_in_identifier = built_in_identifiers.get(transformer_type)
    if built_in_identifier is not None:
        return built_in_identifier

    explicit_identifier = _explicit_custom_identifier(
        transformer,
        built_in_types=frozenset(built_in_identifiers),
    )
    if explicit_identifier is not None:
        return explicit_identifier
    return f"{transformer_type.__module__}.{transformer_type.__qualname__}"


def _explicit_custom_identifier(
    transformer: object,
    *,
    built_in_types: frozenset[type[object]],
) -> str | None:
    instance_attributes = getattr(transformer, "__dict__", None)
    if isinstance(instance_attributes, Mapping):
        identifier = _normalized_identifier(
            instance_attributes.get("stable_identifier")
        )
        if identifier is not None:
            return identifier

    for owner_type in type(transformer).__mro__:
        if owner_type in built_in_types:
            continue
        identifier = _normalized_identifier(
            owner_type.__dict__.get("stable_identifier")
        )
        if identifier is not None:
            return identifier
    return None


def _normalized_identifier(value: object) -> str | None:
    if not isinstance(value, str):
        return None
    normalized = value.strip()
    return normalized or None


__all__ = ["resolve_transformer_identifier"]

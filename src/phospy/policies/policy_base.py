from __future__ import annotations

from enum import Enum
from typing import TypeVar

from phospy.errors.input import PhosPyInputError
from phospy.validation.common.config_values import coerce_policy_enum

_PolicyEnumT = TypeVar("_PolicyEnumT", bound="PolicyEnum")


class PolicyEnum(str, Enum):
    """Base class for stable policy enums with strict parsing helpers."""

    def __str__(self) -> str:
        return self.value

    @classmethod
    def parse(
        cls: type[_PolicyEnumT],
        value: object,
        *,
        field_name: str,
    ) -> _PolicyEnumT:
        return coerce_policy_enum(
            cls,
            value,
            field_name=field_name,
            error_type=PhosPyInputError,
        )


__all__ = ["PolicyEnum"]

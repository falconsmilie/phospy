"""Deprecated science-side dataset build validator route."""

from __future__ import annotations

from phospy.errors.input import PhosPyInputError


class DatasetBuildRequestValidator:
    """Placeholder for the moved validation-owned builder request validator."""

    def run(self, request: object) -> object:
        raise PhosPyInputError(
            "DatasetBuildRequestValidator is owned by "
            "phospy.validation.datasets.builder_request"
        )


__all__ = ["DatasetBuildRequestValidator"]

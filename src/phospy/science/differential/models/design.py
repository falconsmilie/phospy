"""Design, contrast, and computation request models for differential analysis."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

from dataclasses import dataclass, field
from typing import cast

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.science.differential.models.empirical_bayes_config import (
    EmpiricalBayesConfig,
)
from phospy.science.statistics.multiple_testing import (
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    MultipleTestingCorrection,
)
from phospy.validation.common.dataframes import (
    require_dataframe,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)


@dataclass(frozen=True, slots=True)
class DesignMatrix:
    """Validated design matrix with samples on rows and coefficients on columns."""

    frame: pd.DataFrame

    def __post_init__(self) -> None:
        frame = own_dataframe(
            self.frame,
            field_name="differential.design",
            error_type=PhosPyInputError,
        )
        _validate_numeric_matrix(
            frame,
            field_name="differential.design",
        )
        object.__setattr__(self, "frame", frame)

    def to_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self.frame)


@dataclass(frozen=True, slots=True)
class ContrastMatrix:
    """Validated contrast matrix with design coefficients on rows."""

    frame: pd.DataFrame

    def __post_init__(self) -> None:
        frame = own_dataframe(
            self.frame,
            field_name="differential.contrasts",
            error_type=PhosPyInputError,
        )
        _validate_numeric_matrix(
            frame,
            field_name="differential.contrasts",
        )
        object.__setattr__(self, "frame", frame)

    def to_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self.frame)


@dataclass(frozen=True, slots=True)
class DifferentialAnalysisRequest:
    """Request payload for moderated OLS-style differential analysis."""

    matrix: pd.DataFrame
    design: DesignMatrix | pd.DataFrame
    contrasts: ContrastMatrix | pd.DataFrame
    empirical_bayes: EmpiricalBayesConfig = field(default_factory=EmpiricalBayesConfig)
    multiple_testing_method: MultipleTestingCorrection = (
        MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG
    )

    def __post_init__(self) -> None:
        matrix = own_dataframe(
            self.matrix,
            field_name="differential.matrix",
            error_type=PhosPyInputError,
        )
        _validate_numeric_matrix(
            matrix,
            field_name="differential.matrix",
        )
        _validate_no_all_constant_site_rows(
            matrix,
            field_name="differential.matrix",
        )
        design = self.design
        if isinstance(design, pd.DataFrame):
            design = DesignMatrix(design)
        if not isinstance(cast(object, design), DesignMatrix):
            raise PhosPyInputError(
                "differential.design must be a DesignMatrix or pandas DataFrame"
            )
        contrasts = self.contrasts
        if isinstance(contrasts, pd.DataFrame):
            contrasts = ContrastMatrix(contrasts)
        if not isinstance(cast(object, contrasts), ContrastMatrix):
            raise PhosPyInputError(
                "differential.contrasts must be a ContrastMatrix or pandas DataFrame"
            )
        if not isinstance(cast(object, self.empirical_bayes), EmpiricalBayesConfig):
            raise PhosPyInputError(
                "differential.empirical_bayes must be an EmpiricalBayesConfig"
            )
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "design", design)
        object.__setattr__(self, "contrasts", contrasts)


def _validate_numeric_matrix(frame: pd.DataFrame, *, field_name: str) -> None:
    require_dataframe(
        frame,
        field_name=field_name,
        allow_empty=False,
        error_type=PhosPyInputError,
    )
    require_non_empty_dataframe(
        frame,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_unique_index(
        frame,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_unique_columns(
        frame,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_numeric_dataframe(
        frame,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_finite_numeric_dataframe(
        frame,
        field_name=field_name,
        error_type=PhosPyInputError,
        allow_missing=False,
    )


def _validate_no_all_constant_site_rows(
    frame: pd.DataFrame,
    *,
    field_name: str,
) -> None:
    values = frame.to_numpy(dtype=float)
    constant_mask = np.all(values == values[:, [0]], axis=1)
    if not np.any(constant_mask):
        return
    invalid_labels = frame.index[constant_mask].astype(str).tolist()
    preview = ", ".join(invalid_labels[:5])
    suffix = "" if len(invalid_labels) <= 5 else ", ..."
    raise PhosPyInputError(
        f"{field_name} contains all-constant site intensities, which are "
        "unsupported for differential analysis; all_constant_sites="
        f"{preview}{suffix}"
    )


__all__ = [
    "ContrastMatrix",
    "DesignMatrix",
    "DifferentialAnalysisRequest",
]

"""Design, contrast, and computation request models for differential analysis."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import cast

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.comparison import dataframe_equals
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.frames.validation import (
    require_dataframe,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.science.differential.linear_model import (
    DifferentialDesignDecomposition,
    DifferentialDesignDecompositionError,
    decompose_differential_design,
)
from phospy.science.differential.models.empirical_bayes_config import (
    EmpiricalBayesConfig,
)
from phospy.science.statistics.multiple_testing import (
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    MultipleTestingCorrection,
)


@dataclass(frozen=True, slots=True, eq=False)
class DesignMatrix:
    """Validated design matrix with samples on rows and coefficients on columns.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit matrix-content comparison.
    """

    __hash__ = object.__hash__

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

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another design matrix has the same content."""

        if not isinstance(other, DesignMatrix):
            return False
        return dataframe_equals(self.frame, other.frame)


@dataclass(frozen=True, slots=True, eq=False)
class ContrastMatrix:
    """Validated contrast matrix with design coefficients on rows.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit matrix-content comparison.
    """

    __hash__ = object.__hash__

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

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another contrast matrix has the same content."""

        if not isinstance(other, ContrastMatrix):
            return False
        return dataframe_equals(self.frame, other.frame)


@dataclass(frozen=True, slots=True, eq=False)
class DifferentialAnalysisRequest:
    """Request payload for moderated OLS-style differential analysis.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit request-content comparison.
    """

    __hash__ = object.__hash__

    matrix: pd.DataFrame
    design: DesignMatrix | pd.DataFrame
    contrasts: ContrastMatrix | pd.DataFrame
    design_decomposition: DifferentialDesignDecomposition | None = None
    empirical_bayes: EmpiricalBayesConfig = field(default_factory=EmpiricalBayesConfig)
    multiple_testing_method: MultipleTestingCorrection = (
        MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG
    )

    def __post_init__(self) -> None:
        self._validate_and_store(assume_matrix_owned=False)

    @classmethod
    def _from_owned(
        cls,
        *,
        matrix: pd.DataFrame,
        design: DesignMatrix | pd.DataFrame,
        contrasts: ContrastMatrix | pd.DataFrame,
        design_decomposition: DifferentialDesignDecomposition | None = None,
        empirical_bayes: EmpiricalBayesConfig | None = None,
        multiple_testing_method: MultipleTestingCorrection = (
            MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG
        ),
    ) -> DifferentialAnalysisRequest:
        request = object.__new__(cls)
        object.__setattr__(request, "matrix", matrix)
        object.__setattr__(request, "design", design)
        object.__setattr__(request, "contrasts", contrasts)
        object.__setattr__(request, "design_decomposition", design_decomposition)
        object.__setattr__(
            request,
            "empirical_bayes",
            empirical_bayes if empirical_bayes is not None else EmpiricalBayesConfig(),
        )
        object.__setattr__(
            request,
            "multiple_testing_method",
            multiple_testing_method,
        )
        request._validate_and_store(assume_matrix_owned=True)
        return request

    def _validate_and_store(self, *, assume_matrix_owned: bool) -> None:
        matrix = own_dataframe(
            self.matrix,
            field_name="differential.matrix",
            error_type=PhosPyInputError,
            assume_owned=assume_matrix_owned,
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
        design_decomposition = self.design_decomposition
        if design_decomposition is None:
            try:
                design_decomposition = decompose_differential_design(
                    design.frame.to_numpy(dtype=float)
                )
            except DifferentialDesignDecompositionError as error:
                raise PhosPyInputError(
                    "differential.design is not admissible for stable "
                    f"moderated-contrast analysis: {error}"
                ) from error
        if not isinstance(
            cast(object, design_decomposition), DifferentialDesignDecomposition
        ):
            raise PhosPyInputError(
                "differential.design_decomposition must be a "
                "DifferentialDesignDecomposition"
            )
        try:
            design_decomposition.assert_matches_design(
                design.frame.to_numpy(dtype=float),
                field_name="differential.design",
            )
        except DifferentialDesignDecompositionError as error:
            raise PhosPyInputError(
                "differential.design_decomposition must describe "
                f"differential.design exactly: {error}"
            ) from error
        if not isinstance(cast(object, self.empirical_bayes), EmpiricalBayesConfig):
            raise PhosPyInputError(
                "differential.empirical_bayes must be an EmpiricalBayesConfig"
            )
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "design", design)
        object.__setattr__(self, "contrasts", contrasts)
        object.__setattr__(
            self,
            "design_decomposition",
            design_decomposition,
        )

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another request has the same scientific content."""

        if not isinstance(other, DifferentialAnalysisRequest):
            return False
        if not isinstance(self.design, DesignMatrix) or not isinstance(
            other.design,
            DesignMatrix,
        ):
            return False
        if not isinstance(self.contrasts, ContrastMatrix) or not isinstance(
            other.contrasts,
            ContrastMatrix,
        ):
            return False
        return (
            dataframe_equals(self.matrix, other.matrix)
            and self.design.scientifically_equals(other.design)
            and self.contrasts.scientifically_equals(other.contrasts)
            and _design_decomposition_equals(
                self.design_decomposition,
                other.design_decomposition,
            )
            and self.empirical_bayes == other.empirical_bayes
            and self.multiple_testing_method == other.multiple_testing_method
        )


def _design_decomposition_equals(
    left: DifferentialDesignDecomposition | None,
    right: DifferentialDesignDecomposition | None,
) -> bool:
    if left is None or right is None:
        return left is right
    return left.scientifically_equals(right)


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
    values: npt.NDArray[np.float64] = np.asarray(
        frame.to_numpy(dtype=float),
        dtype=np.float64,
    )
    constant_mask: npt.NDArray[np.bool_] = np.asarray(
        np.all(values == values[:, [0]], axis=1),
        dtype=bool,
    )
    if not np.any(constant_mask):
        return
    invalid_positions = np.flatnonzero(constant_mask)
    invalid_labels = [str(frame.index[int(position)]) for position in invalid_positions]
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

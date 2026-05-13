"""Validator for differential-analysis requests."""

from __future__ import annotations

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.differential.models import DifferentialAnalysisRequest


class DifferentialAnalysisRequestValidator:
    """Validate cross-object contracts for differential analysis."""

    def run(self, request: object) -> DifferentialAnalysisRequest:
        if not isinstance(request, DifferentialAnalysisRequest):
            raise PhosPyInputError(
                "differential analysis input must be a DifferentialAnalysisRequest"
            )

        matrix = request.matrix
        design = request.design.frame
        contrasts = request.contrasts.frame

        if matrix.shape[1] != design.shape[0]:
            raise PhosPyInputError(
                "differential.matrix.columns count must match differential.design rows; "
                f"matrix_columns={int(matrix.shape[1])}, design_rows={int(design.shape[0])}"
            )

        matrix_samples = pd.Index(matrix.columns)
        design_samples = pd.Index(design.index)
        if not matrix_samples.equals(design_samples):
            if (
                not matrix_samples.isin(design_samples).all()
                or not design_samples.isin(matrix_samples).all()
            ):
                raise PhosPyInputError(
                    "differential.matrix.columns must match differential.design.index "
                    "exactly as a set of sample labels"
                )

        design_terms = pd.Index(design.columns)
        contrast_terms = pd.Index(contrasts.index)
        if not design_terms.equals(contrast_terms):
            if (
                not design_terms.isin(contrast_terms).all()
                or not contrast_terms.isin(design_terms).all()
            ):
                raise PhosPyInputError(
                    "differential.contrasts.index must match differential.design.columns "
                    "exactly as design-term labels"
                )

        if int(contrasts.shape[1]) < 1:
            raise PhosPyInputError(
                "differential.contrasts must contain at least one contrast column"
            )

        return request

"""Internal validator for differential workflow requests."""

from __future__ import annotations

from typing import cast

import pandas as pd

from phospy.api.requests import DifferentialAnalysisRequest
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.differential.models import ContrastMatrix, DesignMatrix
from phospy.errors.validation import WorkflowValidationError
from phospy.workflows.differential.models import (
    ValidatedDifferentialAnalysisRequest,
)


class DifferentialAnalysisValidator:
    """Validate `DifferentialAnalysisRequest` before interpretation."""

    def run(self, request: object) -> ValidatedDifferentialAnalysisRequest:
        if not isinstance(request, DifferentialAnalysisRequest):
            raise WorkflowValidationError(
                "differential workflow input must be a DifferentialAnalysisRequest"
            )
        if not isinstance(request.dataset, AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                "differential workflow request dataset must be AnalysisReadyPhosphoDataset"
            )
        design = request.design
        if not isinstance(design, DesignMatrix):
            raise WorkflowValidationError(
                "differential workflow request design must resolve to DesignMatrix"
            )
        contrasts = request.contrasts
        if not isinstance(contrasts, ContrastMatrix):
            raise WorkflowValidationError(
                "differential workflow request contrasts must resolve to ContrastMatrix"
            )
        design_terms = pd.Index(design.frame.columns)
        contrast_terms = pd.Index(contrasts.frame.index)
        if not design_terms.equals(contrast_terms):
            if (
                not design_terms.isin(contrast_terms).all()
                or not contrast_terms.isin(design_terms).all()
            ):
                raise WorkflowValidationError(
                    "differential workflow request contrasts.index must match "
                    "differential workflow request design.columns exactly as "
                    "design-term labels"
                )
        if int(contrasts.frame.shape[1]) < 1:
            raise WorkflowValidationError(
                "differential workflow request contrasts must contain at least one "
                "contrast column"
            )
        return ValidatedDifferentialAnalysisRequest(
            dataset=request.dataset,
            design=cast(DesignMatrix, design),
            contrasts=cast(ContrastMatrix, contrasts),
            empirical_bayes=request.empirical_bayes,
            multiple_testing=request.multiple_testing,
        )


__all__ = ["DifferentialAnalysisValidator"]

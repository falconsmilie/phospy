"""Active batch-correction matrix result."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.preparation_records.batch_correction import (
    BatchCorrectionReport as _BatchCorrectionReport,
)
from phospy.provenance.immutability import freeze_json_mapping_with_error_type


@dataclass(frozen=True, slots=True, init=False)
class BatchCorrectionResult:
    """Matrix correction output with public report and engine diagnostics."""

    report: _BatchCorrectionReport
    diagnostics: Mapping[str, object]
    _corrected_matrix: pd.DataFrame = field(init=False, repr=False)

    def __init__(
        self,
        *,
        corrected_matrix: pd.DataFrame,
        report: _BatchCorrectionReport,
        diagnostics: Mapping[str, object],
    ) -> None:
        object.__setattr__(
            self,
            "_corrected_matrix",
            own_dataframe(
                corrected_matrix,
                field_name="batch_correction_result.corrected_matrix",
                error_type=PhosPyInputError,
                assume_owned=False,
            ),
        )
        object.__setattr__(self, "report", report)
        object.__setattr__(
            self,
            "diagnostics",
            freeze_json_mapping_with_error_type(
                diagnostics,
                field_name="batch_correction_result.diagnostics",
                error_type=PhosPyInputError,
            ),
        )

    @property
    def corrected_matrix(self) -> pd.DataFrame:
        """Return a corrected phospho-matrix snapshot."""

        return export_dataframe(self._corrected_matrix)

    @property
    def corrected(self) -> pd.DataFrame:
        """Return the corrected phospho matrix."""

        return export_dataframe(self._corrected_matrix)


__all__ = ["BatchCorrectionResult"]

"""Passive collaborator contracts for core batch-correction integration."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Protocol

import pandas as pd

from phospy.preprocessing.batch_correction.contracts.output import (
    CorrectedPreprocessingOutput,
)
from phospy.preprocessing.control_sites.contracts import ControlSiteSet
from phospy.provenance.models._shared import JsonValue
from phospy.provenance.models.batch_correction import BatchCorrectionProvenance


class BatchDesignMetadataValidatorProtocol(Protocol):
    """Validator required before resolving batch-correction design metadata."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        batch_column: str,
        condition_columns: tuple[str, ...],
        context: str,
    ) -> object: ...


class BatchCorrectionAdequacyValidatorProtocol(Protocol):
    """Validator required before linear residual batch correction."""

    def run(
        self,
        *,
        batch_by_sample: Mapping[str, str],
        condition_by_sample: Mapping[str, str],
        sample_order: tuple[str, ...],
        preserve_condition_effects: bool,
    ) -> object: ...


class SpsRuvStyleBatchCorrectionResult(Protocol):
    """Passive SPS/RUV result shape consumed by capability execution."""

    @property
    def corrected_preprocessing_output(self) -> CorrectedPreprocessingOutput: ...

    @property
    def diagnostics(self) -> Mapping[str, JsonValue]: ...

    @property
    def provenance(self) -> BatchCorrectionProvenance: ...


class SpsRuvStyleBatchCorrectionRunner(Protocol):
    """Capability-facing SPS/RUV workflow adapter contract."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        config: object,
        sample_metadata: pd.DataFrame | None,
        control_site_set: ControlSiteSet,
        missingness_policy: object,
        upstream_observation_mask: pd.DataFrame | None,
        site_metadata: pd.DataFrame,
    ) -> SpsRuvStyleBatchCorrectionResult: ...


__all__ = [
    "BatchCorrectionAdequacyValidatorProtocol",
    "BatchDesignMetadataValidatorProtocol",
    "SpsRuvStyleBatchCorrectionResult",
    "SpsRuvStyleBatchCorrectionRunner",
]

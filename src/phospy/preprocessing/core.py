from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import pandas as pd

from ..datasets.schema import DatasetSchema
from ..errors import InputCompatibilityError
from ..internal.constants import (
    GENE_P_SITE_COLUMN,
    LOCALIZATION_PROB_COLUMN,
    PHOSPHO_GENE_COLUMN,
    TOTAL_GENE_COLUMN,
    ComparisonSpec,
)
from ..internal.defaults import (
    DEFAULT_LOCALIZATION_THRESHOLD,
    DEFAULT_MAX_UNMATCHED_FRACTION,
    DEFAULT_MIN_OBSERVED_VALUES,
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
)
from .services import (
    PhosphoPreprocessor,
    ProteinCorrectionService,
    TotalPreprocessor,
)
from .site_matrix import SiteMatrixBuilder, SiteMatrixPolicy, SiteMatrixResult

"""Core preprocessing seam for validated phosphoproteomics inputs.

The preprocessing domain now revolves around three layers:

1. ``DatasetPreprocessing`` binds a dataset workspace to the core path.
2. ``CoreProcessor`` orchestrates the real preprocessing steps for full or
   phospho-only inputs.
3. The step services and site-matrix builder perform the concrete transforms.

``CoreProcessor.process()`` and ``process_phospho_only()`` are defensive
ownership boundaries. Trusted internal call sites should use
``process_owned()`` and ``process_phospho_only_owned()``.

Most callers should start with ``PhosphoDataset.preprocessing.run()`` or the
analysis-ready builder and only reach for this module when they need explicit
control over the core orchestration.

Ownership rule: public boundaries copy once, then all internal owned paths
reuse mutable frames without taking additional full-frame defensive copies.
"""


@dataclass(frozen=True, slots=True)
class CorePreprocessingConfig:
    localization_threshold: float = DEFAULT_LOCALIZATION_THRESHOLD
    min_observed: int = DEFAULT_MIN_OBSERVED_VALUES
    total_sentinel: float = DEFAULT_TOTAL_SENTINEL
    phospho_sentinel: float = DEFAULT_PHOSPHO_SENTINEL
    max_unmatched_fraction: float = DEFAULT_MAX_UNMATCHED_FRACTION
    site_matrix_policy: SiteMatrixPolicy = field(default_factory=SiteMatrixPolicy)


def resolve_core_preprocessing_config(
    *,
    config: CorePreprocessingConfig | None = None,
    localization_threshold: float = DEFAULT_LOCALIZATION_THRESHOLD,
    min_observed: int = DEFAULT_MIN_OBSERVED_VALUES,
    total_sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
    phospho_sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
    max_unmatched_fraction: float = DEFAULT_MAX_UNMATCHED_FRACTION,
    site_matrix_policy: SiteMatrixPolicy | Mapping[str, object] | None = None,
    context: str,
    config_param_name: str,
) -> CorePreprocessingConfig:
    default_config = CorePreprocessingConfig()
    has_scalar_overrides = any(
        (
            localization_threshold != default_config.localization_threshold,
            min_observed != default_config.min_observed,
            float(total_sentinel) != default_config.total_sentinel,
            float(phospho_sentinel) != default_config.phospho_sentinel,
            max_unmatched_fraction != default_config.max_unmatched_fraction,
            site_matrix_policy is not None,
        )
    )

    if config is not None and not isinstance(config, CorePreprocessingConfig):
        msg = (
            f"{context}: {config_param_name} must be a CorePreprocessingConfig instance"
        )
        raise TypeError(msg)

    if config is not None and has_scalar_overrides:
        msg = (
            f"{context}: pass either {config_param_name} or scalar "
            "preprocessing options, not both."
        )
        raise ValueError(msg)

    if config is not None:
        return config

    resolved_site_matrix_policy = (
        default_config.site_matrix_policy
        if site_matrix_policy is None
        else SiteMatrixPolicy.from_value(site_matrix_policy)
    )

    return CorePreprocessingConfig(
        localization_threshold=localization_threshold,
        min_observed=min_observed,
        total_sentinel=float(total_sentinel),
        phospho_sentinel=float(phospho_sentinel),
        max_unmatched_fraction=max_unmatched_fraction,
        site_matrix_policy=resolved_site_matrix_policy,
    )


@dataclass(slots=True)
class CoreProcessingResult:
    """Core preprocessing tables produced for one dataset run."""

    total_unique: pd.DataFrame
    total_filtered: pd.DataFrame
    phospho_filtered: pd.DataFrame
    phospho_corrected: pd.DataFrame
    site_matrix: SiteMatrixResult

    @classmethod
    def from_phospho_only(
        cls,
        *,
        schema: DatasetSchema,
        phospho_filtered: pd.DataFrame,
        phospho_corrected: pd.DataFrame,
        site_matrix: SiteMatrixResult,
    ) -> CoreProcessingResult:
        """Build a core-style result for phospho-only preprocessing runs."""

        empty_total = pd.DataFrame(columns=[TOTAL_GENE_COLUMN, *schema.total_cols])
        return cls(
            total_unique=empty_total.copy(),
            total_filtered=empty_total,
            phospho_filtered=phospho_filtered,
            phospho_corrected=phospho_corrected,
            site_matrix=site_matrix,
        )


class CoreProcessor:
    """Run the core preprocessing pipeline over validated dataset frames.

    This is the orchestration layer used underneath the dataset-bound
    preprocessing facade and the analysis-ready builder.
    """

    def __init__(
        self,
        *,
        schema: DatasetSchema,
        comparisons: Sequence[ComparisonSpec] | None = None,
        total_preprocessor: TotalPreprocessor | None = None,
        phospho_preprocessor: PhosphoPreprocessor | None = None,
        protein_correction_service: ProteinCorrectionService | None = None,
        site_matrix_builder: SiteMatrixBuilder | None = None,
    ) -> None:
        self.schema = schema
        self.comparisons = tuple(comparisons) if comparisons is not None else None
        self.total_preprocessor = total_preprocessor or TotalPreprocessor(
            schema=self.schema
        )
        self.phospho_preprocessor = phospho_preprocessor or PhosphoPreprocessor(
            schema=self.schema
        )
        self.protein_correction_service = protein_correction_service or (
            ProteinCorrectionService(
                schema=self.schema,
                comparisons=self.comparisons,
            )
        )
        self.site_matrix_builder = site_matrix_builder or SiteMatrixBuilder(
            value_cols=self.schema.corrected_cols
        )
        self._validate_services()
        self._validate_site_matrix_builder()

    def _validate_services(self) -> None:
        if self.total_preprocessor.schema != self.schema:
            msg = "Total preprocessor schema must match CoreProcessor.schema"
            raise InputCompatibilityError(msg)
        if self.phospho_preprocessor.schema != self.schema:
            msg = "Phospho preprocessor schema must match CoreProcessor.schema"
            raise InputCompatibilityError(msg)
        if self.protein_correction_service.schema != self.schema:
            msg = "Protein correction service schema must match CoreProcessor.schema"
            raise InputCompatibilityError(msg)

    def _validate_site_matrix_builder(self) -> None:
        builder_value_cols = tuple(self.site_matrix_builder.value_cols)
        if builder_value_cols != self.schema.corrected_cols:
            msg = "Site matrix builder value columns must match schema.corrected_cols"
            raise InputCompatibilityError(msg)

    def prepare_total(
        self,
        total_df: pd.DataFrame,
        *,
        gene_col: str = TOTAL_GENE_COLUMN,
        sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
        min_observed: int = DEFAULT_MIN_OBSERVED_VALUES,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        return self.total_preprocessor.prepare(
            total_df,
            gene_col=gene_col,
            sentinel=sentinel,
            min_observed=min_observed,
        )

    def prepare_phospho(
        self,
        phospho_df: pd.DataFrame,
        *,
        gene_col: str = PHOSPHO_GENE_COLUMN,
        site_col: str = GENE_P_SITE_COLUMN,
        localization_col: str = LOCALIZATION_PROB_COLUMN,
        localization_threshold: float = DEFAULT_LOCALIZATION_THRESHOLD,
        sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        min_observed: int = DEFAULT_MIN_OBSERVED_VALUES,
    ) -> pd.DataFrame:
        return self.phospho_preprocessor.prepare(
            phospho_df,
            gene_col=gene_col,
            site_col=site_col,
            localization_col=localization_col,
            localization_threshold=localization_threshold,
            sentinel=sentinel,
            min_observed=min_observed,
        )

    def correct_to_protein(
        self,
        phospho_df: pd.DataFrame,
        total_df: pd.DataFrame,
        *,
        phospho_gene_col: str = PHOSPHO_GENE_COLUMN,
        total_gene_col: str = TOTAL_GENE_COLUMN,
        max_unmatched_fraction: float = DEFAULT_MAX_UNMATCHED_FRACTION,
    ) -> pd.DataFrame:
        """Correct phosphosite values against total protein with boundary awareness.

        When inputs come from owned preparation paths, the correction service
        reuses the owned fast path and avoids another full-frame defensive copy.
        """
        return self.protein_correction_service.correct(
            phospho_df,
            total_df,
            phospho_gene_col=phospho_gene_col,
            total_gene_col=total_gene_col,
            max_unmatched_fraction=max_unmatched_fraction,
        )

    def add_pairwise_comparisons(
        self,
        corrected_df: pd.DataFrame,
        *,
        output_prefix: str = "p_",
    ) -> pd.DataFrame:
        return self.protein_correction_service.add_pairwise_comparisons(
            corrected_df,
            output_prefix=output_prefix,
        )

    def _build_site_matrix(
        self,
        phospho_corrected: pd.DataFrame,
        *,
        policy: SiteMatrixPolicy,
    ) -> SiteMatrixResult:
        return self.site_matrix_builder.build_owned(
            phospho_corrected,
            policy=policy,
        )

    def _rename_phospho_to_corrected_columns(
        self,
        phospho_filtered: pd.DataFrame,
    ) -> pd.DataFrame:
        return phospho_filtered.rename(
            columns=dict(
                zip(self.schema.phospho_cols, self.schema.corrected_cols, strict=True)
            )
        )

    def process(
        self,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        *,
        config: CorePreprocessingConfig | None = None,
    ) -> CoreProcessingResult:
        """Run core preprocessing through the external defensive ownership boundary.

        This public method protects callers from unintended mutation by taking
        deep copies once at entry, then delegating to :meth:`process_owned`.
        """
        resolved_config = config or CorePreprocessingConfig()
        return self.process_owned(
            total_df=total_df.copy(deep=True),
            phospho_df=phospho_df.copy(deep=True),
            config=resolved_config,
        )

    def process_owned(
        self,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        *,
        config: CorePreprocessingConfig | None = None,
    ) -> CoreProcessingResult:
        """Run core preprocessing on already-owned mutable tables.

        Call this from trusted internal workflows that already own the input
        frames and do not need another defensive full-frame copy.
        """
        resolved_config = config or CorePreprocessingConfig()
        total_unique, total_filtered = self.total_preprocessor.prepare_owned(
            total_df,
            sentinel=resolved_config.total_sentinel,
            min_observed=resolved_config.min_observed,
        )
        phospho_filtered = self.phospho_preprocessor.prepare_owned(
            phospho_df,
            localization_threshold=resolved_config.localization_threshold,
            sentinel=resolved_config.phospho_sentinel,
            min_observed=resolved_config.min_observed,
        )
        phospho_corrected = self.protein_correction_service.correct_owned(
            phospho_filtered,
            total_filtered,
            max_unmatched_fraction=resolved_config.max_unmatched_fraction,
        )
        phospho_corrected = self.add_pairwise_comparisons(phospho_corrected)
        site_matrix = self._build_site_matrix(
            phospho_corrected,
            policy=resolved_config.site_matrix_policy,
        )
        return CoreProcessingResult(
            total_unique=total_unique,
            total_filtered=total_filtered,
            phospho_filtered=phospho_filtered,
            phospho_corrected=phospho_corrected,
            site_matrix=site_matrix,
        )

    def process_phospho_only(
        self,
        phospho_df: pd.DataFrame,
        *,
        config: CorePreprocessingConfig | None = None,
    ) -> CoreProcessingResult:
        """Run phospho-only preprocessing through the defensive boundary."""
        resolved_config = config or CorePreprocessingConfig()
        return self.process_phospho_only_owned(
            phospho_df=phospho_df.copy(deep=True),
            config=resolved_config,
        )

    def process_phospho_only_owned(
        self,
        phospho_df: pd.DataFrame,
        *,
        config: CorePreprocessingConfig | None = None,
    ) -> CoreProcessingResult:
        """Run phospho-only preprocessing on already-owned mutable tables."""
        resolved_config = config or CorePreprocessingConfig()
        phospho_filtered = self.phospho_preprocessor.prepare_owned(
            phospho_df,
            localization_threshold=resolved_config.localization_threshold,
            sentinel=resolved_config.phospho_sentinel,
            min_observed=resolved_config.min_observed,
        )
        phospho_corrected = self._rename_phospho_to_corrected_columns(phospho_filtered)
        phospho_corrected = self.add_pairwise_comparisons(phospho_corrected)
        site_matrix = self._build_site_matrix(
            phospho_corrected,
            policy=resolved_config.site_matrix_policy,
        )
        return CoreProcessingResult.from_phospho_only(
            schema=self.schema,
            phospho_filtered=phospho_filtered,
            phospho_corrected=phospho_corrected,
            site_matrix=site_matrix,
        )

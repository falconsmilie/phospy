"""Internal convention normalisation for dataset builder inputs."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.science.datasets.builders.normalization_reporter import (
    SAMPLE_LABEL_INDEX_POLICY as _SAMPLE_LABEL_INDEX_POLICY,
)
from phospy.science.datasets.builders.normalization_reporter import (
    SITE_IDENTIFIER_INDEX_POLICY as _REPORTER_SITE_IDENTIFIER_INDEX_POLICY,
)
from phospy.science.datasets.builders.normalization_reporter import (
    DatasetConventionNormalisationReporter,
)
from phospy.science.datasets.builders.normalization_reporter import (
    IndexLabelNormalizationPolicy as _IndexLabelNormalizationPolicy,
)
from phospy.science.datasets.builders.sample_metadata_normalizer import (
    SampleMetadataNormalizer,
)
from phospy.science.datasets.builders.site_metadata_normalizer import (
    SiteMetadataNormalizer,
)
from phospy.science.datasets.builders.total_matrix_normalizer import (
    TotalProteinMatrixNormalizer,
)
from phospy.science.sites.identifiers import (
    SiteIdentifierNormalisationRecord,
    SiteIdentifierNormalisationReport,
)

_DEFAULT_REPORTER = DatasetConventionNormalisationReporter()
_SITE_IDENTIFIER_INDEX_POLICY = _REPORTER_SITE_IDENTIFIER_INDEX_POLICY


@dataclass(frozen=True, slots=True)
class NormalizedDatasetInputs:
    """Normalised tables ready for sequence derivation and execution."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    site_identifier_normalisation: SiteIdentifierNormalisationReport | None = None


class DatasetConventionNormalizer:
    """Apply narrow, documented shaping rules for supported inputs.

    The supported index derivation convention can populate `gene_symbol` and `site`
    only. It never infers `protein_id`.
    """

    def __init__(
        self,
        *,
        site_metadata_normalizer: SiteMetadataNormalizer | None = None,
        sample_metadata_normalizer: SampleMetadataNormalizer | None = None,
        total_matrix_normalizer: TotalProteinMatrixNormalizer | None = None,
        reporter: DatasetConventionNormalisationReporter | None = None,
    ) -> None:
        self._reporter = reporter or DatasetConventionNormalisationReporter()
        self._site_metadata_normalizer = site_metadata_normalizer or (
            SiteMetadataNormalizer()
        )
        self._sample_metadata_normalizer = (
            sample_metadata_normalizer or SampleMetadataNormalizer()
        )
        self._total_matrix_normalizer = (
            total_matrix_normalizer or TotalProteinMatrixNormalizer()
        )

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
    ) -> NormalizedDatasetInputs:
        site_identifier_records: list[SiteIdentifierNormalisationRecord] = []
        normalized_phospho = phospho.copy(deep=True)
        normalized_phospho.index = (
            self._reporter.normalize_supported_site_index_if_present(
                normalized_phospho.index,
                field_name="dataset build request phospho.index",
                site_identifier_records=site_identifier_records,
            )
        )
        normalized_phospho.columns = self._reporter.normalize_index_labels(
            normalized_phospho.columns,
            field_name="dataset build request phospho.columns",
            policy=_SAMPLE_LABEL_INDEX_POLICY,
        )
        normalized_site_metadata = self._site_metadata_normalizer.run(
            site_metadata.copy(deep=True),
            phospho_index=normalized_phospho.index,
            site_identifier_records=site_identifier_records,
        )
        normalized_sample_metadata = self._sample_metadata_normalizer.run(
            None if sample_metadata is None else sample_metadata.copy(deep=True),
            phospho_columns=normalized_phospho.columns,
        )
        normalized_total = self._total_matrix_normalizer.run(
            None if total is None else total.copy(deep=True),
            phospho_columns=normalized_phospho.columns,
        )
        return NormalizedDatasetInputs(
            phospho=normalized_phospho,
            site_metadata=normalized_site_metadata,
            sample_metadata=normalized_sample_metadata,
            total=normalized_total,
            site_identifier_normalisation=self._reporter.build_site_identifier_report(
                site_identifier_records
            ),
        )


def _normalize_index_labels(
    index: pd.Index,
    *,
    field_name: str,
    policy: _IndexLabelNormalizationPolicy,
) -> pd.Index:
    return _DEFAULT_REPORTER.normalize_index_labels(
        index,
        field_name=field_name,
        policy=policy,
    )


def _normalize_supported_site_index_if_present(
    index: pd.Index,
    *,
    field_name: str,
    site_identifier_records: list[SiteIdentifierNormalisationRecord],
) -> pd.Index:
    return _DEFAULT_REPORTER.normalize_supported_site_index_if_present(
        index,
        field_name=field_name,
        site_identifier_records=site_identifier_records,
    )


def _canonicalize_site_index_with_label_validation(
    index: pd.Index,
    *,
    field_name: str,
    site_identifier_records: list[SiteIdentifierNormalisationRecord],
    index_name: str | None = None,
) -> pd.Index:
    return _DEFAULT_REPORTER.canonicalize_site_index_with_label_validation(
        index,
        field_name=field_name,
        site_identifier_records=site_identifier_records,
        index_name=index_name,
    )

"""Internal convention normalisation for dataset builder inputs."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.frames.ownership import own_dataframe, own_optional_dataframe
from phospy.science.datasets.builders.normalization_reporter import (
    SAMPLE_LABEL_INDEX_POLICY as _SAMPLE_LABEL_INDEX_POLICY,
)
from phospy.science.datasets.builders.normalization_reporter import (
    DatasetConventionNormalisationReporter,
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
        normalized_phospho = own_dataframe(
            phospho,
            field_name="dataset build request phospho",
        )
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
            own_dataframe(
                site_metadata,
                field_name="dataset build request site_metadata",
            ),
            phospho_index=normalized_phospho.index,
            site_identifier_records=site_identifier_records,
        )
        normalized_sample_metadata = self._sample_metadata_normalizer.run(
            own_optional_dataframe(
                sample_metadata,
                field_name="dataset build request sample_metadata",
            ),
            phospho_columns=normalized_phospho.columns,
        )
        normalized_total = self._total_matrix_normalizer.run(
            own_optional_dataframe(
                total,
                field_name="dataset build request total",
            ),
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

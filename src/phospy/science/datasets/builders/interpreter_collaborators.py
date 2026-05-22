"""Internal collaborators for dataset build request interpretation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import NoReturn

import pandas as pd

from phospy.contracts.requests import DatasetBuildRequest
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.builders.contracts import DatasetInput
from phospy.science.datasets.builders.normalizer import (
    DatasetConventionNormalizer,
    NormalizedDatasetInputs,
)
from phospy.science.datasets.builders.reader import DatasetInputReader
from phospy.science.datasets.builders.sequence_derivation import SiteSequenceDeriver
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.datasets.preprocessing.policy_models import SiteMatrixPolicy
from phospy.science.evidence.dataset_resolution import (
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
    DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
    PeptideEvidenceDatasetResolver,
    build_multi_site_handling_config_for_dataset_policy,
)
from phospy.science.evidence.models import PeptideEvidenceTable
from phospy.science.references.models import Organism
from phospy.science.sites.identifiers import SiteIdentifierNormalisationReport


@dataclass(frozen=True, slots=True)
class ResolvedDatasetBuildSources:
    """Input tables resolved from a validated dataset build request."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    site_resolution_mode: str
    multi_site_policy: str | None
    peptide_evidence_resolution: dict[str, object] | None
    site_identifier_normalisation: SiteIdentifierNormalisationReport | None


@dataclass(frozen=True, slots=True)
class ResolvedSiteSequenceInputs:
    """Site metadata after sequence derivation and resolution payload collection."""

    site_metadata: pd.DataFrame
    site_sequence_derivation: dict[str, object] | None


class DatasetBuildSourceResolver:
    """Read and normalize input sources for dataset builder execution planning."""

    def __init__(
        self,
        *,
        reader: DatasetInputReader | None = None,
        normalizer: DatasetConventionNormalizer | None = None,
        peptide_evidence_resolver: PeptideEvidenceDatasetResolver | None = None,
    ) -> None:
        self._reader = reader or DatasetInputReader()
        self._normalizer = normalizer or DatasetConventionNormalizer()
        self._peptide_evidence_resolver = (
            peptide_evidence_resolver or PeptideEvidenceDatasetResolver()
        )

    def run(self, request: DatasetBuildRequest) -> ResolvedDatasetBuildSources:
        sample_metadata = self._read_optional(
            request.sample_metadata,
            field_name="sample_metadata",
        )
        total = self._read_optional(request.total, field_name="total")
        site_resolution_mode = str(request.site_resolution_mode).strip()
        multi_site_policy: str | None = None
        peptide_evidence_resolution_payload: dict[str, object] | None = None
        if site_resolution_mode == DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED:
            phospho = self._reader.run(
                _require_dataset_input(
                    request.phospho,
                    field_name="dataset build request phospho",
                ),
                field_name="phospho",
            )
            site_metadata = self._reader.run(
                _require_dataset_input(
                    request.site_metadata,
                    field_name="dataset build request site_metadata",
                ),
                field_name="site_metadata",
            )
        elif site_resolution_mode == DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE:
            (
                phospho,
                site_metadata,
                multi_site_policy,
                peptide_evidence_resolution_payload,
            ) = self._resolve_peptide_evidence_inputs(request)
        else:  # pragma: no cover - validator owns this branch; keep defensive.
            raise PhosPyInputError(
                "dataset build request site_resolution_mode is unsupported after "
                "validation"
            )
        normalized = self._normalize(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
        )
        return ResolvedDatasetBuildSources(
            phospho=normalized.phospho,
            site_metadata=normalized.site_metadata,
            sample_metadata=normalized.sample_metadata,
            total=normalized.total,
            site_resolution_mode=site_resolution_mode,
            multi_site_policy=multi_site_policy,
            peptide_evidence_resolution=peptide_evidence_resolution_payload,
            site_identifier_normalisation=normalized.site_identifier_normalisation,
        )

    def _read_optional(
        self,
        value: DatasetInput | None,
        *,
        field_name: str,
    ) -> pd.DataFrame | None:
        if value is None:
            return None
        return self._reader.run(value, field_name=field_name)

    def _resolve_peptide_evidence_inputs(
        self,
        request: DatasetBuildRequest,
    ) -> tuple[pd.DataFrame, pd.DataFrame, str, dict[str, object]]:
        peptide_evidence = self._reader.run(
            _require_dataset_input(
                request.peptide_evidence,
                field_name="dataset build request peptide_evidence",
            ),
            field_name="peptide_evidence",
        )
        peptide_site_mapping = self._read_optional(
            request.peptide_site_mapping,
            field_name="peptide_site_mapping",
        )
        multi_site_policy = str(request.multi_site_policy)
        try:
            evidence = PeptideEvidenceTable(
                frame=peptide_evidence,
                sample_intensity_columns=(
                    ()
                    if request.peptide_evidence_sample_intensity_columns is None
                    else request.peptide_evidence_sample_intensity_columns
                ),
                site_mapping=peptide_site_mapping,
                multi_site_handling_config=build_multi_site_handling_config_for_dataset_policy(
                    multi_site_policy=multi_site_policy
                ),
            )
            resolved = self._peptide_evidence_resolver.run(
                evidence=evidence,
                multi_site_policy=multi_site_policy,
            )
        except (TypeError, ValueError, KeyError, PhosPyInputError) as exc:
            _raise_wrapped_input_error(
                stage_name="dataset_builder.peptide_evidence_resolution",
                field_name="dataset build request peptide_evidence",
                operation=(
                    "resolving peptide evidence into site-level phospho and "
                    "site_metadata tables"
                ),
                next_action=(
                    "provide peptide_evidence with required evidence columns, valid "
                    "peptide_evidence_sample_intensity_columns, and a compatible "
                    "multi_site_policy/site mapping"
                ),
                original_error=exc,
            )
        return (
            resolved.phospho,
            resolved.site_metadata,
            multi_site_policy,
            resolved.summary.to_payload(),
        )

    def _normalize(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
    ) -> NormalizedDatasetInputs:
        try:
            return self._normalizer.run(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
            )
        except (TypeError, ValueError, KeyError) as exc:
            _raise_wrapped_input_error(
                stage_name="dataset_builder.normalization",
                field_name="dataset build request input tables",
                operation="normalizing input indices and metadata column conventions",
                next_action=(
                    "ensure phospho/site_metadata/sample_metadata/total tables use "
                    "supported rectangular DataFrame shapes, canonical site labels, "
                    "and non-conflicting metadata column conventions"
                ),
                original_error=exc,
            )


class DatasetBuildPreprocessingPlanner:
    """Resolve request preprocessing config into a concrete preprocessing plan."""

    def run(self, request: DatasetBuildRequest) -> PreprocessingPlan:
        return PreprocessingPlan.from_config(request.preprocessing_config)


class DatasetBuildSiteSequenceResolver:
    """Derive and validate site sequence metadata for execution planning."""

    def __init__(self, *, site_sequence_deriver: SiteSequenceDeriver | None = None):
        self._site_sequence_deriver = site_sequence_deriver or SiteSequenceDeriver()

    def run(
        self,
        *,
        site_metadata: pd.DataFrame,
        organism: Organism | None,
        preprocessing_plan: PreprocessingPlan,
    ) -> ResolvedSiteSequenceInputs:
        defer_missing_site_sequence_fill = bool(
            preprocessing_plan.site_sequence_resolution_enabled
        )
        try:
            enriched_site_metadata = self._site_sequence_deriver.run(
                site_metadata,
                organism=organism,
                allow_partial=(
                    preprocessing_plan.site_matrix_policy
                    is SiteMatrixPolicy.BUILD_FROM_METADATA
                    or preprocessing_plan.site_sequence_resolution_enabled
                ),
                derive_missing_from_reference=not defer_missing_site_sequence_fill,
            )
        except (TypeError, ValueError, KeyError) as exc:
            _raise_wrapped_input_error(
                stage_name="dataset_builder.site_sequence_derivation",
                field_name="dataset build request site_metadata",
                operation="deriving and validating site_sequence values",
                next_action=(
                    "ensure site_metadata contains canonical non-empty gene/site "
                    "fields and supported site identifiers for sequence derivation"
                ),
                original_error=exc,
            )
        return ResolvedSiteSequenceInputs(
            site_metadata=enriched_site_metadata,
            site_sequence_derivation=_resolve_site_sequence_derivation_payload(
                self._site_sequence_deriver
            ),
        )


def _resolve_site_sequence_derivation_payload(
    site_sequence_deriver: object,
) -> dict[str, object] | None:
    report = getattr(site_sequence_deriver, "last_report", None)
    if report is None:
        return None
    to_payload = getattr(report, "to_payload", None)
    if callable(to_payload):
        payload = to_payload()
        if isinstance(payload, dict):
            return payload
    return None


def _require_dataset_input(
    value: DatasetInput | None, *, field_name: str
) -> DatasetInput:
    if value is None:
        raise PhosPyInputError(f"{field_name} is required")
    return value


def _raise_wrapped_input_error(
    *,
    stage_name: str,
    field_name: str,
    operation: str,
    next_action: str,
    original_error: Exception,
) -> NoReturn:
    original_message = " ".join(str(original_error).split())
    raise PhosPyInputError(
        f"{stage_name} failed while {operation} for {field_name}. "
        f"Original error: {type(original_error).__name__}: {original_message}. "
        f"Next action: {next_action}"
    ) from original_error

"""Dataset construction planning collaborators."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.datasets.construction.contracts import DatasetBuildRequestProtocol
from phospy.datasets.construction.input_resolution import _raise_wrapped_input_error
from phospy.datasets.construction.sequence_derivation import (
    SiteSequenceDerivationReport,
    SiteSequenceDeriver,
)
from phospy.errors.build import DatasetBuildError
from phospy.preprocessing.contracts import PreprocessingPlan
from phospy.preprocessing.contracts.policies import SiteMatrixPolicy
from phospy.preprocessing.planning.boundary import PreprocessingPlanInterpreter
from phospy.preprocessing.workflow import GroupCoverageFilterMetadataValidator
from phospy.quantitative.models import IntensityScaleKind
from phospy.taxonomy.organisms import Organism


@dataclass(frozen=True, slots=True)
class ResolvedSiteSequenceInputs:
    """Site metadata after sequence derivation and report collection."""

    site_metadata: pd.DataFrame
    site_sequence_derivation: SiteSequenceDerivationReport | None


class DatasetBuildPreprocessingPlanner:
    """Resolve request preprocessing config into a concrete preprocessing plan."""

    def __init__(
        self,
        *,
        plan_interpreter: PreprocessingPlanInterpreter | None = None,
        group_coverage_validator: GroupCoverageFilterMetadataValidator | None = None,
    ) -> None:
        self._plan_interpreter = plan_interpreter or PreprocessingPlanInterpreter()
        self._group_coverage_validator = (
            group_coverage_validator or GroupCoverageFilterMetadataValidator()
        )

    def run(
        self,
        request: DatasetBuildRequestProtocol,
        *,
        phospho: pd.DataFrame | None = None,
        sample_metadata: pd.DataFrame | None = None,
        declared_input_intensity_scale_kind: IntensityScaleKind | None = None,
    ) -> PreprocessingPlan:
        plan = self._plan_interpreter.run(
            request.preprocessing_config,
            declared_input_scale_kind=declared_input_intensity_scale_kind,
        )
        if plan.group_coverage_filter_enabled and phospho is not None:
            self._group_coverage_validator.run(
                phospho=phospho,
                sample_metadata=sample_metadata,
                group_column=plan.group_coverage_filter_group_column,
                min_groups_passing_threshold=(
                    plan.group_coverage_filter_min_groups_passing_threshold
                ),
            )
        return plan


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
                    "ensure site_metadata contains non-empty gene/site "
                    "fields and supported site identifiers for sequence derivation"
                ),
                original_error=exc,
            )
        return ResolvedSiteSequenceInputs(
            site_metadata=enriched_site_metadata,
            site_sequence_derivation=_resolve_site_sequence_derivation_report(
                self._site_sequence_deriver
            ),
        )


def _resolve_site_sequence_derivation_report(
    site_sequence_deriver: object,
) -> SiteSequenceDerivationReport | None:
    report = getattr(site_sequence_deriver, "last_report", None)
    if report is None:
        return None
    if isinstance(report, SiteSequenceDerivationReport):
        return report
    raise DatasetBuildError(
        "dataset builder site-sequence derivation report must be a "
        "SiteSequenceDerivationReport when present; "
        f"got {type(report).__name__}"
    )

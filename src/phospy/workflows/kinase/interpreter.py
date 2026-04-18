"""Internal interpreter for kinase workflow requests."""

from __future__ import annotations

import pandas as pd

from phospy.api.requests import KinaseWorkflowRequest
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.references.resolution import (
    BundledReferenceProvider,
    ReferenceResolver,
    ReferenceResolverContract,
)
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest


class KinaseWorkflowInterpreter:
    """Resolve workflow request defaults and references for execution."""

    _KINASE_COLUMN = "kinase"
    _SUBSTRATE_COLUMN = "substrate_site"

    def __init__(
        self, *, reference_resolver: ReferenceResolverContract | None = None
    ) -> None:
        self._reference_resolver = reference_resolver or ReferenceResolver(
            provider=BundledReferenceProvider()
        )

    def run(self, request: KinaseWorkflowRequest) -> ResolvedKinaseWorkflowRequest:
        references = self._reference_resolver.run(
            request.references,
            dataset_organism=request.dataset.organism,
        )
        kinase_substrate_map = self._normalized_kinase_substrate_map(
            references.kinase_substrate_map
        )
        overlap_counts = self._summarize_overlap(
            dataset=request.dataset.phospho,
            kinase_substrate_map=kinase_substrate_map,
        )
        self._validate_reference_coverage(
            overlap_counts=overlap_counts,
            request=request,
        )
        self._validate_eligible_kinases(
            overlap_counts=overlap_counts,
            request=request,
        )
        return ResolvedKinaseWorkflowRequest(
            dataset=request.dataset,
            references=references,
            kinase_substrate_map=kinase_substrate_map,
            scoring_config=request.scoring_config,
            prediction_config=request.prediction_config,
            activity_config=request.activity_config,
        )

    @classmethod
    def _normalized_kinase_substrate_map(cls, mapping: pd.DataFrame) -> pd.DataFrame:
        cleaned = mapping[[cls._KINASE_COLUMN, cls._SUBSTRATE_COLUMN]].copy(deep=True)
        cleaned.loc[:, cls._KINASE_COLUMN] = (
            cleaned.loc[:, cls._KINASE_COLUMN].astype(str).str.strip()
        )
        cleaned.loc[:, cls._SUBSTRATE_COLUMN] = (
            cleaned.loc[:, cls._SUBSTRATE_COLUMN].astype(str).str.strip()
        )
        return cleaned.drop_duplicates(ignore_index=True)

    @classmethod
    def _summarize_overlap(
        cls,
        *,
        dataset: pd.DataFrame,
        kinase_substrate_map: pd.DataFrame,
    ) -> dict[str, int | pd.Series]:
        dataset_sites = {
            str(site_id).strip()
            for site_id in dataset.index
            if str(site_id).strip() != ""
        }
        reference_sites = set(
            kinase_substrate_map.loc[:, cls._SUBSTRATE_COLUMN].tolist()
        )
        overlapping_sites = dataset_sites.intersection(reference_sites)
        overlapping_map = kinase_substrate_map[
            kinase_substrate_map.loc[:, cls._SUBSTRATE_COLUMN].isin(overlapping_sites)
        ]
        per_kinase_quantified = (
            overlapping_map.groupby(cls._KINASE_COLUMN, sort=False)[
                cls._SUBSTRATE_COLUMN
            ]
            .nunique()
            .astype("int64")
        )
        return {
            "dataset_sites": len(dataset_sites),
            "reference_sites": len(reference_sites),
            "overlap_sites": len(overlapping_sites),
            "reference_kinases": int(
                kinase_substrate_map.loc[:, cls._KINASE_COLUMN].nunique()
            ),
            "kinases_with_overlap": int(per_kinase_quantified.size),
            "max_quantified_sites_per_kinase": int(
                per_kinase_quantified.max() if not per_kinase_quantified.empty else 0
            ),
            "per_kinase_quantified": per_kinase_quantified,
        }

    def _validate_reference_coverage(
        self,
        *,
        overlap_counts: dict[str, int | pd.Series],
        request: KinaseWorkflowRequest,
    ) -> None:
        overlap_sites = int(overlap_counts["overlap_sites"])
        if overlap_sites > 0:
            return
        self._raise_boundary_error(
            seam="kinase.interpreter.reference_coverage",
            next_action=(
                "use references that contain dataset phosphosite IDs or verify site "
                "identifier formatting in dataset.phospho.index"
            ),
            dataset_sites=overlap_counts["dataset_sites"],
            reference_sites=overlap_counts["reference_sites"],
            overlap_sites=overlap_sites,
            scoring_config_min_substrates=request.scoring_config.min_substrates,
        )

    def _validate_eligible_kinases(
        self,
        *,
        overlap_counts: dict[str, int | pd.Series],
        request: KinaseWorkflowRequest,
    ) -> None:
        per_kinase_quantified = overlap_counts["per_kinase_quantified"]
        assert isinstance(per_kinase_quantified, pd.Series)
        eligible_kinases = per_kinase_quantified[
            per_kinase_quantified >= request.scoring_config.min_substrates
        ]
        if not eligible_kinases.empty:
            return
        self._raise_boundary_error(
            seam="kinase.interpreter.eligible_kinases",
            next_action=(
                "lower scoring_config.min_substrates or provide references with "
                "deeper overlap for the current dataset"
            ),
            reference_kinases=overlap_counts["reference_kinases"],
            kinases_with_overlap=overlap_counts["kinases_with_overlap"],
            eligible_kinases=int(eligible_kinases.size),
            max_quantified_sites_per_kinase=overlap_counts[
                "max_quantified_sites_per_kinase"
            ],
            scoring_config_min_substrates=request.scoring_config.min_substrates,
            prediction_config_ensemble_size=request.prediction_config.ensemble_size,
        )

    @staticmethod
    def _raise_boundary_error(
        *,
        seam: str,
        next_action: str,
        **details: int,
    ) -> None:
        details_text = ", ".join(f"{key}={value}" for key, value in details.items())
        raise WorkflowBoundaryError(
            "kinase workflow boundary validation failed at "
            f"seam={seam}; {details_text}; next_action={next_action}"
        )

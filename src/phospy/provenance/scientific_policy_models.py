"""Stable scientific-policy identifiers and serializable metadata records."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum

from phospy.errors.input import PhosPyInputError
from phospy.provenance.immutability import freeze_json_mapping, thaw_json_mapping

ScientificPolicyParameter = str | int | float | bool | None


class ScientificPolicyId(str, Enum):
    """Stable identifiers for emitted scientific-policy records.

    Enum membership preserves payload readability; it is not by itself a public
    support claim for the producing route.
    """

    PROFILE_CORRELATION_SHIFTED_UNIT = "profile_correlation_shifted_unit_v1"
    KINASE_PROFILE_SCORING = "kinase_profile_scoring_v1"
    KINASE_LIBRARY_MOTIF_SCORING = "kinase_library_motif_scoring_v1"
    MOTIF_PROFILE_RANK_FUSION = "motif_profile_rank_fusion_v1"
    CANDIDATE_SUBSTRATE_SELECTION = "candidate_substrate_selection_v1"
    SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY = "simplified_weighted_substrate_activity_v1"
    KSEA_ZSCORE_ACTIVITY = "ksea_zscore_activity_v1"
    SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY = "ssgsea_substrate_enrichment_activity_v1"
    SIGNALOME_MISSING_VALUE_CLUSTERING = "signalome_missing_value_clustering_v1"
    SIGNALOME_SCORE_PRECONDITIONING = "signalome_score_preconditioning_v1"
    PREPROCESSING_STAGE_ORDER = "preprocessing_stage_order_v1"
    SIGNALOME_MODULE_CANDIDATE_SCORE = "signalome_module_candidate_score_v1"
    PROTEIN_MODULE_FROM_SITE_MEMBERSHIP = "protein_module_from_site_membership_v1"
    DUPLICATE_SITE_RESOLUTION = "duplicate_site_resolution_v1"
    ADAPTIVE_PREDICTION_SAMPLING = "adaptive_prediction_sampling_v1"
    SIGNALOME_DOWNSTREAM_SCORE_SELECTION = "signalome_downstream_score_selection_v1"
    SIGNALOME_CANDIDATE_SCORING = "signalome_candidate_scoring_v1"
    SIGNALOME_ASSIGNMENT_POLICY = "signalome_assignment_policy_v1"
    SIGNALOME_NETWORK_POLICY = "signalome_network_policy_v1"
    PEPTIDE_TO_SITE_AGGREGATION = "peptide_to_site_aggregation_v1"


@dataclass(frozen=True, slots=True)
class ScientificPolicyRecord:
    """Serializable metadata for one scientific scoring/derivation policy."""

    id: ScientificPolicyId
    name: str
    version: str
    description: str
    parameters: Mapping[str, object]
    assumptions: tuple[str, ...]
    output_scale: str | None = None
    quantitative_meaning: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "parameters",
            freeze_json_mapping(
                self.parameters,
                field_name="scientific_policy_record.parameters",
            ),
        )
        object.__setattr__(
            self,
            "assumptions",
            tuple(str(value) for value in self.assumptions),
        )

    def to_payload(self) -> dict[str, object]:
        return {
            "id": self.id.value,
            "name": self.name,
            "version": self.version,
            "description": self.description,
            "parameters": thaw_json_mapping(
                self.parameters,
                field_name="scientific_policy_record.parameters",
            ),
            "assumptions": list(self.assumptions),
            "output_scale": self.output_scale,
            "quantitative_meaning": self.quantitative_meaning,
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> ScientificPolicyRecord:
        payload = _require_payload_mapping(
            payload,
            field_name="scientific_policy_record",
        )
        assumptions = payload.get("assumptions", ())
        if isinstance(assumptions, (list, tuple)):
            resolved_assumptions = tuple(str(value) for value in assumptions)
        else:
            resolved_assumptions = ()
        parameters_raw = payload.get("parameters", {})
        if not isinstance(parameters_raw, Mapping):
            raise PhosPyInputError(
                "scientific_policy_record.parameters must be an object"
            )
        output_scale = payload.get("output_scale")
        resolved_output_scale = None if output_scale is None else str(output_scale)
        quantitative_meaning = payload.get("quantitative_meaning")
        resolved_quantitative_meaning = (
            None if quantitative_meaning is None else str(quantitative_meaning)
        )
        return cls(
            id=ScientificPolicyId(str(payload.get("id"))),
            name=str(payload.get("name")),
            version=str(payload.get("version")),
            description=str(payload.get("description")),
            parameters=parameters_raw,
            assumptions=resolved_assumptions,
            output_scale=resolved_output_scale,
            quantitative_meaning=resolved_quantitative_meaning,
        )


def _require_payload_mapping(
    value: Mapping[str, object],
    *,
    field_name: str,
) -> Mapping[str, object]:
    result: dict[str, object] = {}
    for key, item in value.items():
        if not isinstance(key, str):
            raise PhosPyInputError(
                f"{field_name} JSON object keys must be strings; "
                f"got {type(key).__name__}"
            )
        if key in result:
            raise PhosPyInputError(
                f"{field_name} contains duplicate JSON object key {key!r}"
            )
        result[key] = item
    return result


__all__ = [
    "ScientificPolicyId",
    "ScientificPolicyParameter",
    "ScientificPolicyRecord",
]

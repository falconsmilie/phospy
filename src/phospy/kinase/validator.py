"""Internal validator for kinase workflow requests."""

from __future__ import annotations

import math
from typing import cast

import numpy as np
import pandas as pd

from phospy.activities.validation import (
    MethodQuantitativeInputValidator,
)
from phospy.datasets.internal_view import DatasetInternalView
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.validation import WorkflowValidationError
from phospy.evidence.site_metadata import (
    enforce_localisation_requirement,
)
from phospy.frames.validation import require_dataframe
from phospy.kinase.activity.config import KinaseActivityConfig
from phospy.kinase.activity.contracts import (
    kinase_activity_method_quantitative_input_contract,
)
from phospy.kinase.contracts import ValidatedKinaseWorkflowRequest
from phospy.kinase.prediction.config import KinasePredictionConfig
from phospy.kinase.references.models import KinaseLibraryResource
from phospy.kinase.requests import KinaseWorkflowRequest
from phospy.kinase.scoring.config import (
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
    KinaseScoringConfig,
)
from phospy.kinase.scoring.contracts import (
    kinase_scoring_method_quantitative_input_contract,
    kinase_scoring_mode_input_contract,
)
from phospy.kinase.sequence_contracts import (
    dataset_sequence_source_label,
    kinase_sequence_context_contract,
)
from phospy.kinase.site_sequence_policy import (
    resolve_site_sequence_conflict_policy,
)
from phospy.phosphosites.contracts import WORKFLOW_SEQUENCE_CONTEXT_IDENTITY_CONTRACT
from phospy.phosphosites.sequence_context import (
    enforce_site_sequence_context_contract,
)
from phospy.phosphosites.workflow_identity import (
    WorkflowIdentityContract,
    enforce_workflow_site_identity_contract,
)
from phospy.quantitative.validation import (
    WorkflowQuantitativeInputContract,
    phosphosite_abundance_workflow_input_contract,
)
from phospy.references.config import ReferenceContextCompatibilityPolicy
from phospy.references.context_compatibility import (
    validate_reference_context_compatibility,
)
from phospy.references.models import (
    ReferenceBundle,
    ReferencePreset,
)
from phospy.references.provenance.context import ReferenceContext

KINASE_IDENTITY_CONTRACT = WorkflowIdentityContract(
    workflow_name="kinase workflow request",
    contract_id="sty_site_identity_plus_sequence_context",
    identity_contract=WORKFLOW_SEQUENCE_CONTEXT_IDENTITY_CONTRACT,
    allow_gapped_sequence_context=True,
)


def kinase_profile_scoring_workflow_input_contract(
    *,
    allow_mixed_total_protein_quantitative_meaning: bool = False,
) -> WorkflowQuantitativeInputContract:
    """Return the quantitative contract for Kinase profile-scoring inputs."""

    return phosphosite_abundance_workflow_input_contract(
        allow_mixed_total_protein_quantitative_meaning=(
            allow_mixed_total_protein_quantitative_meaning
        )
    )


class KinaseWorkflowConfigValidator:
    """Validate Kinase workflow config objects and local invariants."""

    def run(
        self,
        *,
        scoring_config: object,
        prediction_config: object,
        activity_config: object | None,
    ) -> tuple[
        KinaseScoringConfig,
        KinasePredictionConfig,
        KinaseActivityConfig | None,
    ]:
        validated_scoring = self._validate_scoring(scoring_config)
        validated_prediction = self._validate_prediction(prediction_config)
        validated_activity = self._validate_activity(activity_config)
        return validated_scoring, validated_prediction, validated_activity

    @staticmethod
    def _validate_scoring(config: object) -> KinaseScoringConfig:
        if config is None:
            raise WorkflowValidationError(
                "kinase workflow request scoring_config is required; choose "
                "KinaseScoringConfig.exploratory(), "
                "KinaseScoringConfig.production(...), or direct construction with "
                "reliability_profile='custom'"
            )
        if not isinstance(config, KinaseScoringConfig):
            raise WorkflowValidationError(
                "kinase workflow request scoring_config must be KinaseScoringConfig"
            )
        return config

    @staticmethod
    def _validate_prediction(config: object) -> KinasePredictionConfig:
        if not isinstance(config, KinasePredictionConfig):
            raise WorkflowValidationError(
                "kinase workflow request prediction_config must be "
                "KinasePredictionConfig"
            )
        return config

    @staticmethod
    def _validate_activity(config: object | None) -> KinaseActivityConfig | None:
        if config is None:
            return None
        if not isinstance(config, KinaseActivityConfig):
            raise WorkflowValidationError(
                "kinase workflow request activity_config must be "
                "KinaseActivityConfig or None"
            )
        return config


class KinaseWorkflowValidator:
    """Validate `KinaseWorkflowRequest` before interpretation."""

    def __init__(
        self,
        *,
        config_validator: KinaseWorkflowConfigValidator | None = None,
        method_quantitative_input_validator: (
            MethodQuantitativeInputValidator | None
        ) = None,
    ) -> None:
        self._config_validator = config_validator or KinaseWorkflowConfigValidator()
        self._method_quantitative_input_validator = (
            method_quantitative_input_validator or MethodQuantitativeInputValidator()
        )

    def run(self, request: object) -> ValidatedKinaseWorkflowRequest:
        if not isinstance(request, KinaseWorkflowRequest):
            raise WorkflowValidationError(
                "kinase workflow input must be a KinaseWorkflowRequest"
            )
        dataset = cast(object, request.dataset)
        if not isinstance(dataset, AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                "kinase workflow request dataset must be AnalysisReadyPhosphoDataset"
            )
        dataset_view = DatasetInternalView(dataset)
        references = cast(object, request.references)
        if not isinstance(references, (ReferencePreset, ReferenceBundle)):
            raise WorkflowValidationError(
                "kinase workflow request references must be ReferencePreset or ReferenceBundle"
            )
        scoring_config, _, activity_config = self._config_validator.run(
            scoring_config=request.scoring_config,
            prediction_config=request.prediction_config,
            activity_config=request.activity_config,
        )
        if isinstance(references, ReferenceBundle):
            validate_reference_context_compatibility(
                dataset.reference_context,
                _reference_bundle_context(references),
                operation="kinase workflow request dataset/reference bundle",
                allow_unknown=(
                    scoring_config.reference_context_compatibility_policy
                    is ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
                error_type=WorkflowValidationError,
            )
        if (
            request.reference_display_ambiguity_policy
            not in KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES
        ):
            supported = ", ".join(sorted(KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES))
            raise WorkflowValidationError(
                "kinase workflow request reference_display_ambiguity_policy "
                f"must be one of: {supported}"
            )
        site_sequence_conflict_policy = resolve_site_sequence_conflict_policy(
            request.site_sequence_conflict_policy,
            field_name="kinase workflow request site_sequence_conflict_policy",
            error_type=WorkflowValidationError,
        )
        mode_contract = kinase_scoring_mode_input_contract(scoring_config.scoring_mode)
        if mode_contract.requires_kinase_library_resource:
            if request.kinase_library_resource is None:
                raise WorkflowValidationError(
                    "kinase workflow request kinase_library_resource is required "
                    f"when scoring_config.scoring_mode={scoring_config.scoring_mode!r}"
                )
            if not isinstance(request.kinase_library_resource, KinaseLibraryResource):
                raise WorkflowValidationError(
                    "kinase workflow request kinase_library_resource must be "
                    "KinaseLibraryResource when Kinase Library scoring is selected"
                )
        self._method_quantitative_input_validator.run(
            intensity_scale_state=dataset.intensity_scale_state,
            has_total_matrix=dataset_view.total is not None,
            contract=kinase_scoring_method_quantitative_input_contract(
                scoring_config.scoring_mode,
                allow_mixed_total_protein_quantitative_meaning=(
                    scoring_config.allow_mixed_total_protein_quantitative_meaning
                ),
            ),
            context="kinase workflow request dataset",
        )
        if activity_config is not None and activity_config.enabled:
            self._method_quantitative_input_validator.run(
                intensity_scale_state=dataset.intensity_scale_state,
                has_total_matrix=dataset_view.total is not None,
                contract=kinase_activity_method_quantitative_input_contract(
                    activity_config.method
                ),
                context="kinase activity request dataset",
            )
        site_metadata = require_dataframe(
            dataset_view.site_metadata,
            field_name="kinase workflow request dataset.site_metadata",
            allow_empty=False,
            error_type=WorkflowValidationError,
        )
        enforce_workflow_site_identity_contract(
            site_metadata=site_metadata,
            expected_index=dataset_view.phospho.index,
            expected_index_field_name="kinase workflow request dataset.phospho.index",
            field_name="kinase workflow request dataset.site_metadata",
            contract=KINASE_IDENTITY_CONTRACT,
            error_type=WorkflowValidationError,
            allow_opaque_site_values=dataset.opaque_site_values_allowed,
            scoring_mode=scoring_config.scoring_mode,
        )
        sequence_contract = kinase_sequence_context_contract(
            scoring_mode=scoring_config.scoring_mode,
            kinase_library_resource=request.kinase_library_resource,
        )
        if sequence_contract is not None and isinstance(references, ReferenceBundle):
            selected_site_sequences, sequence_source_by_site = (
                _selected_explicit_reference_sequence_context(
                    dataset=dataset,
                    site_metadata=site_metadata,
                    references=references,
                    conflict_policy=site_sequence_conflict_policy,
                )
            )
            enforce_site_sequence_context_contract(
                site_metadata=selected_site_sequences,
                field_name="kinase workflow request selected_site_sequences",
                workflow_name="kinase workflow request",
                scoring_mode=scoring_config.scoring_mode,
                contract=sequence_contract,
                error_type=WorkflowValidationError,
                sequence_source_by_site=sequence_source_by_site,
                allow_unknown_site_residue=False,
            )
        enforce_localisation_requirement(
            site_metadata=site_metadata,
            field_name="kinase workflow request dataset.site_metadata",
            workflow_name="kinase workflow request",
            requirement=scoring_config.localisation_requirement,
            error_type=WorkflowValidationError,
        )
        return ValidatedKinaseWorkflowRequest(
            request=request,
            dataset_view=dataset_view,
        )


def _selected_explicit_reference_sequence_context(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    site_metadata: pd.DataFrame,
    references: ReferenceBundle,
    conflict_policy: object,
) -> tuple[pd.DataFrame, dict[str, str]]:
    resolved_policy = resolve_site_sequence_conflict_policy(
        conflict_policy,
        field_name="kinase workflow request site_sequence_conflict_policy",
        error_type=WorkflowValidationError,
    )
    source_by_site: dict[str, str] = {}
    reference_sequences = {
        str(display_id).strip(): sequence
        for display_id, sequence in (
            (
                display_id,
                _normalise_sequence_value(value),
            )
            for display_id, value in references.site_sequences.loc[
                :, "site_sequence"
            ].items()
        )
        if str(display_id).strip() != "" and sequence is not None
    }
    dataset_source_label = dataset_sequence_source_label(dataset)
    selected_sequence_values: list[object] = []
    site_ids = site_metadata.index.tolist()
    display_ids = site_metadata.loc[:, "display_id"].tolist()
    dataset_sequences = site_metadata.loc[:, "site_sequence"].tolist()
    for site_id, raw_display_id, raw_dataset_sequence in zip(
        site_ids,
        display_ids,
        dataset_sequences,
        strict=True,
    ):
        site_key = str(site_id)
        display_id = str(raw_display_id).strip()
        dataset_sequence_text = _normalise_sequence_value(raw_dataset_sequence)
        reference_sequence_text = reference_sequences.get(display_id)
        if reference_sequence_text is None:
            if dataset_sequence_text is not None:
                source_by_site[site_key] = dataset_source_label or "unknown"
                selected_sequence_values.append(dataset_sequence_text)
            else:
                selected_sequence_values.append(raw_dataset_sequence)
            continue
        if dataset_sequence_text is None:
            source_by_site[site_key] = "reference"
            selected_sequence_values.append(reference_sequence_text)
            continue
        if reference_sequence_text == dataset_sequence_text:
            source_by_site[site_key] = "reference"
            selected_sequence_values.append(reference_sequence_text)
            continue
        if resolved_policy == KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR:
            raise WorkflowValidationError(
                "kinase workflow request selected sequence context failed; "
                "dataset sequence and reference sequence are incompatible and no "
                "conflict policy resolves them; "
                f"site_key={site_key!r}; display_id={display_id!r}; "
                f"dataset_sequence={dataset_sequence_text!r}; "
                f"reference_sequence={reference_sequence_text!r}; "
                f"conflict_policy={str(resolved_policy)!r}"
            )
        if resolved_policy == KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET:
            source_by_site[site_key] = dataset_source_label or "unknown"
            selected_sequence_values.append(dataset_sequence_text)
            continue
        if resolved_policy == KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE:
            source_by_site[site_key] = "reference"
            selected_sequence_values.append(reference_sequence_text)
            continue
        selected_sequence_values.append(raw_dataset_sequence)
    # `assign` constructs the selected metadata frame without mutating a shallow
    # copy that may still share blocks with the dataset's internal view.
    selected_sequences = pd.Series(
        selected_sequence_values,
        index=site_metadata.index,
        dtype=object,
    )
    selected = site_metadata.assign(site_sequence=selected_sequences)
    return selected, source_by_site


def _normalise_sequence_value(value: object) -> str | None:
    if _is_missing(value):
        return None
    text = str(value).strip()
    return text or None


def _reference_bundle_context(
    references: ReferenceBundle,
) -> ReferenceContext | None:
    provenance = references.provenance
    if provenance is None:
        return None
    return provenance.reference_context


def _is_missing(value: object) -> bool:
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, np.floating):
        scalar_value: object = value
        return str(scalar_value).lower() == "nan"
    if isinstance(value, (np.datetime64, np.timedelta64)):
        temporal_value: object = value
        return str(temporal_value) == "NaT"
    return False

"""Internal validator for signalome workflow requests."""

from __future__ import annotations

import pandas as pd

from phospy.datasets.internal_view import DatasetInternalView
from phospy.errors.validation import WorkflowValidationError
from phospy.evidence.config import (
    LOCALISATION_PRODUCTION_MINIMUM_PROBABILITY,
    LocalisationRequirement,
)
from phospy.evidence.site_metadata import enforce_localisation_requirement
from phospy.kinase.result_boundary import DownstreamScoreSource
from phospy.kinase.results import KinaseWorkflowResult
from phospy.phosphosites.contracts import (
    WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT,
    WORKFLOW_PROTEIN_SEQUENCE_CONTEXT_IDENTITY_CONTRACT,
)
from phospy.phosphosites.workflow_identity import (
    WorkflowIdentityContract,
    enforce_workflow_site_identity_contract,
)
from phospy.quantitative.models import QuantitativeMeaning
from phospy.quantitative.validation import (
    WorkflowQuantitativeInputContract,
    WorkflowQuantitativeInputValidator,
    is_mixed_quantitative_meaning,
    phosphosite_abundance_workflow_input_contract,
)
from phospy.references.config import ReferenceContextCompatibilityPolicy
from phospy.references.context_compatibility import (
    validate_reference_context_compatibility,
)
from phospy.references.provenance.context import ReferenceContext
from phospy.signalome.config import (
    SIGNALOME_MODE_PRODUCTION,
    SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT,
    SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_FLOOR,
    SignalomeConfig,
)
from phospy.signalome.constants import (
    LEGACY_PROTEIN_GROUP_ID_COLUMN,
    PROTEIN_GROUP_ID_COLUMN,
)
from phospy.signalome.contracts import ValidatedSignalomeWorkflowRequest
from phospy.signalome.requests import SignalomeWorkflowRequest
from phospy.validation.common.dataframes import (
    require_dataframe,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.validation.frame_columns import enforce_required_non_empty_string_column

SIGNALOME_IDENTITY_CONTRACT = WorkflowIdentityContract(
    workflow_name="signalome workflow request",
    contract_id="protein_scoped_site_identity",
    identity_contract=WORKFLOW_PROTEIN_SEQUENCE_CONTEXT_IDENTITY_CONTRACT,
    allow_gapped_sequence_context=True,
)

SIGNALOME_PROTEIN_GROUP_ID_COLUMN = PROTEIN_GROUP_ID_COLUMN
SIGNALOME_LEGACY_PROTEIN_GROUP_ID_COLUMN = LEGACY_PROTEIN_GROUP_ID_COLUMN

SIGNALOME_PROTEIN_GROUPING_METADATA_NOTE = (
    "Signalome uses dataset.site_metadata.protein_group_id as explicit "
    "signalome-specific protein grouping metadata for protein-level module and "
    "protein-site context summaries. The legacy dataset.site_metadata.protein_id "
    "column is accepted only as a migration alias. This grouping value is not "
    "core protein identity. "
    "Core protein identity remains the dataset-level organism, "
    "protein_namespace, and protein_identifier metadata carried by the "
    "site_key-indexed row identity contract. "
    "Do not replace protein_namespace or protein_identifier with "
    "protein_group_id. "
    "Signalome does not infer protein_group_id from gene_symbol or display_id and "
    "does not repair invalid site_key identity."
)

SIGNALOME_SITE_IDENTITY_CONTRACT_NOTE = (
    "Base dataset phosphosite identity must already be valid and site_key-indexed "
    "before signalome validation; signalome does not repair weak dataset identity."
)


def signalome_workflow_input_contract(
    *,
    allow_mixed_total_protein_quantitative_meaning: bool = False,
) -> WorkflowQuantitativeInputContract:
    """Return the quantitative contract for Signalome dataset inputs."""

    return phosphosite_abundance_workflow_input_contract(
        allow_mixed_total_protein_quantitative_meaning=(
            allow_mixed_total_protein_quantitative_meaning
        )
    )


_MIXED_QUANTITATIVE_MEANING = QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE.value


class SignalomeConfigValidator:
    """Validate Signalome workflow config objects and local invariants."""

    def run(self, config: object) -> SignalomeConfig:
        if not isinstance(config, SignalomeConfig):
            raise WorkflowValidationError(
                "signalome workflow request config must be SignalomeConfig"
            )
        network_minimum = config.output.network_min_paired_finite_observations
        if (
            network_minimum is not None
            and int(network_minimum)
            < SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_FLOOR
        ):
            raise WorkflowValidationError(
                "signalome workflow request config.output."
                "network_min_paired_finite_observations must be at least "
                f"{SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_FLOOR}; "
                f"got {int(network_minimum)}. Historical signalome bundles with "
                "threshold 2 remain readable, but replay/re-execution must migrate "
                "to config.output.network_min_paired_finite_observations >= 3."
            )
        if config.mode == SIGNALOME_MODE_PRODUCTION:
            effective_minimum = (
                SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT
                if network_minimum is None
                else int(network_minimum)
            )
            if (
                effective_minimum
                < SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT
            ):
                raise WorkflowValidationError(
                    "signalome workflow request config.mode='production' requires "
                    "config.output.network_min_paired_finite_observations to "
                    "resolve to at least "
                    f"{SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT}; "
                    f"got {effective_minimum}"
                )
            localisation = config.validation.localisation_requirement
            if (
                not localisation.require_present
                or localisation.minimum_probability is None
                or float(localisation.minimum_probability)
                < LOCALISATION_PRODUCTION_MINIMUM_PROBABILITY
            ):
                raise WorkflowValidationError(
                    "signalome workflow request config.mode='production' requires "
                    "site-level localisation with minimum_probability >= "
                    f"{LOCALISATION_PRODUCTION_MINIMUM_PROBABILITY}"
                )
        return config


class SignalomeWorkflowValidator:
    """Validate `SignalomeWorkflowRequest` before interpretation."""

    def __init__(
        self,
        *,
        config_validator: SignalomeConfigValidator | None = None,
        quantitative_input_validator: WorkflowQuantitativeInputValidator | None = None,
    ) -> None:
        self._config_validator = config_validator or SignalomeConfigValidator()
        self._quantitative_input_validator = (
            quantitative_input_validator or WorkflowQuantitativeInputValidator()
        )

    def run(
        self, request: SignalomeWorkflowRequest
    ) -> ValidatedSignalomeWorkflowRequest:
        if not isinstance(request, SignalomeWorkflowRequest):
            raise WorkflowValidationError(
                "signalome workflow input must be a SignalomeWorkflowRequest"
            )
        if not isinstance(request.kinase_result, KinaseWorkflowResult):
            raise WorkflowValidationError(
                "signalome workflow request kinase_result must be KinaseWorkflowResult"
            )
        config = self._config_validator.run(request.config)
        dataset = request.kinase_result.dataset
        _require_reference_context_compatibility(
            request.kinase_result,
            policy=config.validation.reference_context_compatibility_policy,
        )
        dataset_view = DatasetInternalView(dataset)
        prediction_result = request.kinase_result.prediction_result
        scoring_result = request.kinase_result.scoring_result
        if (
            is_mixed_quantitative_meaning(dataset.intensity_scale_state.quantity)
            and not config.validation.allow_mixed_total_protein_quantitative_meaning
        ):
            context = "signalome workflow request kinase_result.dataset"
            contract = phosphosite_abundance_workflow_input_contract(
                allow_mixed_total_protein_quantitative_meaning=False
            )
            expected = ", ".join(
                repr(value.value)
                for value in sorted(
                    contract.allowed_meanings,
                    key=lambda item: item.value,
                )
            )
            diagnostics = dataset.processing_state.total_protein_correction.diagnostics
            uncorrected_row_count = None
            unmatched_policy = None
            if diagnostics is not None:
                uncorrected_row_count = diagnostics.get("uncorrected_row_count")
                unmatched_policy = diagnostics.get("unmatched_policy")
            raise WorkflowValidationError(
                f"{context} requires quantitative meaning in {{{expected}}}; "
                f"got {_MIXED_QUANTITATIVE_MEANING!r}. "
                f"{context} received a dataset with mixed total-protein "
                f"quantitative meaning ({_MIXED_QUANTITATIVE_MEANING}). "
                f"uncorrected_rows={uncorrected_row_count!r}, "
                f"unmatched_policy={unmatched_policy!r}. "
                "This usually happens when total-protein correction is applied "
                "with unmatched_policy='allow_uncorrected'. "
                "Recommended actions: set unmatched_policy='error' or complete "
                "the phosphosite-to-total mapping. "
                "If you intentionally want this mixed dataset, set the workflow "
                "allow_mixed_total_protein_quantitative_meaning option to True."
            )
        self._quantitative_input_validator.run(
            intensity_scale_state=dataset.intensity_scale_state,
            has_total_matrix=dataset_view.total is not None,
            contract=signalome_workflow_input_contract(
                allow_mixed_total_protein_quantitative_meaning=(
                    config.validation.allow_mixed_total_protein_quantitative_meaning
                )
            ),
            context="signalome workflow request kinase_result.dataset",
        )

        prediction_matrix = require_dataframe(
            prediction_result.pred_mat,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            allow_empty=True,
            error_type=WorkflowValidationError,
        )
        require_non_empty_dataframe(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            error_type=WorkflowValidationError,
        )
        require_numeric_dataframe(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            error_type=WorkflowValidationError,
        )
        require_finite_numeric_dataframe(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            allow_missing=True,
            error_type=WorkflowValidationError,
        )
        prediction_matrix = require_unique_index(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            error_type=WorkflowValidationError,
        )
        prediction_matrix = require_unique_columns(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            error_type=WorkflowValidationError,
        )
        if prediction_matrix.shape[1] == 0:
            raise WorkflowValidationError(
                "signalome workflow request kinase_result.prediction_result.pred_mat "
                "must contain at least one kinase column"
            )

        downstream_score_source = DownstreamScoreSource.parse(
            scoring_result.score_source,
            field_name="signalome workflow request downstream score source",
        )
        downstream_score_matrix = scoring_result.authoritative_scores
        score_field_name = (
            "signalome workflow request kinase_result.scoring_result."
            f"{downstream_score_source}"
        )
        score_matrix = require_dataframe(
            downstream_score_matrix,
            field_name=score_field_name,
            allow_empty=True,
            error_type=WorkflowValidationError,
        )
        require_non_empty_dataframe(
            score_matrix,
            field_name=score_field_name,
            error_type=WorkflowValidationError,
        )
        require_numeric_dataframe(
            score_matrix,
            field_name=score_field_name,
            error_type=WorkflowValidationError,
        )
        require_finite_numeric_dataframe(
            score_matrix,
            field_name=score_field_name,
            # Correlation-based kinase scoring can legitimately emit missing
            # values (for example, zero-variance denominator collapse).
            # Missingness is preconditioned downstream by the interpreter.
            allow_missing=True,
            error_type=WorkflowValidationError,
        )
        score_matrix = require_unique_index(
            score_matrix,
            field_name=score_field_name,
            error_type=WorkflowValidationError,
        )
        score_matrix = require_unique_columns(
            score_matrix,
            field_name=score_field_name,
            error_type=WorkflowValidationError,
        )
        if score_matrix.shape[1] == 0:
            raise WorkflowValidationError(
                f"{score_field_name} must contain at least one kinase column"
            )
        self._require_network_min_paired_observations_compatible(
            config=config,
            score_matrix=score_matrix,
            score_field_name=score_field_name,
        )
        self._require_site_identity_and_protein_grouping_metadata(
            site_metadata=dataset_view.site_metadata,
            dataset_sites=dataset_view.phospho.index,
            prediction_sites=prediction_matrix.index,
            score_sites=score_matrix.index,
            score_field_name=score_field_name,
            localisation_requirement=config.validation.localisation_requirement,
            allow_opaque_site_values=dataset.opaque_site_values_allowed,
        )
        return ValidatedSignalomeWorkflowRequest(
            request=request,
            dataset_view=dataset_view,
            prediction_matrix=prediction_matrix,
            authoritative_score_matrix=score_matrix,
        )

    @staticmethod
    def _require_network_min_paired_observations_compatible(
        *,
        config: SignalomeConfig,
        score_matrix: pd.DataFrame,
        score_field_name: str,
    ) -> None:
        threshold = config.output.network_min_paired_finite_observations
        if threshold is None:
            return
        available_observations = int(score_matrix.shape[0])
        if int(threshold) <= available_observations:
            return
        raise WorkflowValidationError(
            "signalome workflow request config.output."
            "network_min_paired_finite_observations "
            f"({int(threshold)}) cannot exceed available downstream score "
            f"observations ({available_observations}) in {score_field_name}"
        )

    @staticmethod
    def _require_site_identity_and_protein_grouping_metadata(
        *,
        site_metadata: object,
        dataset_sites: pd.Index,
        prediction_sites: pd.Index,
        score_sites: pd.Index,
        score_field_name: str,
        localisation_requirement: LocalisationRequirement,
        allow_opaque_site_values: bool,
    ) -> None:
        field_name = "signalome workflow request kinase_result.dataset.site_metadata"
        site_metadata_frame = require_dataframe(
            site_metadata,
            field_name=field_name,
            allow_empty=False,
            error_type=WorkflowValidationError,
        )
        try:
            enforce_workflow_site_identity_contract(
                site_metadata=site_metadata_frame,
                expected_index=dataset_sites,
                expected_index_field_name=(
                    "signalome workflow request kinase_result.dataset.phospho.index"
                ),
                field_name=field_name,
                contract=SIGNALOME_IDENTITY_CONTRACT,
                error_type=WorkflowValidationError,
                allow_opaque_site_values=allow_opaque_site_values,
                sequence_context_contract=WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT,
                scoring_mode="signalome",
            )
        except WorkflowValidationError as exc:
            raise WorkflowValidationError(
                f"{exc}. {SIGNALOME_SITE_IDENTITY_CONTRACT_NOTE}"
            ) from exc
        enforce_localisation_requirement(
            site_metadata=site_metadata_frame,
            field_name=field_name,
            workflow_name="signalome workflow request",
            requirement=localisation_requirement,
            error_type=WorkflowValidationError,
        )
        SignalomeWorkflowValidator._require_signalome_protein_grouping_metadata(
            site_metadata=site_metadata_frame,
            field_name=field_name,
        )
        require_exact_index_match(
            left=prediction_sites,
            right=site_metadata_frame.index,
            left_name=(
                "signalome workflow request "
                "kinase_result.prediction_result.pred_mat.index"
            ),
            right_name=f"{field_name}.index",
            error_type=WorkflowValidationError,
        )
        require_exact_index_match(
            left=score_sites,
            right=site_metadata_frame.index,
            left_name=f"{score_field_name}.index",
            right_name=f"{field_name}.index",
            error_type=WorkflowValidationError,
        )

    @staticmethod
    def _require_signalome_protein_grouping_metadata(
        *,
        site_metadata: pd.DataFrame,
        field_name: str,
    ) -> None:
        try:
            column_name = _resolve_signalome_protein_grouping_column(
                site_metadata=site_metadata,
                field_name=field_name,
            )
            enforce_required_non_empty_string_column(
                site_metadata=site_metadata,
                field_name=field_name,
                workflow_name="signalome protein grouping metadata",
                column_name=column_name,
                error_type=WorkflowValidationError,
            )
        except WorkflowValidationError as exc:
            if str(exc).startswith("Conflicting signalome protein grouping metadata"):
                raise
            raise WorkflowValidationError(
                "Missing signalome protein grouping metadata: protein_group_id "
                "(legacy alias: protein_id); "
                "signalome protein grouping metadata requirement failed: "
                f"{exc}. {SIGNALOME_PROTEIN_GROUPING_METADATA_NOTE}"
            ) from exc


def _unchecked_signalome_workflow_request_boundary(
    request: object,
) -> ValidatedSignalomeWorkflowRequest:
    """Wrap an internal interpreter-only request without rerunning validation.

    Ordinary `SignalomeWorkflow.run(...)` calls must use
    `SignalomeWorkflowValidator.run(...)`. This compatibility seam exists for
    focused interpreter tests and internal collaborators that intentionally
    exercise interpreter boundary errors before full workflow validation would
    accept the public request.
    """

    if not isinstance(request, SignalomeWorkflowRequest):
        raise WorkflowValidationError(
            "signalome workflow input must be a SignalomeWorkflowRequest"
        )
    if not isinstance(request.kinase_result, KinaseWorkflowResult):
        raise WorkflowValidationError(
            "signalome workflow request kinase_result must be KinaseWorkflowResult"
        )
    return ValidatedSignalomeWorkflowRequest(
        request=request,
        dataset_view=DatasetInternalView(request.kinase_result.dataset),
        prediction_matrix=request.kinase_result.prediction_result.pred_mat,
        authoritative_score_matrix=(
            request.kinase_result.scoring_result.authoritative_scores
        ),
    )


def _resolve_signalome_protein_grouping_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
) -> str:
    has_new = SIGNALOME_PROTEIN_GROUP_ID_COLUMN in site_metadata.columns
    has_legacy = SIGNALOME_LEGACY_PROTEIN_GROUP_ID_COLUMN in site_metadata.columns
    if has_new and has_legacy:
        new_values = _normalised_grouping_values(
            site_metadata.loc[:, SIGNALOME_PROTEIN_GROUP_ID_COLUMN]
        )
        legacy_values = _normalised_grouping_values(
            site_metadata.loc[:, SIGNALOME_LEGACY_PROTEIN_GROUP_ID_COLUMN]
        )
        mismatch_mask = new_values.ne(legacy_values)
        if bool(mismatch_mask.any()):
            mismatch_index = site_metadata.index[mismatch_mask.to_numpy()]
            mismatch_rows = [
                str(row_id) for row_id in mismatch_index.astype(str).tolist()
            ]
            preview = ", ".join(mismatch_rows[:5])
            suffix = "" if len(mismatch_rows) <= 5 else " ..."
            raise WorkflowValidationError(
                "Conflicting signalome protein grouping metadata: "
                f"{field_name}.protein_group_id and legacy alias "
                f"{field_name}.protein_id contain different values; "
                f"mismatch_site_keys=[{preview}{suffix}]. "
                f"{SIGNALOME_PROTEIN_GROUPING_METADATA_NOTE}"
            )
        return SIGNALOME_PROTEIN_GROUP_ID_COLUMN
    if has_new:
        return SIGNALOME_PROTEIN_GROUP_ID_COLUMN
    if has_legacy:
        return SIGNALOME_LEGACY_PROTEIN_GROUP_ID_COLUMN
    raise WorkflowValidationError(
        f"{field_name} is missing required columns: protein_group_id "
        "(legacy alias: protein_id). "
        f"{SIGNALOME_PROTEIN_GROUPING_METADATA_NOTE}"
    )


def _normalised_grouping_values(values: pd.Series) -> pd.Series:
    return values.fillna("").astype(str).str.strip()


def _require_reference_context_compatibility(
    kinase_result: KinaseWorkflowResult,
    *,
    policy: ReferenceContextCompatibilityPolicy,
) -> None:
    dataset_context = kinase_result.dataset.reference_context
    result_provenance = kinase_result.provenance
    allow_unknown = (
        policy is ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
    )
    validate_reference_context_compatibility(
        dataset_context,
        None if result_provenance is None else result_provenance.reference_context,
        operation="signalome workflow request dataset/upstream kinase result",
        allow_unknown=allow_unknown,
        error_type=WorkflowValidationError,
    )
    validate_reference_context_compatibility(
        dataset_context,
        _reference_bundle_context(kinase_result),
        operation="signalome workflow request dataset/upstream kinase reference",
        allow_unknown=allow_unknown,
        error_type=WorkflowValidationError,
    )


def _reference_bundle_context(
    kinase_result: KinaseWorkflowResult,
) -> ReferenceContext | None:
    provenance = kinase_result.references.provenance
    if provenance is None:
        return None
    return provenance.reference_context

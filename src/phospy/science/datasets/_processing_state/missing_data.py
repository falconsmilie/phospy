"""Missing-data diagnostics models and parsing contracts."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field

from phospy.errors.input import PhosPyInputError
from phospy.science.configs.preprocessing import (
    DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICIES,
    DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_VERSION,
)
from phospy.science.datasets.preprocessing.imputation_scale_policy import (
    IMPUTATION_INPUT_SCALE_SOURCE_CALLER_SELECTED,
    IMPUTATION_INPUT_SCALE_SOURCE_METHOD_REQUIRED,
    IMPUTATION_OPERATION_ORDERS,
)
from phospy.science.datasets.preprocessing.policy_models import MissingDataPolicy

from .json_contracts import (
    MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
    V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS,
    FrozenJsonMapping,
    JsonValue,
    require_frozen_json_mapping,
    require_int,
    require_mapping,
    require_optional_bool,
    require_optional_frozen_json_mapping,
    require_optional_frozen_string_to_float_mapping,
    require_optional_int,
    require_optional_non_negative_int,
    require_optional_str,
    require_optional_string_tuple,
    require_required_non_negative_int,
    require_required_str,
    require_required_string_tuple,
    require_string_keys,
    set_optional_payload_value,
    thaw_frozen_json_mapping,
)


class MissingDataDiagnostics(Mapping[str, JsonValue]):
    """Typed diagnostics contract for missing-data preprocessing state."""

    diagnostics_schema_version: int

    @classmethod
    def from_payload(
        cls,
        payload: object,
        *,
        field_name: str,
    ) -> MissingDataDiagnostics:
        mapping = require_mapping(payload, field_name=field_name)
        return MissingDataDiagnosticsV1.from_mapping(mapping, field_name=field_name)

    def to_payload(self) -> dict[str, JsonValue]:
        """Return normalized diagnostics payload suitable for bundle JSON."""

        raise NotImplementedError


@dataclass(frozen=True, slots=True, eq=False)
class MissingDataDiagnosticsV1(MissingDataDiagnostics):
    """Versioned missing-data diagnostics payload (schema v1)."""

    missing_data_policy: str
    input_missing_cell_count: int
    output_missing_cell_count: int
    imputed_cell_count: int
    affected_row_count: int
    affected_column_count: int
    affected_row_ids: tuple[str, ...]
    affected_column_ids: tuple[str, ...]
    imputed_row_ids: tuple[str, ...]
    imputed_column_ids: tuple[str, ...]
    dropped_row_ids: tuple[str, ...]
    method_parameters: Mapping[str, object]
    stage_order: tuple[str, ...]
    missingness_mask_hash: str
    rows_not_imputable: tuple[str, ...]
    row_medians_used: Mapping[str, object] = field(default_factory=dict)
    imputed_row_count: int | None = None
    imputed_column_count: int | None = None
    dropped_row_count: int | None = None
    imputation_mask_hash: str | None = None
    imputation_method_id: str | None = None
    imputation_method_family: str | None = None
    random_seed: int | None = None
    matrix_scale_requirement: str | None = None
    imputation_input_scale: str | None = None
    imputation_input_scale_source: str | None = None
    imputation_operation_order: str | None = None
    left_censored_assumption: bool | None = None
    per_column_distribution_parameters: Mapping[str, object] | None = None
    dropped_rows_above_max_missing_fraction: tuple[str, ...] | None = None
    neighbour_count: int | None = None
    distance_metric: str | None = None
    knn_no_overlap_policy: str | None = None
    knn_no_overlap_policy_version: int | None = None
    knn_nearest_neighbour_imputed_cell_count: int | None = None
    knn_nearest_neighbour_imputed_row_ids: tuple[str, ...] | None = None
    knn_nearest_neighbour_imputed_column_ids: tuple[str, ...] | None = None
    knn_column_mean_fallback_imputed_cell_count: int | None = None
    knn_column_mean_fallback_row_ids: tuple[str, ...] | None = None
    knn_column_mean_fallback_column_ids: tuple[str, ...] | None = None
    knn_nearest_neighbour_imputation_mask_hash: str | None = None
    knn_column_mean_fallback_imputation_mask_hash: str | None = None
    knn_fully_column_mean_fallback_row_ids: tuple[str, ...] | None = None
    diagnostic_caveat_codes: tuple[str, ...] | None = None
    diagnostics_schema_version: int = MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1

    @classmethod
    def from_mapping(
        cls,
        payload: Mapping[str, object],
        *,
        field_name: str,
    ) -> MissingDataDiagnosticsV1:
        if "diagnostics_schema_version" not in payload:
            raise PhosPyInputError(
                f"{field_name}.diagnostics_schema_version is required"
            )
        return cls._from_versioned_payload(payload, field_name=field_name)

    @classmethod
    def _from_versioned_payload(
        cls,
        payload: Mapping[str, object],
        *,
        field_name: str,
    ) -> MissingDataDiagnosticsV1:
        require_string_keys(payload, field_name=field_name)
        version = require_int(
            payload.get("diagnostics_schema_version"),
            field_name=f"{field_name}.diagnostics_schema_version",
        )
        if version != MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1:
            raise PhosPyInputError(
                f"{field_name}.diagnostics_schema_version={version!r} is unsupported; "
                f"expected {MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1}"
            )
        unknown_fields = sorted(
            key
            for key in payload
            if key not in V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS
        )
        if unknown_fields:
            raise PhosPyInputError(
                f"{field_name} contains unsupported field(s): "
                + ", ".join(unknown_fields)
            )
        return cls(
            diagnostics_schema_version=version,
            missing_data_policy=require_required_str(
                payload.get("missing_data_policy"),
                field_name=f"{field_name}.missing_data_policy",
            ),
            imputation_method_id=require_optional_str(
                payload.get("imputation_method_id"),
                field_name=f"{field_name}.imputation_method_id",
            ),
            imputation_method_family=require_optional_str(
                payload.get("imputation_method_family"),
                field_name=f"{field_name}.imputation_method_family",
            ),
            input_missing_cell_count=require_required_non_negative_int(
                payload.get("input_missing_cell_count"),
                field_name=f"{field_name}.input_missing_cell_count",
            ),
            output_missing_cell_count=require_required_non_negative_int(
                payload.get("output_missing_cell_count"),
                field_name=f"{field_name}.output_missing_cell_count",
            ),
            imputed_cell_count=require_required_non_negative_int(
                payload.get("imputed_cell_count"),
                field_name=f"{field_name}.imputed_cell_count",
            ),
            affected_row_count=require_required_non_negative_int(
                payload.get("affected_row_count"),
                field_name=f"{field_name}.affected_row_count",
            ),
            affected_column_count=require_required_non_negative_int(
                payload.get("affected_column_count"),
                field_name=f"{field_name}.affected_column_count",
            ),
            affected_row_ids=require_required_string_tuple(
                payload.get("affected_row_ids"),
                field_name=f"{field_name}.affected_row_ids",
            ),
            affected_column_ids=require_required_string_tuple(
                payload.get("affected_column_ids"),
                field_name=f"{field_name}.affected_column_ids",
            ),
            imputed_row_ids=require_required_string_tuple(
                payload.get("imputed_row_ids"),
                field_name=f"{field_name}.imputed_row_ids",
            ),
            imputed_column_ids=require_required_string_tuple(
                payload.get("imputed_column_ids"),
                field_name=f"{field_name}.imputed_column_ids",
            ),
            dropped_row_ids=require_required_string_tuple(
                payload.get("dropped_row_ids"),
                field_name=f"{field_name}.dropped_row_ids",
            ),
            imputed_row_count=require_optional_non_negative_int(
                payload.get("imputed_row_count"),
                field_name=f"{field_name}.imputed_row_count",
            ),
            imputed_column_count=require_optional_non_negative_int(
                payload.get("imputed_column_count"),
                field_name=f"{field_name}.imputed_column_count",
            ),
            dropped_row_count=require_optional_non_negative_int(
                payload.get("dropped_row_count"),
                field_name=f"{field_name}.dropped_row_count",
            ),
            random_seed=require_optional_int(
                payload.get("random_seed"),
                field_name=f"{field_name}.random_seed",
            ),
            method_parameters=require_frozen_json_mapping(
                payload.get("method_parameters"),
                field_name=f"{field_name}.method_parameters",
            ),
            matrix_scale_requirement=require_optional_str(
                payload.get("matrix_scale_requirement"),
                field_name=f"{field_name}.matrix_scale_requirement",
            ),
            imputation_input_scale=require_optional_str(
                payload.get("imputation_input_scale"),
                field_name=f"{field_name}.imputation_input_scale",
            ),
            imputation_input_scale_source=require_optional_str(
                payload.get("imputation_input_scale_source"),
                field_name=f"{field_name}.imputation_input_scale_source",
            ),
            imputation_operation_order=require_optional_str(
                payload.get("imputation_operation_order"),
                field_name=f"{field_name}.imputation_operation_order",
            ),
            stage_order=require_required_string_tuple(
                payload.get("stage_order"),
                field_name=f"{field_name}.stage_order",
            ),
            missingness_mask_hash=require_required_str(
                payload.get("missingness_mask_hash"),
                field_name=f"{field_name}.missingness_mask_hash",
            ),
            imputation_mask_hash=require_optional_str(
                payload.get("imputation_mask_hash"),
                field_name=f"{field_name}.imputation_mask_hash",
            ),
            left_censored_assumption=require_optional_bool(
                payload.get("left_censored_assumption"),
                field_name=f"{field_name}.left_censored_assumption",
            ),
            rows_not_imputable=require_required_string_tuple(
                payload.get("rows_not_imputable"),
                field_name=f"{field_name}.rows_not_imputable",
            ),
            row_medians_used=require_optional_frozen_string_to_float_mapping(
                payload.get("row_medians_used"),
                field_name=f"{field_name}.row_medians_used",
            )
            or FrozenJsonMapping(),
            per_column_distribution_parameters=require_optional_frozen_json_mapping(
                payload.get("per_column_distribution_parameters"),
                field_name=f"{field_name}.per_column_distribution_parameters",
            ),
            dropped_rows_above_max_missing_fraction=require_optional_string_tuple(
                payload.get("dropped_rows_above_max_missing_fraction"),
                field_name=f"{field_name}.dropped_rows_above_max_missing_fraction",
            ),
            neighbour_count=require_optional_int(
                payload.get("neighbour_count"),
                field_name=f"{field_name}.neighbour_count",
            ),
            distance_metric=require_optional_str(
                payload.get("distance_metric"),
                field_name=f"{field_name}.distance_metric",
            ),
            knn_no_overlap_policy=require_optional_str(
                payload.get("knn_no_overlap_policy"),
                field_name=f"{field_name}.knn_no_overlap_policy",
            ),
            knn_no_overlap_policy_version=require_optional_non_negative_int(
                payload.get("knn_no_overlap_policy_version"),
                field_name=f"{field_name}.knn_no_overlap_policy_version",
            ),
            knn_nearest_neighbour_imputed_cell_count=(
                require_optional_non_negative_int(
                    payload.get("knn_nearest_neighbour_imputed_cell_count"),
                    field_name=(
                        f"{field_name}.knn_nearest_neighbour_imputed_cell_count"
                    ),
                )
            ),
            knn_nearest_neighbour_imputed_row_ids=require_optional_string_tuple(
                payload.get("knn_nearest_neighbour_imputed_row_ids"),
                field_name=f"{field_name}.knn_nearest_neighbour_imputed_row_ids",
            ),
            knn_nearest_neighbour_imputed_column_ids=require_optional_string_tuple(
                payload.get("knn_nearest_neighbour_imputed_column_ids"),
                field_name=(f"{field_name}.knn_nearest_neighbour_imputed_column_ids"),
            ),
            knn_column_mean_fallback_imputed_cell_count=(
                require_optional_non_negative_int(
                    payload.get("knn_column_mean_fallback_imputed_cell_count"),
                    field_name=(
                        f"{field_name}.knn_column_mean_fallback_imputed_cell_count"
                    ),
                )
            ),
            knn_column_mean_fallback_row_ids=require_optional_string_tuple(
                payload.get("knn_column_mean_fallback_row_ids"),
                field_name=f"{field_name}.knn_column_mean_fallback_row_ids",
            ),
            knn_column_mean_fallback_column_ids=require_optional_string_tuple(
                payload.get("knn_column_mean_fallback_column_ids"),
                field_name=f"{field_name}.knn_column_mean_fallback_column_ids",
            ),
            knn_nearest_neighbour_imputation_mask_hash=require_optional_str(
                payload.get("knn_nearest_neighbour_imputation_mask_hash"),
                field_name=f"{field_name}.knn_nearest_neighbour_imputation_mask_hash",
            ),
            knn_column_mean_fallback_imputation_mask_hash=require_optional_str(
                payload.get("knn_column_mean_fallback_imputation_mask_hash"),
                field_name=(
                    f"{field_name}.knn_column_mean_fallback_imputation_mask_hash"
                ),
            ),
            knn_fully_column_mean_fallback_row_ids=require_optional_string_tuple(
                payload.get("knn_fully_column_mean_fallback_row_ids"),
                field_name=f"{field_name}.knn_fully_column_mean_fallback_row_ids",
            ),
            diagnostic_caveat_codes=require_optional_string_tuple(
                payload.get("diagnostic_caveat_codes"),
                field_name=f"{field_name}.diagnostic_caveat_codes",
            ),
        )

    def __post_init__(self) -> None:
        if (
            self.diagnostics_schema_version
            != MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1
        ):
            raise PhosPyInputError(
                "dataset processing state missing_data diagnostics schema version "
                f"must be {MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1}"
            )
        missing_data_policy = require_required_str(
            self.missing_data_policy,
            field_name="dataset processing state missing_data.diagnostics.missing_data_policy",
        )
        missing_data_policy = MissingDataPolicy.parse(
            missing_data_policy,
            field_name="dataset processing state missing_data.diagnostics.missing_data_policy",
        ).value
        imputation_method_id = require_optional_str(
            self.imputation_method_id,
            field_name="dataset processing state missing_data.diagnostics.imputation_method_id",
        )
        imputation_method_family = require_optional_str(
            self.imputation_method_family,
            field_name="dataset processing state missing_data.diagnostics.imputation_method_family",
        )
        input_missing_cell_count = require_required_non_negative_int(
            self.input_missing_cell_count,
            field_name="dataset processing state missing_data.diagnostics.input_missing_cell_count",
        )
        output_missing_cell_count = require_required_non_negative_int(
            self.output_missing_cell_count,
            field_name="dataset processing state missing_data.diagnostics.output_missing_cell_count",
        )
        imputed_cell_count = require_required_non_negative_int(
            self.imputed_cell_count,
            field_name="dataset processing state missing_data.diagnostics.imputed_cell_count",
        )
        affected_row_count = require_required_non_negative_int(
            self.affected_row_count,
            field_name="dataset processing state missing_data.diagnostics.affected_row_count",
        )
        affected_column_count = require_required_non_negative_int(
            self.affected_column_count,
            field_name="dataset processing state missing_data.diagnostics.affected_column_count",
        )
        affected_row_ids = require_required_string_tuple(
            self.affected_row_ids,
            field_name="dataset processing state missing_data.diagnostics.affected_row_ids",
        )
        affected_column_ids = require_required_string_tuple(
            self.affected_column_ids,
            field_name="dataset processing state missing_data.diagnostics.affected_column_ids",
        )
        imputed_row_ids = require_required_string_tuple(
            self.imputed_row_ids,
            field_name="dataset processing state missing_data.diagnostics.imputed_row_ids",
        )
        imputed_column_ids = require_required_string_tuple(
            self.imputed_column_ids,
            field_name="dataset processing state missing_data.diagnostics.imputed_column_ids",
        )
        dropped_row_ids = require_required_string_tuple(
            self.dropped_row_ids,
            field_name="dataset processing state missing_data.diagnostics.dropped_row_ids",
        )
        imputed_row_count = require_optional_non_negative_int(
            self.imputed_row_count,
            field_name="dataset processing state missing_data.diagnostics.imputed_row_count",
        )
        imputed_column_count = require_optional_non_negative_int(
            self.imputed_column_count,
            field_name=(
                "dataset processing state missing_data.diagnostics.imputed_column_count"
            ),
        )
        dropped_row_count = require_optional_non_negative_int(
            self.dropped_row_count,
            field_name="dataset processing state missing_data.diagnostics.dropped_row_count",
        )
        random_seed = require_optional_int(
            self.random_seed,
            field_name="dataset processing state missing_data.diagnostics.random_seed",
        )
        method_parameters = require_frozen_json_mapping(
            self.method_parameters,
            field_name="dataset processing state missing_data.diagnostics.method_parameters",
        )
        matrix_scale_requirement = require_optional_str(
            self.matrix_scale_requirement,
            field_name="dataset processing state missing_data.diagnostics.matrix_scale_requirement",
        )
        imputation_input_scale = _require_optional_imputation_input_scale(
            self.imputation_input_scale
        )
        imputation_input_scale_source = _require_optional_imputation_input_scale_source(
            self.imputation_input_scale_source
        )
        imputation_operation_order = _require_optional_imputation_operation_order(
            self.imputation_operation_order
        )
        stage_order = require_required_string_tuple(
            self.stage_order,
            field_name="dataset processing state missing_data.diagnostics.stage_order",
        )
        missingness_mask_hash = require_required_str(
            self.missingness_mask_hash,
            field_name="dataset processing state missing_data.diagnostics.missingness_mask_hash",
        )
        imputation_mask_hash = require_optional_str(
            self.imputation_mask_hash,
            field_name="dataset processing state missing_data.diagnostics.imputation_mask_hash",
        )
        left_censored_assumption = require_optional_bool(
            self.left_censored_assumption,
            field_name="dataset processing state missing_data.diagnostics.left_censored_assumption",
        )
        rows_not_imputable = require_required_string_tuple(
            self.rows_not_imputable,
            field_name="dataset processing state missing_data.diagnostics.rows_not_imputable",
        )
        row_medians_used = (
            require_optional_frozen_string_to_float_mapping(
                self.row_medians_used,
                field_name="dataset processing state missing_data.diagnostics.row_medians_used",
            )
            or FrozenJsonMapping()
        )
        if imputation_method_id != "row_median":
            row_medians_used = FrozenJsonMapping()
        per_column_distribution_parameters = require_optional_frozen_json_mapping(
            self.per_column_distribution_parameters,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "per_column_distribution_parameters"
            ),
        )
        dropped_rows_above_max_missing_fraction = require_optional_string_tuple(
            self.dropped_rows_above_max_missing_fraction,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "dropped_rows_above_max_missing_fraction"
            ),
        )
        neighbour_count = require_optional_int(
            self.neighbour_count,
            field_name="dataset processing state missing_data.diagnostics.neighbour_count",
        )
        if neighbour_count is not None and neighbour_count < 1:
            raise PhosPyInputError(
                "dataset processing state missing_data.diagnostics."
                "neighbour_count must be >= 1"
            )
        if imputed_row_count is not None and imputed_row_count != int(
            len(imputed_row_ids)
        ):
            raise PhosPyInputError(
                "dataset processing state missing_data.diagnostics.imputed_row_count "
                "must match len(imputed_row_ids)"
            )
        if imputed_column_count is not None and imputed_column_count != int(
            len(imputed_column_ids)
        ):
            raise PhosPyInputError(
                "dataset processing state missing_data.diagnostics."
                "imputed_column_count must match len(imputed_column_ids)"
            )
        if dropped_row_count is not None and dropped_row_count != int(
            len(dropped_row_ids)
        ):
            raise PhosPyInputError(
                "dataset processing state missing_data.diagnostics.dropped_row_count "
                "must match len(dropped_row_ids)"
            )
        if (
            imputation_method_id in {"row_median", "knn", "minprob"}
            and imputation_mask_hash is None
        ):
            raise PhosPyInputError(
                "dataset processing state missing_data.diagnostics."
                "imputation_mask_hash is required for imputation methods"
            )
        distance_metric = require_optional_str(
            self.distance_metric,
            field_name="dataset processing state missing_data.diagnostics.distance_metric",
        )
        knn_no_overlap_policy = require_optional_str(
            self.knn_no_overlap_policy,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "knn_no_overlap_policy"
            ),
        )
        if (
            knn_no_overlap_policy is not None
            and knn_no_overlap_policy
            not in DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICIES
        ):
            supported = ", ".join(sorted(DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICIES))
            raise PhosPyInputError(
                "dataset processing state missing_data.diagnostics."
                f"knn_no_overlap_policy must be one of: {supported}"
            )
        knn_no_overlap_policy_version = require_optional_non_negative_int(
            self.knn_no_overlap_policy_version,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "knn_no_overlap_policy_version"
            ),
        )
        if (
            knn_no_overlap_policy_version is not None
            and knn_no_overlap_policy_version
            != DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_VERSION
        ):
            raise PhosPyInputError(
                "dataset processing state missing_data.diagnostics."
                "knn_no_overlap_policy_version is unsupported; expected "
                f"{DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_VERSION}"
            )
        knn_nearest_neighbour_imputed_cell_count = require_optional_non_negative_int(
            self.knn_nearest_neighbour_imputed_cell_count,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "knn_nearest_neighbour_imputed_cell_count"
            ),
        )
        knn_column_mean_fallback_imputed_cell_count = require_optional_non_negative_int(
            self.knn_column_mean_fallback_imputed_cell_count,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "knn_column_mean_fallback_imputed_cell_count"
            ),
        )
        if (
            knn_nearest_neighbour_imputed_cell_count is not None
            and knn_column_mean_fallback_imputed_cell_count is not None
            and (
                knn_nearest_neighbour_imputed_cell_count
                + knn_column_mean_fallback_imputed_cell_count
                != imputed_cell_count
            )
        ):
            raise PhosPyInputError(
                "dataset processing state missing_data.diagnostics KNN mechanism "
                "cell counts must sum to imputed_cell_count"
            )
        knn_nearest_neighbour_imputed_row_ids = require_optional_string_tuple(
            self.knn_nearest_neighbour_imputed_row_ids,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "knn_nearest_neighbour_imputed_row_ids"
            ),
        )
        knn_nearest_neighbour_imputed_column_ids = require_optional_string_tuple(
            self.knn_nearest_neighbour_imputed_column_ids,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "knn_nearest_neighbour_imputed_column_ids"
            ),
        )
        knn_column_mean_fallback_row_ids = require_optional_string_tuple(
            self.knn_column_mean_fallback_row_ids,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "knn_column_mean_fallback_row_ids"
            ),
        )
        knn_column_mean_fallback_column_ids = require_optional_string_tuple(
            self.knn_column_mean_fallback_column_ids,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "knn_column_mean_fallback_column_ids"
            ),
        )
        knn_nearest_neighbour_imputation_mask_hash = require_optional_str(
            self.knn_nearest_neighbour_imputation_mask_hash,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "knn_nearest_neighbour_imputation_mask_hash"
            ),
        )
        knn_column_mean_fallback_imputation_mask_hash = require_optional_str(
            self.knn_column_mean_fallback_imputation_mask_hash,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "knn_column_mean_fallback_imputation_mask_hash"
            ),
        )
        knn_fully_column_mean_fallback_row_ids = require_optional_string_tuple(
            self.knn_fully_column_mean_fallback_row_ids,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "knn_fully_column_mean_fallback_row_ids"
            ),
        )
        diagnostic_caveat_codes = require_optional_string_tuple(
            self.diagnostic_caveat_codes,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "diagnostic_caveat_codes"
            ),
        )
        object.__setattr__(self, "missing_data_policy", missing_data_policy)
        object.__setattr__(self, "imputation_method_id", imputation_method_id)
        object.__setattr__(self, "imputation_method_family", imputation_method_family)
        object.__setattr__(self, "input_missing_cell_count", input_missing_cell_count)
        object.__setattr__(self, "output_missing_cell_count", output_missing_cell_count)
        object.__setattr__(self, "imputed_cell_count", imputed_cell_count)
        object.__setattr__(self, "affected_row_count", affected_row_count)
        object.__setattr__(self, "affected_column_count", affected_column_count)
        object.__setattr__(self, "affected_row_ids", affected_row_ids)
        object.__setattr__(self, "affected_column_ids", affected_column_ids)
        object.__setattr__(self, "imputed_row_ids", imputed_row_ids)
        object.__setattr__(self, "imputed_column_ids", imputed_column_ids)
        object.__setattr__(self, "dropped_row_ids", dropped_row_ids)
        object.__setattr__(self, "imputed_row_count", imputed_row_count)
        object.__setattr__(self, "imputed_column_count", imputed_column_count)
        object.__setattr__(self, "dropped_row_count", dropped_row_count)
        object.__setattr__(self, "random_seed", random_seed)
        object.__setattr__(self, "method_parameters", method_parameters)
        object.__setattr__(self, "matrix_scale_requirement", matrix_scale_requirement)
        object.__setattr__(self, "imputation_input_scale", imputation_input_scale)
        object.__setattr__(
            self,
            "imputation_input_scale_source",
            imputation_input_scale_source,
        )
        object.__setattr__(
            self,
            "imputation_operation_order",
            imputation_operation_order,
        )
        object.__setattr__(self, "stage_order", stage_order)
        object.__setattr__(self, "missingness_mask_hash", missingness_mask_hash)
        object.__setattr__(self, "imputation_mask_hash", imputation_mask_hash)
        object.__setattr__(self, "left_censored_assumption", left_censored_assumption)
        object.__setattr__(self, "rows_not_imputable", rows_not_imputable)
        object.__setattr__(self, "row_medians_used", row_medians_used)
        object.__setattr__(
            self,
            "per_column_distribution_parameters",
            per_column_distribution_parameters,
        )
        object.__setattr__(
            self,
            "dropped_rows_above_max_missing_fraction",
            dropped_rows_above_max_missing_fraction,
        )
        object.__setattr__(self, "neighbour_count", neighbour_count)
        object.__setattr__(self, "distance_metric", distance_metric)
        object.__setattr__(self, "knn_no_overlap_policy", knn_no_overlap_policy)
        object.__setattr__(
            self,
            "knn_no_overlap_policy_version",
            knn_no_overlap_policy_version,
        )
        object.__setattr__(
            self,
            "knn_nearest_neighbour_imputed_cell_count",
            knn_nearest_neighbour_imputed_cell_count,
        )
        object.__setattr__(
            self,
            "knn_nearest_neighbour_imputed_row_ids",
            knn_nearest_neighbour_imputed_row_ids,
        )
        object.__setattr__(
            self,
            "knn_nearest_neighbour_imputed_column_ids",
            knn_nearest_neighbour_imputed_column_ids,
        )
        object.__setattr__(
            self,
            "knn_column_mean_fallback_imputed_cell_count",
            knn_column_mean_fallback_imputed_cell_count,
        )
        object.__setattr__(
            self,
            "knn_column_mean_fallback_row_ids",
            knn_column_mean_fallback_row_ids,
        )
        object.__setattr__(
            self,
            "knn_column_mean_fallback_column_ids",
            knn_column_mean_fallback_column_ids,
        )
        object.__setattr__(
            self,
            "knn_nearest_neighbour_imputation_mask_hash",
            knn_nearest_neighbour_imputation_mask_hash,
        )
        object.__setattr__(
            self,
            "knn_column_mean_fallback_imputation_mask_hash",
            knn_column_mean_fallback_imputation_mask_hash,
        )
        object.__setattr__(
            self,
            "knn_fully_column_mean_fallback_row_ids",
            knn_fully_column_mean_fallback_row_ids,
        )
        object.__setattr__(self, "diagnostic_caveat_codes", diagnostic_caveat_codes)

    def __getitem__(self, key: str) -> JsonValue:
        return self.to_payload()[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self.to_payload())

    def __len__(self) -> int:
        return len(self.to_payload())

    def to_payload(self) -> dict[str, JsonValue]:
        payload: dict[str, JsonValue] = {
            "diagnostics_schema_version": int(self.diagnostics_schema_version),
            "missing_data_policy": self.missing_data_policy,
            "input_missing_cell_count": int(self.input_missing_cell_count),
            "output_missing_cell_count": int(self.output_missing_cell_count),
            "imputed_cell_count": int(self.imputed_cell_count),
            "affected_row_count": int(self.affected_row_count),
            "affected_column_count": int(self.affected_column_count),
            "affected_row_ids": list(self.affected_row_ids),
            "affected_column_ids": list(self.affected_column_ids),
            "imputed_row_ids": list(self.imputed_row_ids),
            "imputed_column_ids": list(self.imputed_column_ids),
            "dropped_row_ids": list(self.dropped_row_ids),
            "method_parameters": thaw_frozen_json_mapping(
                self.method_parameters,
                field_name="dataset processing state missing_data.diagnostics.method_parameters",
            ),
            "stage_order": list(self.stage_order),
            "missingness_mask_hash": self.missingness_mask_hash,
            "rows_not_imputable": list(self.rows_not_imputable),
            "row_medians_used": thaw_frozen_json_mapping(
                self.row_medians_used,
                field_name="dataset processing state missing_data.diagnostics.row_medians_used",
            ),
        }
        set_optional_payload_value(payload, "imputed_row_count", self.imputed_row_count)
        set_optional_payload_value(
            payload, "imputed_column_count", self.imputed_column_count
        )
        set_optional_payload_value(payload, "dropped_row_count", self.dropped_row_count)
        set_optional_payload_value(
            payload, "imputation_method_id", self.imputation_method_id
        )
        set_optional_payload_value(
            payload, "imputation_method_family", self.imputation_method_family
        )
        set_optional_payload_value(payload, "random_seed", self.random_seed)
        set_optional_payload_value(
            payload, "matrix_scale_requirement", self.matrix_scale_requirement
        )
        set_optional_payload_value(
            payload, "imputation_input_scale", self.imputation_input_scale
        )
        set_optional_payload_value(
            payload,
            "imputation_input_scale_source",
            self.imputation_input_scale_source,
        )
        set_optional_payload_value(
            payload,
            "imputation_operation_order",
            self.imputation_operation_order,
        )
        set_optional_payload_value(
            payload, "left_censored_assumption", self.left_censored_assumption
        )
        set_optional_payload_value(
            payload,
            "imputation_mask_hash",
            self.imputation_mask_hash,
        )
        if self.per_column_distribution_parameters is not None:
            payload["per_column_distribution_parameters"] = thaw_frozen_json_mapping(
                self.per_column_distribution_parameters,
                field_name=(
                    "dataset processing state missing_data.diagnostics."
                    "per_column_distribution_parameters"
                ),
            )
        if self.dropped_rows_above_max_missing_fraction is not None:
            payload["dropped_rows_above_max_missing_fraction"] = list(
                self.dropped_rows_above_max_missing_fraction
            )
        set_optional_payload_value(payload, "neighbour_count", self.neighbour_count)
        set_optional_payload_value(payload, "distance_metric", self.distance_metric)
        set_optional_payload_value(
            payload,
            "knn_no_overlap_policy",
            self.knn_no_overlap_policy,
        )
        set_optional_payload_value(
            payload,
            "knn_no_overlap_policy_version",
            self.knn_no_overlap_policy_version,
        )
        set_optional_payload_value(
            payload,
            "knn_nearest_neighbour_imputed_cell_count",
            self.knn_nearest_neighbour_imputed_cell_count,
        )
        if self.knn_nearest_neighbour_imputed_row_ids is not None:
            payload["knn_nearest_neighbour_imputed_row_ids"] = list(
                self.knn_nearest_neighbour_imputed_row_ids
            )
        if self.knn_nearest_neighbour_imputed_column_ids is not None:
            payload["knn_nearest_neighbour_imputed_column_ids"] = list(
                self.knn_nearest_neighbour_imputed_column_ids
            )
        set_optional_payload_value(
            payload,
            "knn_column_mean_fallback_imputed_cell_count",
            self.knn_column_mean_fallback_imputed_cell_count,
        )
        if self.knn_column_mean_fallback_row_ids is not None:
            payload["knn_column_mean_fallback_row_ids"] = list(
                self.knn_column_mean_fallback_row_ids
            )
        if self.knn_column_mean_fallback_column_ids is not None:
            payload["knn_column_mean_fallback_column_ids"] = list(
                self.knn_column_mean_fallback_column_ids
            )
        set_optional_payload_value(
            payload,
            "knn_nearest_neighbour_imputation_mask_hash",
            self.knn_nearest_neighbour_imputation_mask_hash,
        )
        set_optional_payload_value(
            payload,
            "knn_column_mean_fallback_imputation_mask_hash",
            self.knn_column_mean_fallback_imputation_mask_hash,
        )
        if self.knn_fully_column_mean_fallback_row_ids is not None:
            payload["knn_fully_column_mean_fallback_row_ids"] = list(
                self.knn_fully_column_mean_fallback_row_ids
            )
        if self.diagnostic_caveat_codes is not None:
            payload["diagnostic_caveat_codes"] = list(self.diagnostic_caveat_codes)
        return payload


def _require_optional_imputation_input_scale(value: object) -> str | None:
    parsed = require_optional_str(
        value,
        field_name=(
            "dataset processing state missing_data.diagnostics.imputation_input_scale"
        ),
    )
    if parsed is None or parsed in {"linear", "log2"}:
        return parsed
    raise PhosPyInputError(
        "dataset processing state missing_data.diagnostics."
        "imputation_input_scale must be one of: linear, log2"
    )


def _require_optional_imputation_input_scale_source(value: object) -> str | None:
    parsed = require_optional_str(
        value,
        field_name=(
            "dataset processing state missing_data.diagnostics."
            "imputation_input_scale_source"
        ),
    )
    if parsed is None or parsed in {
        IMPUTATION_INPUT_SCALE_SOURCE_CALLER_SELECTED,
        IMPUTATION_INPUT_SCALE_SOURCE_METHOD_REQUIRED,
    }:
        return parsed
    raise PhosPyInputError(
        "dataset processing state missing_data.diagnostics."
        "imputation_input_scale_source must be one of: caller_selected, "
        "method_required"
    )


def _require_optional_imputation_operation_order(value: object) -> str | None:
    parsed = require_optional_str(
        value,
        field_name=(
            "dataset processing state missing_data.diagnostics."
            "imputation_operation_order"
        ),
    )
    if parsed is None or parsed in IMPUTATION_OPERATION_ORDERS:
        return parsed
    supported = ", ".join(sorted(IMPUTATION_OPERATION_ORDERS))
    raise PhosPyInputError(
        "dataset processing state missing_data.diagnostics."
        f"imputation_operation_order must be one of: {supported}"
    )

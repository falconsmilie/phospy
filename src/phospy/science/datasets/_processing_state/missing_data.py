"""Missing-data diagnostics models and parsing contracts."""

from __future__ import annotations

import math
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field, fields
from typing import ClassVar, cast

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
    MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V2,
    V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS,
    V2_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS,
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
    require_required_label_tuple,
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
    missing_data_policy: str

    @classmethod
    def from_payload(
        cls,
        payload: object,
        *,
        field_name: str,
    ) -> MissingDataDiagnostics:
        mapping = require_mapping(payload, field_name=field_name)
        version = require_int(
            mapping.get("diagnostics_schema_version"),
            field_name=f"{field_name}.diagnostics_schema_version",
        )
        if version == MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1:
            return MissingDataDiagnosticsV1.from_mapping(mapping, field_name=field_name)
        if version == MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V2:
            return MissingDataDiagnosticsV2.from_mapping(mapping, field_name=field_name)
        raise PhosPyInputError(
            f"{field_name}.diagnostics_schema_version={version!r} is unsupported; "
            "expected one of: 1, 2"
        )

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
    _expected_schema_version: ClassVar[int] = MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1

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
            affected_column_ids=require_required_label_tuple(
                payload.get("affected_column_ids"),
                field_name=f"{field_name}.affected_column_ids",
            ),
            imputed_row_ids=require_required_string_tuple(
                payload.get("imputed_row_ids"),
                field_name=f"{field_name}.imputed_row_ids",
            ),
            imputed_column_ids=require_required_label_tuple(
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
        if self.diagnostics_schema_version != self._expected_schema_version:
            raise PhosPyInputError(
                "dataset processing state missing_data diagnostics schema version "
                f"must be {self._expected_schema_version}"
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
        affected_column_ids = require_required_label_tuple(
            self.affected_column_ids,
            field_name="dataset processing state missing_data.diagnostics.affected_column_ids",
        )
        imputed_row_ids = require_required_string_tuple(
            self.imputed_row_ids,
            field_name="dataset processing state missing_data.diagnostics.imputed_row_ids",
        )
        imputed_column_ids = require_required_label_tuple(
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


_GROUP_AWARE_ROUTE_CATEGORIES = (
    "partial_observation_knn",
    "asymmetric_absence_minprob",
    "unsupported_partial",
    "unsupported_absence",
)
_GROUP_AWARE_KNOWN_FIELDS = frozenset(
    {
        "group_column",
        "observed_group_sizes",
        "min_partial_observed_fraction",
        "min_reference_observed_fraction",
        "retained_row_count",
        "dropped_unsupported_row_count",
        "knn_routed_cell_count",
        "minprob_routed_cell_count",
        "knn_imputed_cell_count",
        "minprob_imputed_cell_count",
        "knn_target_mask_hash",
        "minprob_target_mask_hash",
        "knn_imputation_mask_hash",
        "minprob_imputation_mask_hash",
        "unsupported_partial_group_count",
        "unsupported_absence_group_count",
        "unsupported_partial_row_ids",
        "unsupported_absence_row_ids",
        "routed_rows",
        "rejected_rows",
        "minprob_left_censored_assumption",
        "route_categories",
        "mechanism_input",
    }
)
_REJECTED_ROW_KNOWN_FIELDS = frozenset(
    {
        "row_id",
        "unsupported_partial_groups",
        "unsupported_absence_groups",
        "observed_finite_count_by_group",
        "observed_fraction_by_group",
    }
)
_ROUTED_ROW_KNOWN_FIELDS = frozenset(
    {
        "row_id",
        "knn_imputed_columns",
        "minprob_imputed_columns",
        "knn_group_labels",
        "minprob_group_labels",
    }
)


@dataclass(frozen=True, slots=True)
class GroupAwareRoutedRowRecord:
    """Mechanism attribution for every imputed cell in one retained row."""

    row_id: str
    knn_imputed_columns: tuple[str, ...] = ()
    minprob_imputed_columns: tuple[str, ...] = ()
    knn_group_labels: tuple[str, ...] = ()
    minprob_group_labels: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        prefix = (
            "dataset processing state missing_data.diagnostics.group_aware.routed_row"
        )
        object.__setattr__(
            self,
            "row_id",
            require_required_str(self.row_id, field_name=f"{prefix}.row_id"),
        )
        for name in ("knn_imputed_columns", "minprob_imputed_columns"):
            values = require_required_label_tuple(
                getattr(self, name), field_name=f"{prefix}.{name}"
            )
            if len(set(values)) != len(values):
                raise PhosPyInputError(f"{prefix}.{name} must contain unique values")
            object.__setattr__(self, name, values)
        for name in ("knn_group_labels", "minprob_group_labels"):
            values = require_required_string_tuple(
                getattr(self, name), field_name=f"{prefix}.{name}"
            )
            if len(set(values)) != len(values):
                raise PhosPyInputError(f"{prefix}.{name} must contain unique values")
            object.__setattr__(self, name, values)
        if not self.knn_imputed_columns and not self.minprob_imputed_columns:
            raise PhosPyInputError(
                f"{prefix} must identify at least one mechanism-attributed cell"
            )
        if set(self.knn_imputed_columns).intersection(self.minprob_imputed_columns):
            raise PhosPyInputError(
                f"{prefix} must not attribute one cell to both mechanisms"
            )
        if bool(self.knn_imputed_columns) != bool(self.knn_group_labels):
            raise PhosPyInputError(
                f"{prefix} KNN columns and group labels must either both be present or absent"
            )
        if bool(self.minprob_imputed_columns) != bool(self.minprob_group_labels):
            raise PhosPyInputError(
                f"{prefix} MinProb columns and group labels must either both be present or absent"
            )
        if set(self.knn_group_labels).intersection(self.minprob_group_labels):
            raise PhosPyInputError(
                f"{prefix} must not route one group to both mechanisms"
            )

    @classmethod
    def from_payload(
        cls, payload: object, *, field_name: str
    ) -> GroupAwareRoutedRowRecord:
        mapping = require_mapping(payload, field_name=field_name)
        require_string_keys(mapping, field_name=field_name)
        unknown = sorted(set(mapping) - _ROUTED_ROW_KNOWN_FIELDS)
        if unknown:
            raise PhosPyInputError(
                f"{field_name} contains unsupported field(s): " + ", ".join(unknown)
            )
        return cls(
            row_id=require_required_str(
                mapping.get("row_id"), field_name=f"{field_name}.row_id"
            ),
            knn_imputed_columns=require_required_label_tuple(
                mapping.get("knn_imputed_columns"),
                field_name=f"{field_name}.knn_imputed_columns",
            ),
            minprob_imputed_columns=require_required_label_tuple(
                mapping.get("minprob_imputed_columns"),
                field_name=f"{field_name}.minprob_imputed_columns",
            ),
            knn_group_labels=require_required_string_tuple(
                mapping.get("knn_group_labels"),
                field_name=f"{field_name}.knn_group_labels",
            ),
            minprob_group_labels=require_required_string_tuple(
                mapping.get("minprob_group_labels"),
                field_name=f"{field_name}.minprob_group_labels",
            ),
        )

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "row_id": self.row_id,
            "knn_imputed_columns": list(self.knn_imputed_columns),
            "minprob_imputed_columns": list(self.minprob_imputed_columns),
            "knn_group_labels": list(self.knn_group_labels),
            "minprob_group_labels": list(self.minprob_group_labels),
        }


@dataclass(frozen=True, slots=True)
class GroupAwareRejectedRowRecord:
    """Typed unsupported-route facts for one rejected input row."""

    row_id: str
    unsupported_partial_groups: tuple[str, ...] = ()
    unsupported_absence_groups: tuple[str, ...] = ()
    observed_finite_count_by_group: Mapping[str, object] = field(default_factory=dict)
    observed_fraction_by_group: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        prefix = (
            "dataset processing state missing_data.diagnostics.group_aware.rejected_row"
        )
        object.__setattr__(
            self,
            "row_id",
            require_required_str(self.row_id, field_name=f"{prefix}.row_id"),
        )
        for name in ("unsupported_partial_groups", "unsupported_absence_groups"):
            object.__setattr__(
                self,
                name,
                require_required_string_tuple(
                    getattr(self, name), field_name=f"{prefix}.{name}"
                ),
            )
        if len(set(self.unsupported_partial_groups)) != len(
            self.unsupported_partial_groups
        ):
            raise PhosPyInputError(
                f"{prefix}.unsupported_partial_groups must contain unique groups"
            )
        if len(set(self.unsupported_absence_groups)) != len(
            self.unsupported_absence_groups
        ):
            raise PhosPyInputError(
                f"{prefix}.unsupported_absence_groups must contain unique groups"
            )
        if set(self.unsupported_partial_groups).intersection(
            self.unsupported_absence_groups
        ):
            raise PhosPyInputError(
                f"{prefix} must not list a group in both unsupported categories"
            )
        if not self.unsupported_partial_groups and not self.unsupported_absence_groups:
            raise PhosPyInputError(
                f"{prefix} must identify at least one unsupported group"
            )
        count_mapping = require_mapping(
            self.observed_finite_count_by_group,
            field_name=f"{prefix}.observed_finite_count_by_group",
        )
        fraction_mapping = require_mapping(
            self.observed_fraction_by_group,
            field_name=f"{prefix}.observed_fraction_by_group",
        )
        require_string_keys(
            count_mapping, field_name=f"{prefix}.observed_finite_count_by_group"
        )
        require_string_keys(
            fraction_mapping, field_name=f"{prefix}.observed_fraction_by_group"
        )
        counts: list[tuple[str, int]] = []
        fractions: list[tuple[str, float]] = []
        for group, value in count_mapping.items():
            counts.append(
                (
                    require_required_str(
                        group,
                        field_name=f"{prefix}.observed_finite_count_by_group.<key>",
                    ),
                    require_required_non_negative_int(
                        value,
                        field_name=f"{prefix}.observed_finite_count_by_group.{group}",
                    ),
                )
            )
        for group, value in fraction_mapping.items():
            parsed_group = require_required_str(
                group, field_name=f"{prefix}.observed_fraction_by_group.<key>"
            )
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise PhosPyInputError(
                    f"{prefix}.observed_fraction_by_group.{group} must be a float"
                )
            parsed_fraction = float(value)
            if not math.isfinite(parsed_fraction) or not 0.0 <= parsed_fraction <= 1.0:
                raise PhosPyInputError(
                    f"{prefix}.observed_fraction_by_group.{group} must satisfy 0 <= value <= 1"
                )
            fractions.append((parsed_group, parsed_fraction))
        if {group for group, _ in counts} != {group for group, _ in fractions}:
            raise PhosPyInputError(
                f"{prefix} observed count and fraction groups must match"
            )
        object.__setattr__(
            self,
            "observed_finite_count_by_group",
            FrozenJsonMapping(
                counts, field_name=f"{prefix}.observed_finite_count_by_group"
            ),
        )
        object.__setattr__(
            self,
            "observed_fraction_by_group",
            FrozenJsonMapping(
                fractions, field_name=f"{prefix}.observed_fraction_by_group"
            ),
        )

    @classmethod
    def from_payload(
        cls, payload: object, *, field_name: str
    ) -> GroupAwareRejectedRowRecord:
        mapping = require_mapping(payload, field_name=field_name)
        require_string_keys(mapping, field_name=field_name)
        unknown = sorted(set(mapping) - _REJECTED_ROW_KNOWN_FIELDS)
        if unknown:
            raise PhosPyInputError(
                f"{field_name} contains unsupported field(s): " + ", ".join(unknown)
            )
        record = cls(
            row_id=require_required_str(
                mapping.get("row_id"), field_name=f"{field_name}.row_id"
            ),
            unsupported_partial_groups=require_required_string_tuple(
                mapping.get("unsupported_partial_groups"),
                field_name=f"{field_name}.unsupported_partial_groups",
            ),
            unsupported_absence_groups=require_required_string_tuple(
                mapping.get("unsupported_absence_groups"),
                field_name=f"{field_name}.unsupported_absence_groups",
            ),
            observed_finite_count_by_group=require_mapping(
                mapping.get("observed_finite_count_by_group"),
                field_name=f"{field_name}.observed_finite_count_by_group",
            ),
            observed_fraction_by_group=require_mapping(
                mapping.get("observed_fraction_by_group"),
                field_name=f"{field_name}.observed_fraction_by_group",
            ),
        )
        return record

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "row_id": self.row_id,
            "unsupported_partial_groups": list(self.unsupported_partial_groups),
            "unsupported_absence_groups": list(self.unsupported_absence_groups),
            "observed_finite_count_by_group": thaw_frozen_json_mapping(
                self.observed_finite_count_by_group,
                field_name="group-aware rejected-row observed counts",
            ),
            "observed_fraction_by_group": thaw_frozen_json_mapping(
                self.observed_fraction_by_group,
                field_name="group-aware rejected-row observed fractions",
            ),
        }


@dataclass(frozen=True, slots=True)
class GroupAwareMissingDataDiagnostics:
    """Typed group-aware routing and mechanism provenance for diagnostics v2."""

    group_column: str
    observed_group_sizes: Mapping[str, object]
    min_partial_observed_fraction: float
    min_reference_observed_fraction: float
    retained_row_count: int
    dropped_unsupported_row_count: int
    knn_routed_cell_count: int
    minprob_routed_cell_count: int
    knn_imputed_cell_count: int
    minprob_imputed_cell_count: int
    knn_target_mask_hash: str
    minprob_target_mask_hash: str
    knn_imputation_mask_hash: str
    minprob_imputation_mask_hash: str
    unsupported_partial_group_count: int
    unsupported_absence_group_count: int
    unsupported_partial_row_ids: tuple[str, ...]
    unsupported_absence_row_ids: tuple[str, ...]
    routed_rows: tuple[GroupAwareRoutedRowRecord, ...]
    rejected_rows: tuple[GroupAwareRejectedRowRecord, ...]
    minprob_left_censored_assumption: bool
    route_categories: tuple[str, ...] = _GROUP_AWARE_ROUTE_CATEGORIES
    mechanism_input: str = "original_retained_matrix"

    def __post_init__(self) -> None:
        prefix = "dataset processing state missing_data.diagnostics.group_aware"
        object.__setattr__(
            self,
            "group_column",
            require_required_str(
                self.group_column, field_name=f"{prefix}.group_column"
            ),
        )
        group_sizes_mapping = require_mapping(
            self.observed_group_sizes, field_name=f"{prefix}.observed_group_sizes"
        )
        require_string_keys(
            group_sizes_mapping, field_name=f"{prefix}.observed_group_sizes"
        )
        group_sizes: list[tuple[str, int]] = []
        for group, size in group_sizes_mapping.items():
            parsed_group = require_required_str(
                group, field_name=f"{prefix}.observed_group_sizes.<key>"
            )
            parsed_size = require_required_non_negative_int(
                size, field_name=f"{prefix}.observed_group_sizes.{parsed_group}"
            )
            if parsed_size < 1:
                raise PhosPyInputError(
                    f"{prefix}.observed_group_sizes.{parsed_group} must be >= 1"
                )
            group_sizes.append((parsed_group, parsed_size))
        if not group_sizes:
            raise PhosPyInputError(f"{prefix}.observed_group_sizes must not be empty")
        object.__setattr__(
            self,
            "observed_group_sizes",
            FrozenJsonMapping(group_sizes, field_name=f"{prefix}.observed_group_sizes"),
        )
        for name in (
            "retained_row_count",
            "dropped_unsupported_row_count",
            "knn_routed_cell_count",
            "minprob_routed_cell_count",
            "knn_imputed_cell_count",
            "minprob_imputed_cell_count",
            "unsupported_partial_group_count",
            "unsupported_absence_group_count",
        ):
            object.__setattr__(
                self,
                name,
                require_required_non_negative_int(
                    getattr(self, name), field_name=f"{prefix}.{name}"
                ),
            )
        for name in (
            "min_partial_observed_fraction",
            "min_reference_observed_fraction",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise PhosPyInputError(f"{prefix}.{name} must be a float")
            resolved = float(value)
            if not (0.0 < resolved <= 1.0):
                raise PhosPyInputError(f"{prefix}.{name} must satisfy 0 < value <= 1")
            object.__setattr__(self, name, resolved)
        for name in (
            "knn_target_mask_hash",
            "minprob_target_mask_hash",
            "knn_imputation_mask_hash",
            "minprob_imputation_mask_hash",
        ):
            object.__setattr__(
                self,
                name,
                require_required_str(
                    getattr(self, name), field_name=f"{prefix}.{name}"
                ),
            )
        for name in ("unsupported_partial_row_ids", "unsupported_absence_row_ids"):
            object.__setattr__(
                self,
                name,
                require_required_string_tuple(
                    getattr(self, name), field_name=f"{prefix}.{name}"
                ),
            )
        records: list[GroupAwareRejectedRowRecord] = []
        if not isinstance(self.rejected_rows, (list, tuple)):
            raise PhosPyInputError(f"{prefix}.rejected_rows must be an array")
        for position, value in enumerate(self.rejected_rows):
            records.append(
                value
                if isinstance(value, GroupAwareRejectedRowRecord)
                else GroupAwareRejectedRowRecord.from_payload(
                    value, field_name=f"{prefix}.rejected_rows[{position}]"
                )
            )
        observed_groups = set(self.observed_group_sizes)
        routed_records: list[GroupAwareRoutedRowRecord] = []
        if not isinstance(self.routed_rows, (list, tuple)):
            raise PhosPyInputError(f"{prefix}.routed_rows must be an array")
        for position, value in enumerate(self.routed_rows):
            routed_records.append(
                value
                if isinstance(value, GroupAwareRoutedRowRecord)
                else GroupAwareRoutedRowRecord.from_payload(
                    value, field_name=f"{prefix}.routed_rows[{position}]"
                )
            )
        routed_row_ids = tuple(record.row_id for record in routed_records)
        if len(set(routed_row_ids)) != len(routed_row_ids):
            raise PhosPyInputError(
                f"{prefix}.routed_rows must contain unique row_id values"
            )
        for record in routed_records:
            unknown_groups = sorted(
                set(record.knn_group_labels)
                .union(record.minprob_group_labels)
                .difference(observed_groups)
            )
            if unknown_groups:
                raise PhosPyInputError(
                    f"{prefix}.routed_rows group(s) must exist in observed_group_sizes: "
                    + ", ".join(unknown_groups)
                )
        object.__setattr__(self, "routed_rows", tuple(routed_records))
        for record in records:
            rejected_groups = set(record.unsupported_partial_groups).union(
                record.unsupported_absence_groups
            )
            unobserved_groups = sorted(rejected_groups - observed_groups)
            if unobserved_groups:
                raise PhosPyInputError(
                    f"{prefix}.rejected_rows group(s) must exist in "
                    "observed_group_sizes: " + ", ".join(unobserved_groups)
                )
            evidence_groups = set(record.observed_finite_count_by_group)
            if evidence_groups != observed_groups:
                raise PhosPyInputError(
                    f"{prefix}.rejected_rows observed evidence groups must match "
                    "observed_group_sizes"
                )
            for group, group_size in self.observed_group_sizes.items():
                parsed_group_size = cast(int, group_size)
                observed_count = cast(int, record.observed_finite_count_by_group[group])
                observed_fraction = cast(
                    float, record.observed_fraction_by_group[group]
                )
                if observed_count > parsed_group_size:
                    raise PhosPyInputError(
                        f"{prefix}.rejected_rows observed count cannot exceed group size"
                    )
                expected_fraction = observed_count / float(parsed_group_size)
                if not math.isclose(
                    observed_fraction, expected_fraction, rel_tol=0.0, abs_tol=1e-12
                ):
                    raise PhosPyInputError(
                        f"{prefix}.rejected_rows observed fractions must match counts"
                    )
            for group in record.unsupported_partial_groups:
                count = cast(int, record.observed_finite_count_by_group[group])
                if not 0 < count < cast(int, self.observed_group_sizes[group]):
                    raise PhosPyInputError(
                        f"{prefix}.rejected_rows unsupported partial groups must be partially observed"
                    )
                if (
                    cast(float, record.observed_fraction_by_group[group])
                    >= self.min_partial_observed_fraction
                ):
                    raise PhosPyInputError(
                        f"{prefix}.rejected_rows unsupported partial groups must be below the threshold"
                    )
            for group in record.unsupported_absence_groups:
                if cast(int, record.observed_finite_count_by_group[group]) != 0:
                    raise PhosPyInputError(
                        f"{prefix}.rejected_rows unsupported absence groups must be fully missing"
                    )
                if any(
                    other_group != group
                    and cast(float, record.observed_fraction_by_group[other_group])
                    >= self.min_reference_observed_fraction
                    for other_group in observed_groups
                ):
                    raise PhosPyInputError(
                        f"{prefix}.rejected_rows unsupported absence groups must lack a sufficiently observed reference"
                    )
        object.__setattr__(self, "rejected_rows", tuple(records))
        if len(records) != self.dropped_unsupported_row_count:
            raise PhosPyInputError(
                f"{prefix}.dropped_unsupported_row_count must match len(rejected_rows)"
            )
        if self.knn_imputed_cell_count != self.knn_routed_cell_count:
            raise PhosPyInputError(
                f"{prefix} KNN routed/imputed cell counts must match"
            )
        if self.minprob_imputed_cell_count != self.minprob_routed_cell_count:
            raise PhosPyInputError(
                f"{prefix} MinProb routed/imputed cell counts must match"
            )
        if sum(len(record.knn_imputed_columns) for record in routed_records) != (
            self.knn_imputed_cell_count
        ):
            raise PhosPyInputError(
                f"{prefix}.routed_rows KNN cells must match knn_imputed_cell_count"
            )
        if sum(len(record.minprob_imputed_columns) for record in routed_records) != (
            self.minprob_imputed_cell_count
        ):
            raise PhosPyInputError(
                f"{prefix}.routed_rows MinProb cells must match minprob_imputed_cell_count"
            )
        if set(routed_row_ids).intersection(record.row_id for record in records):
            raise PhosPyInputError(
                f"{prefix}.routed_rows and rejected_rows must be disjoint"
            )
        route_categories = require_required_string_tuple(
            self.route_categories, field_name=f"{prefix}.route_categories"
        )
        if route_categories != _GROUP_AWARE_ROUTE_CATEGORIES:
            raise PhosPyInputError(
                f"{prefix}.route_categories must use the governed neutral categories"
            )
        object.__setattr__(self, "route_categories", route_categories)
        mechanism_input = require_required_str(
            self.mechanism_input, field_name=f"{prefix}.mechanism_input"
        )
        if mechanism_input != "original_retained_matrix":
            raise PhosPyInputError(
                f"{prefix}.mechanism_input must be 'original_retained_matrix'"
            )
        object.__setattr__(self, "mechanism_input", mechanism_input)
        if not isinstance(self.minprob_left_censored_assumption, bool):
            raise PhosPyInputError(
                f"{prefix}.minprob_left_censored_assumption must be a bool"
            )
        if not self.minprob_left_censored_assumption:
            raise PhosPyInputError(
                f"{prefix}.minprob_left_censored_assumption must be true"
            )
        partial_row_ids = tuple(
            record.row_id for record in records if record.unsupported_partial_groups
        )
        absence_row_ids = tuple(
            record.row_id for record in records if record.unsupported_absence_groups
        )
        partial_group_count = sum(
            len(record.unsupported_partial_groups) for record in records
        )
        absence_group_count = sum(
            len(record.unsupported_absence_groups) for record in records
        )
        if self.unsupported_partial_group_count != partial_group_count:
            raise PhosPyInputError(
                f"{prefix}.unsupported_partial_group_count must match rejected_rows"
            )
        if self.unsupported_absence_group_count != absence_group_count:
            raise PhosPyInputError(
                f"{prefix}.unsupported_absence_group_count must match rejected_rows"
            )
        if self.unsupported_partial_row_ids != partial_row_ids:
            raise PhosPyInputError(
                f"{prefix}.unsupported_partial_row_ids must match rejected_rows"
            )
        if self.unsupported_absence_row_ids != absence_row_ids:
            raise PhosPyInputError(
                f"{prefix}.unsupported_absence_row_ids must match rejected_rows"
            )

    @classmethod
    def from_payload(
        cls, payload: object, *, field_name: str
    ) -> GroupAwareMissingDataDiagnostics:
        mapping = require_mapping(payload, field_name=field_name)
        require_string_keys(mapping, field_name=field_name)
        unknown = sorted(set(mapping) - _GROUP_AWARE_KNOWN_FIELDS)
        if unknown:
            raise PhosPyInputError(
                f"{field_name} contains unsupported field(s): " + ", ".join(unknown)
            )
        rejected = mapping.get("rejected_rows")
        if not isinstance(rejected, (list, tuple)):
            raise PhosPyInputError(f"{field_name}.rejected_rows must be an array")
        routed = mapping.get("routed_rows")
        if not isinstance(routed, (list, tuple)):
            raise PhosPyInputError(f"{field_name}.routed_rows must be an array")
        return cls(
            group_column=require_required_str(
                mapping.get("group_column"), field_name=f"{field_name}.group_column"
            ),
            observed_group_sizes=require_mapping(
                mapping.get("observed_group_sizes"),
                field_name=f"{field_name}.observed_group_sizes",
            ),
            min_partial_observed_fraction=_require_required_fraction(
                mapping.get("min_partial_observed_fraction"),
                field_name=f"{field_name}.min_partial_observed_fraction",
            ),
            min_reference_observed_fraction=_require_required_fraction(
                mapping.get("min_reference_observed_fraction"),
                field_name=f"{field_name}.min_reference_observed_fraction",
            ),
            retained_row_count=require_required_non_negative_int(
                mapping.get("retained_row_count"),
                field_name=f"{field_name}.retained_row_count",
            ),
            dropped_unsupported_row_count=require_required_non_negative_int(
                mapping.get("dropped_unsupported_row_count"),
                field_name=f"{field_name}.dropped_unsupported_row_count",
            ),
            knn_routed_cell_count=require_required_non_negative_int(
                mapping.get("knn_routed_cell_count"),
                field_name=f"{field_name}.knn_routed_cell_count",
            ),
            minprob_routed_cell_count=require_required_non_negative_int(
                mapping.get("minprob_routed_cell_count"),
                field_name=f"{field_name}.minprob_routed_cell_count",
            ),
            knn_imputed_cell_count=require_required_non_negative_int(
                mapping.get("knn_imputed_cell_count"),
                field_name=f"{field_name}.knn_imputed_cell_count",
            ),
            minprob_imputed_cell_count=require_required_non_negative_int(
                mapping.get("minprob_imputed_cell_count"),
                field_name=f"{field_name}.minprob_imputed_cell_count",
            ),
            knn_target_mask_hash=require_required_str(
                mapping.get("knn_target_mask_hash"),
                field_name=f"{field_name}.knn_target_mask_hash",
            ),
            minprob_target_mask_hash=require_required_str(
                mapping.get("minprob_target_mask_hash"),
                field_name=f"{field_name}.minprob_target_mask_hash",
            ),
            knn_imputation_mask_hash=require_required_str(
                mapping.get("knn_imputation_mask_hash"),
                field_name=f"{field_name}.knn_imputation_mask_hash",
            ),
            minprob_imputation_mask_hash=require_required_str(
                mapping.get("minprob_imputation_mask_hash"),
                field_name=f"{field_name}.minprob_imputation_mask_hash",
            ),
            unsupported_partial_group_count=require_required_non_negative_int(
                mapping.get("unsupported_partial_group_count"),
                field_name=f"{field_name}.unsupported_partial_group_count",
            ),
            unsupported_absence_group_count=require_required_non_negative_int(
                mapping.get("unsupported_absence_group_count"),
                field_name=f"{field_name}.unsupported_absence_group_count",
            ),
            unsupported_partial_row_ids=require_required_string_tuple(
                mapping.get("unsupported_partial_row_ids"),
                field_name=f"{field_name}.unsupported_partial_row_ids",
            ),
            unsupported_absence_row_ids=require_required_string_tuple(
                mapping.get("unsupported_absence_row_ids"),
                field_name=f"{field_name}.unsupported_absence_row_ids",
            ),
            routed_rows=tuple(
                GroupAwareRoutedRowRecord.from_payload(
                    value, field_name=f"{field_name}.routed_rows[{position}]"
                )
                for position, value in enumerate(routed)
            ),
            rejected_rows=tuple(
                GroupAwareRejectedRowRecord.from_payload(
                    value, field_name=f"{field_name}.rejected_rows[{position}]"
                )
                for position, value in enumerate(rejected)
            ),
            minprob_left_censored_assumption=_require_required_bool(
                mapping.get("minprob_left_censored_assumption"),
                field_name=f"{field_name}.minprob_left_censored_assumption",
            ),
            route_categories=require_required_string_tuple(
                mapping.get("route_categories"),
                field_name=f"{field_name}.route_categories",
            ),
            mechanism_input=require_required_str(
                mapping.get("mechanism_input"),
                field_name=f"{field_name}.mechanism_input",
            ),
        )

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "group_column": self.group_column,
            "observed_group_sizes": thaw_frozen_json_mapping(
                self.observed_group_sizes,
                field_name="dataset processing state missing_data.diagnostics.group_aware.observed_group_sizes",
            ),
            "min_partial_observed_fraction": self.min_partial_observed_fraction,
            "min_reference_observed_fraction": self.min_reference_observed_fraction,
            "retained_row_count": self.retained_row_count,
            "dropped_unsupported_row_count": self.dropped_unsupported_row_count,
            "knn_routed_cell_count": self.knn_routed_cell_count,
            "minprob_routed_cell_count": self.minprob_routed_cell_count,
            "knn_imputed_cell_count": self.knn_imputed_cell_count,
            "minprob_imputed_cell_count": self.minprob_imputed_cell_count,
            "knn_target_mask_hash": self.knn_target_mask_hash,
            "minprob_target_mask_hash": self.minprob_target_mask_hash,
            "knn_imputation_mask_hash": self.knn_imputation_mask_hash,
            "minprob_imputation_mask_hash": self.minprob_imputation_mask_hash,
            "unsupported_partial_group_count": self.unsupported_partial_group_count,
            "unsupported_absence_group_count": self.unsupported_absence_group_count,
            "unsupported_partial_row_ids": list(self.unsupported_partial_row_ids),
            "unsupported_absence_row_ids": list(self.unsupported_absence_row_ids),
            "routed_rows": [record.to_payload() for record in self.routed_rows],
            "rejected_rows": [record.to_payload() for record in self.rejected_rows],
            "minprob_left_censored_assumption": self.minprob_left_censored_assumption,
            "route_categories": list(self.route_categories),
            "mechanism_input": self.mechanism_input,
        }


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class MissingDataDiagnosticsV2(MissingDataDiagnosticsV1):
    """Diagnostics schema v2 with typed group-aware mechanism provenance."""

    group_aware: GroupAwareMissingDataDiagnostics
    diagnostics_schema_version: int = MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V2
    _expected_schema_version: ClassVar[int] = MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V2

    def __post_init__(self) -> None:
        group_aware = self.group_aware
        if not isinstance(group_aware, GroupAwareMissingDataDiagnostics):
            group_aware = GroupAwareMissingDataDiagnostics.from_payload(
                group_aware,
                field_name="dataset processing state missing_data.diagnostics.group_aware",
            )
        object.__setattr__(self, "group_aware", group_aware)
        super(MissingDataDiagnosticsV2, self).__post_init__()
        if self.missing_data_policy != MissingDataPolicy.IMPUTE_GROUP_AWARE.value:
            raise PhosPyInputError(
                "missing-data diagnostics schema v2 requires "
                "missing_data_policy='impute_group_aware'"
            )
        if self.imputation_method_id != "group_aware_knn_minprob":
            raise PhosPyInputError(
                "missing-data diagnostics v2 requires "
                "imputation_method_id='group_aware_knn_minprob'"
            )
        if self.imputation_method_family != "group_aware_mixed_mechanism":
            raise PhosPyInputError(
                "missing-data diagnostics v2 requires "
                "imputation_method_family='group_aware_mixed_mechanism'"
            )
        if self.imputed_cell_count != (
            group_aware.knn_imputed_cell_count + group_aware.minprob_imputed_cell_count
        ):
            raise PhosPyInputError(
                "missing-data diagnostics v2 mechanism imputation counts must sum "
                "to imputed_cell_count"
            )
        if (
            self.dropped_row_count is None
            or self.dropped_row_count != group_aware.dropped_unsupported_row_count
        ):
            raise PhosPyInputError(
                "missing-data diagnostics v2 dropped counts must agree"
            )
        if self.left_censored_assumption is not True:
            raise PhosPyInputError(
                "missing-data diagnostics v2 must record the configured MinProb "
                "left-censored missingness assumption"
            )
        if self.random_seed is None:
            raise PhosPyInputError("missing-data diagnostics v2 requires random_seed")
        if self.random_seed < 0:
            raise PhosPyInputError(
                "missing-data diagnostics v2 random_seed must be >= 0"
            )
        if (
            self.neighbour_count is None
            or self.distance_metric != "nan_euclidean"
            or self.knn_no_overlap_policy != "error"
            or self.knn_no_overlap_policy_version
            != DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_VERSION
        ):
            raise PhosPyInputError(
                "missing-data diagnostics v2 requires KNN parameters "
                "neighbour_count, distance_metric='nan_euclidean', "
                "knn_no_overlap_policy='error', and the current policy version"
            )
        if self.imputation_mask_hash is None:
            raise PhosPyInputError(
                "missing-data diagnostics v2 requires imputation_mask_hash"
            )
        if self.matrix_scale_requirement != "log2":
            raise PhosPyInputError(
                "missing-data diagnostics v2 requires matrix_scale_requirement='log2'"
            )
        if self.imputation_input_scale != "log2":
            raise PhosPyInputError(
                "missing-data diagnostics v2 requires imputation_input_scale='log2'"
            )
        if (
            self.imputation_input_scale_source
            != IMPUTATION_INPUT_SCALE_SOURCE_METHOD_REQUIRED
        ):
            raise PhosPyInputError(
                "missing-data diagnostics v2 requires "
                "imputation_input_scale_source='method_required'"
            )
        _validate_group_aware_stage_order(
            stage_order=self.stage_order,
            operation_order=self.imputation_operation_order,
        )
        rejected_row_ids = tuple(record.row_id for record in group_aware.rejected_rows)
        if len(set(rejected_row_ids)) != len(rejected_row_ids):
            raise PhosPyInputError(
                "missing-data diagnostics v2 rejected_rows must contain unique row_id values"
            )
        imputed_rejected_row_ids = sorted(
            set(self.imputed_row_ids).intersection(rejected_row_ids)
        )
        if imputed_rejected_row_ids:
            raise PhosPyInputError(
                "missing-data diagnostics v2 imputed_row_ids must be disjoint from "
                "rejected/dropped rows: " + ", ".join(imputed_rejected_row_ids)
            )
        if self.dropped_row_ids != rejected_row_ids:
            raise PhosPyInputError(
                "missing-data diagnostics v2 dropped_row_ids must match "
                "group_aware.rejected_rows"
            )
        if self.rows_not_imputable != rejected_row_ids:
            raise PhosPyInputError(
                "missing-data diagnostics v2 rows_not_imputable must match "
                "group_aware.rejected_rows"
            )
        if self.dropped_rows_above_max_missing_fraction != ():
            raise PhosPyInputError(
                "missing-data diagnostics v2 requires "
                "dropped_rows_above_max_missing_fraction to be empty"
            )
        routed_row_ids = tuple(record.row_id for record in group_aware.routed_rows)
        routed_column_ids = tuple(
            dict.fromkeys(
                column
                for record in group_aware.routed_rows
                for column in (
                    *record.knn_imputed_columns,
                    *record.minprob_imputed_columns,
                )
            )
        )
        if self.imputed_row_ids != routed_row_ids:
            raise PhosPyInputError(
                "missing-data diagnostics v2 imputed_row_ids must match "
                "group_aware.routed_rows"
            )
        if set(self.imputed_column_ids) != set(routed_column_ids):
            raise PhosPyInputError(
                "missing-data diagnostics v2 imputed_column_ids must match "
                "group_aware.routed_rows"
            )
        _validate_group_aware_method_parameters(
            method_parameters=self.method_parameters,
            group_aware=group_aware,
            random_seed=self.random_seed,
            neighbour_count=self.neighbour_count,
            distance_metric=self.distance_metric,
            no_overlap_policy=self.knn_no_overlap_policy,
            no_overlap_policy_version=self.knn_no_overlap_policy_version,
            imputation_input_scale=self.imputation_input_scale,
            imputation_operation_order=self.imputation_operation_order,
        )
        _validate_group_aware_distribution_parameters(
            parameters=self.per_column_distribution_parameters,
            method_parameters=self.method_parameters,
        )

    @classmethod
    def from_mapping(
        cls, payload: Mapping[str, object], *, field_name: str
    ) -> MissingDataDiagnosticsV2:
        require_string_keys(payload, field_name=field_name)
        version = require_int(
            payload.get("diagnostics_schema_version"),
            field_name=f"{field_name}.diagnostics_schema_version",
        )
        if version != MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V2:
            raise PhosPyInputError(
                f"{field_name}.diagnostics_schema_version={version!r} is unsupported; expected 2"
            )
        unknown = sorted(set(payload) - V2_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS)
        if unknown:
            raise PhosPyInputError(
                f"{field_name} contains unsupported field(s): " + ", ".join(unknown)
            )
        base_payload = {
            key: value
            for key, value in payload.items()
            if key in V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS
        }
        base_payload["diagnostics_schema_version"] = (
            MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1
        )
        base = MissingDataDiagnosticsV1.from_mapping(
            base_payload, field_name=field_name
        )
        base_values = {
            item.name: getattr(base, item.name)
            for item in fields(MissingDataDiagnosticsV1)
            if item.name != "diagnostics_schema_version"
        }
        return cls(
            **base_values,
            diagnostics_schema_version=version,
            group_aware=GroupAwareMissingDataDiagnostics.from_payload(
                payload.get("group_aware"), field_name=f"{field_name}.group_aware"
            ),
        )

    def to_payload(self) -> dict[str, JsonValue]:
        payload = super(MissingDataDiagnosticsV2, self).to_payload()
        payload["group_aware"] = self.group_aware.to_payload()
        return payload


def _require_required_fraction(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PhosPyInputError(f"{field_name} must be a float")
    parsed = float(value)
    if not (0.0 < parsed <= 1.0):
        raise PhosPyInputError(f"{field_name} must satisfy 0 < value <= 1")
    return parsed


def _require_required_bool(value: object, *, field_name: str) -> bool:
    parsed = require_optional_bool(value, field_name=field_name)
    if parsed is None:
        raise PhosPyInputError(f"{field_name} is required")
    return parsed


def _validate_group_aware_stage_order(
    *, stage_order: tuple[str, ...], operation_order: str | None
) -> None:
    prefix = "missing-data diagnostics v2"
    if "missing_data" not in stage_order:
        raise PhosPyInputError(f"{prefix} stage_order must contain 'missing_data'")
    if operation_order not in {
        "after_intensity_transform",
        "no_intensity_transform",
    }:
        raise PhosPyInputError(
            f"{prefix} requires a group-aware imputation_operation_order"
        )
    has_transform = "intensity_transform" in stage_order
    if operation_order == "after_intensity_transform":
        if not has_transform or stage_order.index(
            "intensity_transform"
        ) > stage_order.index("missing_data"):
            raise PhosPyInputError(
                f"{prefix} stage_order must place intensity_transform before missing_data"
            )
    elif has_transform:
        raise PhosPyInputError(
            f"{prefix} stage_order contradicts imputation_operation_order"
        )


def _validate_group_aware_method_parameters(
    *,
    method_parameters: Mapping[str, object],
    group_aware: GroupAwareMissingDataDiagnostics,
    random_seed: int,
    neighbour_count: int,
    distance_metric: str,
    no_overlap_policy: str,
    no_overlap_policy_version: int,
    imputation_input_scale: str,
    imputation_operation_order: str | None,
) -> None:
    prefix = "dataset processing state missing_data.diagnostics.method_parameters"
    _require_bounded_number(
        method_parameters.get("q"),
        field_name=f"{prefix}.q",
        lower_exclusive=0.0,
        upper_exclusive=0.5,
    )
    _require_bounded_number(
        method_parameters.get("width"),
        field_name=f"{prefix}.width",
        lower_exclusive=0.0,
        upper_inclusive=1.0,
    )
    _require_bounded_number(
        method_parameters.get("min_partial_observed_fraction"),
        field_name=f"{prefix}.min_partial_observed_fraction",
        lower_exclusive=0.0,
        upper_inclusive=1.0,
    )
    _require_bounded_number(
        method_parameters.get("min_reference_observed_fraction"),
        field_name=f"{prefix}.min_reference_observed_fraction",
        lower_exclusive=0.0,
        upper_inclusive=1.0,
    )
    for count_name in ("knn_target_cell_count", "minprob_target_cell_count"):
        require_required_non_negative_int(
            method_parameters.get(count_name), field_name=f"{prefix}.{count_name}"
        )
    expected: tuple[tuple[str, object], ...] = (
        ("seed", random_seed),
        ("k", neighbour_count),
        ("distance", distance_metric),
        ("no_overlap_policy", no_overlap_policy),
        ("no_overlap_policy_version", no_overlap_policy_version),
        ("group_column", group_aware.group_column),
        (
            "min_partial_observed_fraction",
            group_aware.min_partial_observed_fraction,
        ),
        (
            "min_reference_observed_fraction",
            group_aware.min_reference_observed_fraction,
        ),
        ("knn_target_cell_count", group_aware.knn_routed_cell_count),
        ("minprob_target_cell_count", group_aware.minprob_routed_cell_count),
        ("knn_target_mask_hash", group_aware.knn_target_mask_hash),
        ("minprob_target_mask_hash", group_aware.minprob_target_mask_hash),
        ("mechanism_input", group_aware.mechanism_input),
        ("input_scale", imputation_input_scale),
        ("imputation_operation_order", imputation_operation_order),
    )
    for name, expected_value in expected:
        if name not in method_parameters:
            raise PhosPyInputError(f"{prefix}.{name} is required")
        if method_parameters[name] != expected_value:
            raise PhosPyInputError(f"{prefix}.{name} must agree with typed diagnostics")
    if isinstance(method_parameters["seed"], bool) or not isinstance(
        method_parameters["seed"], int
    ):
        raise PhosPyInputError(f"{prefix}.seed must be an int")
    if isinstance(method_parameters["k"], bool) or not isinstance(
        method_parameters["k"], int
    ):
        raise PhosPyInputError(f"{prefix}.k must be an int")
    if isinstance(
        method_parameters["no_overlap_policy_version"], bool
    ) or not isinstance(method_parameters["no_overlap_policy_version"], int):
        raise PhosPyInputError(f"{prefix}.no_overlap_policy_version must be an int")
    resolved_groups = require_mapping(
        method_parameters.get("resolved_group_samples"),
        field_name=f"{prefix}.resolved_group_samples",
    )
    require_string_keys(resolved_groups, field_name=f"{prefix}.resolved_group_samples")
    if set(resolved_groups) != set(group_aware.observed_group_sizes):
        raise PhosPyInputError(
            f"{prefix}.resolved_group_samples groups must match observed_group_sizes"
        )
    all_samples: list[str] = []
    for group, raw_samples in resolved_groups.items():
        if not isinstance(raw_samples, (list, tuple)):
            raise PhosPyInputError(
                f"{prefix}.resolved_group_samples.{group} must be an array"
            )
        samples = require_required_label_tuple(
            raw_samples, field_name=f"{prefix}.resolved_group_samples.{group}"
        )
        if len(samples) != group_aware.observed_group_sizes[group]:
            raise PhosPyInputError(
                f"{prefix}.resolved_group_samples.{group} size must match "
                "observed_group_sizes"
            )
        all_samples.extend(samples)
    if len(set(all_samples)) != len(all_samples):
        raise PhosPyInputError(
            f"{prefix}.resolved_group_samples must not assign a sample twice"
        )
    sample_to_group = {
        sample: group
        for group, raw_samples in resolved_groups.items()
        for sample in require_required_label_tuple(
            raw_samples, field_name=f"{prefix}.resolved_group_samples.{group}"
        )
    }
    for record in group_aware.routed_rows:
        for mechanism, columns, groups in (
            (
                "KNN",
                record.knn_imputed_columns,
                record.knn_group_labels,
            ),
            (
                "MinProb",
                record.minprob_imputed_columns,
                record.minprob_group_labels,
            ),
        ):
            unknown_columns = sorted(set(columns) - set(sample_to_group))
            if unknown_columns:
                raise PhosPyInputError(
                    f"{prefix} {mechanism} routed columns must exist in "
                    "resolved_group_samples: " + ", ".join(unknown_columns)
                )
            expected_groups = tuple(
                dict.fromkeys(sample_to_group[column] for column in columns)
            )
            if groups != expected_groups:
                raise PhosPyInputError(
                    f"{prefix} {mechanism} routed groups must match routed columns"
                )
            columns_by_group = {
                group: {
                    column for column in columns if sample_to_group[column] == group
                }
                for group in groups
            }
            for group, routed_columns in columns_by_group.items():
                group_columns = set(
                    require_required_label_tuple(
                        resolved_groups[group],
                        field_name=f"{prefix}.resolved_group_samples.{group}",
                    )
                )
                if mechanism == "KNN" and routed_columns == group_columns:
                    raise PhosPyInputError(
                        f"{prefix} KNN routed groups must remain partially observed"
                    )
                if mechanism == "MinProb" and routed_columns != group_columns:
                    raise PhosPyInputError(
                        f"{prefix} MinProb routed groups must include every group sample"
                    )


def _validate_group_aware_distribution_parameters(
    *,
    parameters: Mapping[str, object] | None,
    method_parameters: Mapping[str, object],
) -> None:
    prefix = (
        "dataset processing state missing_data.diagnostics."
        "per_column_distribution_parameters"
    )
    if parameters is None:
        raise PhosPyInputError(
            "missing-data diagnostics v2 requires per-column MinProb distribution provenance"
        )
    resolved_groups = require_mapping(
        method_parameters.get("resolved_group_samples"),
        field_name=(
            "dataset processing state missing_data.diagnostics.method_parameters."
            "resolved_group_samples"
        ),
    )
    expected_columns = tuple(
        sample
        for samples in resolved_groups.values()
        for sample in require_required_label_tuple(
            samples, field_name="resolved_group_samples.<group>"
        )
    )
    if len(parameters) != len(expected_columns) or set(parameters) != set(
        expected_columns
    ):
        raise PhosPyInputError(f"{prefix} columns must match resolved_group_samples")
    q = _require_bounded_number(
        method_parameters.get("q"),
        field_name="missing-data diagnostics v2 method_parameters.q",
        lower_exclusive=0.0,
        upper_exclusive=0.5,
    )
    width = _require_bounded_number(
        method_parameters.get("width"),
        field_name="missing-data diagnostics v2 method_parameters.width",
        lower_exclusive=0.0,
        upper_inclusive=1.0,
    )
    required_numeric = (
        "lower_q_quantile",
        "lower_tail_mean",
        "observed_sd",
        "imputation_mean",
        "imputation_sd",
    )
    for column, raw_values in parameters.items():
        values = require_mapping(raw_values, field_name=f"{prefix}.{column}")
        parsed_numeric: dict[str, float] = {}
        for count_name in ("observed_count", "missing_count"):
            require_required_non_negative_int(
                values.get(count_name), field_name=f"{prefix}.{column}.{count_name}"
            )
        column_q = _require_bounded_number(
            values.get("q"),
            field_name=f"{prefix}.{column}.q",
            lower_exclusive=0.0,
            upper_exclusive=0.5,
        )
        column_width = _require_bounded_number(
            values.get("width"),
            field_name=f"{prefix}.{column}.width",
            lower_exclusive=0.0,
            upper_inclusive=1.0,
        )
        if column_q != q:
            raise PhosPyInputError(f"{prefix}.{column}.q must agree with method q")
        if column_width != width:
            raise PhosPyInputError(
                f"{prefix}.{column}.width must agree with method width"
            )
        for name in required_numeric:
            value = values.get(name)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise PhosPyInputError(f"{prefix}.{column}.{name} must be a number")
            parsed_value = float(value)
            if not math.isfinite(parsed_value):
                raise PhosPyInputError(f"{prefix}.{column}.{name} must be finite")
            parsed_numeric[name] = parsed_value
        if (
            parsed_numeric["observed_sd"] <= 0.0
            or parsed_numeric["imputation_sd"] <= 0.0
        ):
            raise PhosPyInputError(
                f"{prefix}.{column} distribution standard deviations must be > 0"
            )


def _require_bounded_number(
    value: object,
    *,
    field_name: str,
    lower_exclusive: float,
    upper_exclusive: float | None = None,
    upper_inclusive: float | None = None,
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PhosPyInputError(f"{field_name} must be a number")
    parsed = float(value)
    valid = math.isfinite(parsed) and parsed > lower_exclusive
    if upper_exclusive is not None:
        valid = valid and parsed < upper_exclusive
    if upper_inclusive is not None:
        valid = valid and parsed <= upper_inclusive
    if not valid:
        raise PhosPyInputError(f"{field_name} is outside its supported domain")
    return parsed


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

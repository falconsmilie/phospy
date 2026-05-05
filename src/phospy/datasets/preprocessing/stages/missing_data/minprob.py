"""MinProb missing-data policy implementation."""

from __future__ import annotations

import datetime as dt
import hashlib
import json

import numpy as np

from phospy.datasets.preprocessing.models import PreprocessingState
from phospy.datasets.processing_state import JsonValue
from phospy.errors.input import PhosPyInputError

from .models import MinProbPolicyOutcome, RowImputationRecord


def _normalise_column_label_value(label: object) -> object:
    """Return a deterministic, JSON-serialisable representation for a label."""

    if isinstance(label, np.generic):
        return _normalise_column_label_value(label.item())
    if label is None:
        return {"kind": "none"}
    if isinstance(label, bool):
        return {"kind": "bool", "value": bool(label)}
    if isinstance(label, int):
        return {"kind": "int", "value": int(label)}
    if isinstance(label, float):
        value = float(label)
        if np.isnan(value):
            return {"kind": "float", "value": "nan"}
        if np.isposinf(value):
            return {"kind": "float", "value": "inf"}
        if np.isneginf(value):
            return {"kind": "float", "value": "-inf"}
        return {"kind": "float", "value": value.hex()}
    if isinstance(label, str):
        return {"kind": "str", "value": label}
    if isinstance(label, bytes):
        return {"kind": "bytes", "value": label.hex()}
    if isinstance(label, dt.datetime):
        return {"kind": "datetime", "value": label.isoformat()}
    if isinstance(label, dt.date):
        return {"kind": "date", "value": label.isoformat()}
    if isinstance(label, dt.time):
        return {"kind": "time", "value": label.isoformat()}
    if isinstance(label, dt.timedelta):
        return {"kind": "timedelta", "value": str(label)}
    if isinstance(label, tuple):
        return {
            "kind": "tuple",
            "value": [_normalise_column_label_value(item) for item in label],
        }
    return {
        "kind": "fallback",
        "type": f"{type(label).__module__}.{type(label).__qualname__}",
        "repr": repr(label),
    }


def _stable_column_label_seed(base_seed: int, column_label: object) -> int:
    """Derive a deterministic per-column seed from base seed + label digest."""

    label_payload = _normalise_column_label_value(column_label)
    label_bytes = json.dumps(
        label_payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    digest = hashlib.blake2b(
        label_bytes,
        digest_size=16,
        person=b"phospy-minprob",
    ).digest()
    return int(base_seed) + int.from_bytes(digest, byteorder="big", signed=False)


def run_minprob_policy(state: PreprocessingState) -> MinProbPolicyOutcome:
    """Apply minprob policy numerical transformation."""

    if state.plan.intensity_transform_policy != "log2":
        raise PhosPyInputError(
            "dataset build request preprocessing_config.missing_data.policy="
            "'impute_minprob' requires log2-scale values. Configure "
            "preprocessing_config.intensity_transform.policy='log2'."
        )

    q = state.plan.missing_data_q
    width = state.plan.missing_data_width
    seed = state.plan.missing_data_seed
    max_missing_fraction_per_row = state.plan.missing_data_max_missing_fraction_per_row
    if (
        q is None
        or width is None
        or seed is None
        or max_missing_fraction_per_row is None
    ):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.missing_data.policy="
            "'impute_minprob' requires q, width, seed, and "
            "max_missing_fraction_per_row"
        )
    q_value = float(q)
    width_value = float(width)
    seed_value = int(seed)
    max_missing_fraction_value = float(max_missing_fraction_per_row)

    missing_fraction = state.phospho.isna().mean(axis=1)
    retained_mask = missing_fraction <= max_missing_fraction_value
    dropped_missing_fraction = missing_fraction.loc[~retained_mask]
    filtered_phospho = state.phospho.loc[retained_mask].copy(deep=True)
    filtered_site_metadata = state.site_metadata.loc[filtered_phospho.index]

    eps = float(np.finfo(float).eps)
    per_column_distribution_parameters: dict[str, dict[str, JsonValue]] = {}
    for column_name in filtered_phospho.columns:
        column_label = str(column_name)
        column = filtered_phospho.loc[:, column_name]
        missing_mask = column.isna().to_numpy(dtype=bool, copy=False)
        missing_count = int(missing_mask.sum())
        observed_values = column.dropna().to_numpy(dtype=float, copy=False)
        observed_count = int(observed_values.size)
        if observed_count == 0 and missing_count > 0:
            raise PhosPyInputError(
                "dataset preprocessing stage 'missing_data' cannot apply "
                "missing_data.policy='impute_minprob' because column "
                f"{column_label!r} has no observed values after row filtering; "
                "adjust missing_data.max_missing_fraction_per_row or input data."
            )

        quantile_value = (
            float(np.quantile(observed_values, q_value)) if observed_count > 0 else 0.0
        )
        if observed_count > 1:
            observed_sd = float(np.std(observed_values, ddof=1))
        elif observed_count == 1:
            observed_sd = 0.0
        else:
            observed_sd = 0.0
        if not np.isfinite(observed_sd) or observed_sd <= 0.0:
            observed_sd = (
                float(np.std(observed_values, ddof=0)) if observed_count > 0 else 0.0
            )
        if not np.isfinite(observed_sd) or observed_sd <= 0.0:
            observed_sd = eps

        imputation_sd = max(observed_sd * width_value, eps)
        imputation_mean = float(quantile_value - (1.8 * imputation_sd))
        lower_tail = observed_values[observed_values <= quantile_value]
        lower_tail_mean = (
            float(np.mean(lower_tail))
            if int(lower_tail.size) > 0
            else float(quantile_value)
        )

        per_column_distribution_parameters[column_label] = {
            "observed_count": int(observed_count),
            "missing_count": int(missing_count),
            "q": float(q_value),
            "width": float(width_value),
            "lower_q_quantile": float(quantile_value),
            "lower_tail_mean": float(lower_tail_mean),
            "observed_sd": float(observed_sd),
            "imputation_mean": float(imputation_mean),
            "imputation_sd": float(imputation_sd),
        }

        if missing_count == 0:
            continue
        column_rng = np.random.default_rng(
            _stable_column_label_seed(seed_value, column_name)
        )
        draws = column_rng.normal(
            loc=imputation_mean,
            scale=imputation_sd,
            size=missing_count,
        )
        missing_index = filtered_phospho.index[missing_mask]
        filtered_phospho.loc[missing_index, column_name] = draws

    if filtered_phospho.empty:
        imputed_mask = filtered_phospho.isna() & filtered_phospho.notna()
    else:
        imputed_mask = (
            state.phospho.loc[retained_mask].isna() & filtered_phospho.notna()
        )
    imputed_cell_count = int(imputed_mask.to_numpy().sum())
    imputed_row_ids = (
        tuple(
            str(row_id)
            for row_id in filtered_phospho.index[
                imputed_mask.any(axis=1).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not filtered_phospho.empty
        else ()
    )
    imputed_column_ids = (
        tuple(
            str(column_name)
            for column_name in filtered_phospho.columns[
                imputed_mask.any(axis=0).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not filtered_phospho.empty
        else ()
    )
    rows_not_imputable = (
        tuple(
            str(row_id)
            for row_id in filtered_phospho.index[
                filtered_phospho.isna().any(axis=1).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not filtered_phospho.empty
        else ()
    )
    dropped_row_ids = tuple(
        str(row_id) for row_id in dropped_missing_fraction.index.tolist()
    )
    output_missing_cell_count = int(filtered_phospho.isna().to_numpy().sum())
    imputed_rows = (
        tuple(
            RowImputationRecord(
                row_id=str(row_id),
                imputed_columns=tuple(
                    str(column_name)
                    for column_name in filtered_phospho.columns[
                        imputed_mask.loc[row_id]
                    ].tolist()
                ),
                imputed_cell_count=int(imputed_mask.loc[row_id].sum()),
            )
            for row_id in filtered_phospho.index[
                imputed_mask.any(axis=1).to_numpy(dtype=bool, copy=False)
            ]
        )
        if not filtered_phospho.empty
        else ()
    )
    dropped_rows_missing_fraction = tuple(
        (str(row_id), float(missing_fraction_value))
        for row_id, missing_fraction_value in dropped_missing_fraction.items()
    )
    return MinProbPolicyOutcome(
        phospho=filtered_phospho,
        site_metadata=filtered_site_metadata,
        q=q_value,
        width=width_value,
        seed=seed_value,
        max_missing_fraction_per_row=max_missing_fraction_value,
        dropped_row_ids=dropped_row_ids,
        dropped_rows_missing_fraction=dropped_rows_missing_fraction,
        imputed_cell_count=imputed_cell_count,
        imputed_row_ids=imputed_row_ids,
        imputed_column_ids=imputed_column_ids,
        output_missing_cell_count=output_missing_cell_count,
        rows_not_imputable=rows_not_imputable,
        per_column_distribution_parameters=per_column_distribution_parameters,
        imputed_rows=imputed_rows,
    )

"""Private preprocessing stage fingerprint service."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

import pandas as pd

from phospy.errors.build import DatasetBuildError
from phospy.provenance.hashing import (
    fingerprint_optional_table,
    hash_table_tolerance,
)
from phospy.provenance.models import TableFingerprint
from phospy.science.datasets.preprocessing.models import (
    PreprocessingState,
    PreprocessingStateTableKey,
)


@dataclass(frozen=True, slots=True)
class _StageFingerprintBundle:
    input_hash: str
    output_hash: str
    phospho_input_hash: str
    phospho_output_hash: str
    consumed_input_tables: tuple[TableFingerprint, ...]
    produced_output_tables: tuple[TableFingerprint, ...]


_StageTableHashCacheKey = tuple[
    str,
    str,
    tuple[tuple[str, int, int, str, str, str, str], ...],
]


class _StageFingerprintService:
    """Compute stage table and phospho hashes for one preprocessing run."""

    def __init__(self) -> None:
        self._stage_table_hash_cache: dict[_StageTableHashCacheKey, str] = {}

    def run(
        self,
        *,
        stage_key: str,
        previous: PreprocessingState,
        current: PreprocessingState,
        consumed_input_tables: tuple[PreprocessingStateTableKey, ...],
        produced_output_tables: tuple[PreprocessingStateTableKey, ...],
    ) -> _StageFingerprintBundle:
        consumed_fingerprints = _collect_stage_table_fingerprints(
            state=previous,
            table_names=consumed_input_tables,
        )
        produced_fingerprints = _collect_stage_table_fingerprints(
            state=current,
            table_names=produced_output_tables,
        )
        return _StageFingerprintBundle(
            input_hash=self._hash_stage_table_fingerprints(
                stage_key=stage_key,
                direction="input",
                table_fingerprints=consumed_fingerprints,
            ),
            output_hash=self._hash_stage_table_fingerprints(
                stage_key=stage_key,
                direction="output",
                table_fingerprints=produced_fingerprints,
            ),
            phospho_input_hash=hash_table_tolerance(
                previous.phospho,
                name=f"{stage_key}.input.phospho",
            ),
            phospho_output_hash=hash_table_tolerance(
                current.phospho,
                name=f"{stage_key}.output.phospho",
            ),
            consumed_input_tables=consumed_fingerprints,
            produced_output_tables=produced_fingerprints,
        )

    def _hash_stage_table_fingerprints(
        self,
        *,
        stage_key: str,
        direction: str,
        table_fingerprints: tuple[TableFingerprint, ...],
    ) -> str:
        cache_key = (
            stage_key,
            direction,
            tuple(_fingerprint_cache_token(item) for item in table_fingerprints),
        )
        cached = self._stage_table_hash_cache.get(cache_key)
        if cached is not None:
            return cached
        digest = _hash_stage_table_fingerprints(
            stage_key=stage_key,
            direction=direction,
            table_fingerprints=table_fingerprints,
        )
        self._stage_table_hash_cache[cache_key] = digest
        return digest


def _fingerprint_cache_token(
    fingerprint: TableFingerprint,
) -> tuple[str, int, int, str, str, str, str]:
    return (
        fingerprint.name,
        int(fingerprint.rows),
        int(fingerprint.columns),
        fingerprint.exact_hash_algorithm,
        fingerprint.exact_hash_value,
        fingerprint.tolerance_hash_algorithm,
        fingerprint.tolerance_hash_value,
    )


def _hash_stage_table_fingerprints(
    *,
    stage_key: str,
    direction: str,
    table_fingerprints: tuple[TableFingerprint, ...],
) -> str:
    payload = {
        "stage": stage_key,
        "direction": direction,
        "tables": [
            {
                "name": item.name,
                "rows": int(item.rows),
                "columns": int(item.columns),
                "exact_hash_algorithm": item.exact_hash_algorithm,
                "exact_hash_value": item.exact_hash_value,
                "tolerance_hash_algorithm": item.tolerance_hash_algorithm,
                "tolerance_hash_value": item.tolerance_hash_value,
            }
            for item in table_fingerprints
        ],
    }
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _collect_stage_table_fingerprints(
    *,
    state: PreprocessingState,
    table_names: tuple[PreprocessingStateTableKey, ...],
) -> tuple[TableFingerprint, ...]:
    fingerprints: list[TableFingerprint] = []
    for table_name in table_names:
        table = _resolve_state_table(state=state, table_name=table_name)
        fingerprint = fingerprint_optional_table(table, name=table_name.value)
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _resolve_state_table(
    *,
    state: PreprocessingState,
    table_name: PreprocessingStateTableKey | str,
) -> pd.DataFrame | None:
    try:
        key = (
            table_name
            if isinstance(table_name, PreprocessingStateTableKey)
            else PreprocessingStateTableKey(str(table_name))
        )
    except ValueError as exc:
        supported = ", ".join(item.value for item in PreprocessingStateTableKey)
        raise DatasetBuildError(
            "dataset preprocessing stage metadata contains unknown table key: "
            f"{table_name!r}; supported tables: {supported}"
        ) from exc
    if key is PreprocessingStateTableKey.DATASET_PHOSPHO:
        return state.phospho
    if key is PreprocessingStateTableKey.DATASET_SITE_METADATA:
        return state.site_metadata
    if key is PreprocessingStateTableKey.DATASET_SAMPLE_METADATA:
        return state.sample_metadata
    if key is PreprocessingStateTableKey.DATASET_TOTAL:
        return state.total
    if key is PreprocessingStateTableKey.DATASET_COMPARISONS:
        return state.comparisons
    if key is PreprocessingStateTableKey.DATASET_IMPUTATION_OBSERVATION_MASK:
        return state.imputation_observation_mask
    if key is PreprocessingStateTableKey.REPORT_COMPARISON_GROUP_STATS:
        return state.comparison_group_stats
    if key is PreprocessingStateTableKey.REPORT_COMPARISON_PAIR_STATS:
        return state.comparison_pair_stats
    if key is PreprocessingStateTableKey.REPORT_DUPLICATE_SITE_RESOLUTION:
        return state.duplicate_site_resolution
    if key is PreprocessingStateTableKey.REPORT_METADATA_CONFLICTS:
        return state.metadata_conflicts
    if key is PreprocessingStateTableKey.REPORT_ROW_AUDIT:
        return state.row_audit
    raise DatasetBuildError(
        "dataset preprocessing stage metadata references an unsupported table key: "
        f"{key.value!r}"
    )

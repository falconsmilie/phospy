"""Pure Kinase Library-style sequence motif scoring."""

from __future__ import annotations

import math
import re
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, cast

import numpy as np
import pandas as pd

from phospy.science.prediction.motif_scoring.models import (
    AMINO_ACIDS,
    DEFAULT_MOTIF_FLANK_SIZE,
    KINASE_LIBRARY_MATRIX_STATUS_DUPLICATE,
    KINASE_LIBRARY_MATRIX_STATUS_FILTERED_RESIDUE_CLASS,
    KINASE_LIBRARY_MATRIX_STATUS_INVALID,
    KINASE_LIBRARY_MATRIX_STATUS_UNSUPPORTED_WINDOW,
    KINASE_LIBRARY_MATRIX_STATUS_VALID,
    KINASE_LIBRARY_RESIDUE_CLASS_SER_THR,
    KINASE_LIBRARY_RESIDUE_CLASS_TYR,
    KINASE_LIBRARY_RESIDUE_CLASSES,
    KINASE_LIBRARY_SITE_STATUS_MISSING_SEQUENCE,
    KINASE_LIBRARY_SITE_STATUS_UNSUPPORTED_RESIDUE,
    KINASE_LIBRARY_SITE_STATUS_UNSUPPORTED_SEQUENCE_LENGTH,
    KINASE_LIBRARY_SITE_STATUS_VALID_SCORED_SITE,
    KINASE_LIBRARY_SITE_STATUS_VALID_UNSCORED_SITE,
    KINASE_LIBRARY_SITE_STATUS_WRONG_CENTRAL_RESIDUE,
    KINASE_LIBRARY_SITE_STATUS_WRONG_RESIDUE_CLASS,
    SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
    SEQUENCE_SEMANTICS_CENTRED_WINDOW,
    KinaseLibraryMotifMatrix,
    KinaseLibraryMotifScoringResult,
    KinaseLibraryScoreScaleMetadata,
    KinaseLibraryWindowConfig,
    SequenceSemantics,
    _is_missing_scalar,
    _normalize_sequence_value,
    normalize_kinase_library_residue_class,
)

_SITE_IDENTITY_PATTERN = re.compile(r"^\s*[^;]+?\s*;\s*(?P<site>[^;]+?)\s*;\s*$")
_PHOSPHO_RESIDUE_CLASSES = {
    "S": KINASE_LIBRARY_RESIDUE_CLASS_SER_THR,
    "T": KINASE_LIBRARY_RESIDUE_CLASS_SER_THR,
    "Y": KINASE_LIBRARY_RESIDUE_CLASS_TYR,
}
_SUPPORTED_AMINO_ACID_SET = frozenset(AMINO_ACIDS)
_RAW_SCORE_FORMULA = "sum(matrix[amino_acid_at_relative_position, relative_position])"
_PERCENTILE_METHOD_HIGHER_IS_BETTER = "100 * count(reference_score <= site_score) / n"
_PERCENTILE_METHOD_LOWER_IS_BETTER = "100 * count(reference_score >= site_score) / n"
_RANK_METHOD_HIGHER_IS_BETTER = "1 + count(reference_score > site_score)"
_RANK_METHOD_LOWER_IS_BETTER = "1 + count(reference_score < site_score)"


@dataclass(frozen=True, slots=True)
class _MatrixLane:
    kinase: str
    residue_class: str
    score_table: pd.DataFrame
    kinase_family: str | None
    kinase_group: str | None
    source_order: int


@dataclass(frozen=True, slots=True)
class _SiteValidation:
    site_id: str
    sequence: str | None
    sequence_window: str | None
    sequence_length: int | None
    observed_central_residue: str | None
    expected_central_residue: str | None
    residue_class: str | None
    status: str
    reason: str | None


class KinaseLibraryMotifScorer:
    """Run pure Kinase Library-style motif scoring without workflow coupling."""

    def run(
        self,
        *,
        site_sequences: Mapping[str, object] | Sequence[object] | pd.Series,
        matrices: object,
        residue_classes: Sequence[object] | None = None,
        window: KinaseLibraryWindowConfig | object | None = None,
        upstream_residues: int | None = None,
        downstream_residues: int | None = None,
        flank_size: int | None = None,
        sequence_semantics: SequenceSemantics | None = None,
        site_index: Sequence[object] | None = None,
        site_identities: Mapping[str, object] | pd.Series | None = None,
        matrix_residue_classes: Mapping[str, object] | None = None,
        reference_distributions: Mapping[object, object] | None = None,
        score_scale: str | None = None,
        higher_is_better: bool = True,
    ) -> KinaseLibraryMotifScoringResult:
        """Score site sequence windows against Kinase Library-style matrices.

        Raw scores are provider-scale sums over the configured relative positions.
        Invalid or residue-class-incompatible sites remain missing; no neutral
        scores are imputed.
        """

        explicit_classes = _coerce_residue_classes(residue_classes)
        default_matrix_class = (
            explicit_classes[0]
            if explicit_classes and len(explicit_classes) == 1
            else None
        )
        matrix_lanes, resource_window, resource_score_scale = _coerce_matrix_lanes(
            matrices,
            matrix_residue_classes=matrix_residue_classes,
            default_residue_class=default_matrix_class,
        )
        window_config = _resolve_window_config(
            window=window,
            upstream_residues=upstream_residues,
            downstream_residues=downstream_residues,
            flank_size=flank_size,
            sequence_semantics=sequence_semantics,
            matrix_lanes=matrix_lanes,
            resource_window=resource_window,
        )
        allowed_classes = explicit_classes or _classes_from_lanes(matrix_lanes)
        if not allowed_classes:
            allowed_classes = tuple(KINASE_LIBRARY_RESIDUE_CLASSES)

        raw_sequences = _coerce_sequence_series(site_sequences, site_index=site_index)
        resolved_site_identities = _coerce_site_identity_series(
            site_identities,
            site_index=list(raw_sequences.index.astype(str)),
        )
        site_validations = [
            _validate_site_sequence(
                site_id=str(site_id),
                sequence_input=sequence,
                site_identity=resolved_site_identities.get(str(site_id)),
                window_config=window_config,
                allowed_residue_classes=allowed_classes,
            )
            for site_id, sequence in raw_sequences.items()
        ]
        sequence_windows = pd.Series(
            {site.site_id: site.sequence_window for site in site_validations},
            index=raw_sequences.index.copy(),
            dtype=object,
            name="sequence_window",
        )

        valid_lanes, kinase_diagnostics = _validate_matrix_lanes(
            matrix_lanes=matrix_lanes,
            window_config=window_config,
            allowed_residue_classes=allowed_classes,
        )
        kinase_columns = _stable_kinase_columns(valid_lanes)
        raw_scores = pd.DataFrame(
            np.nan,
            index=raw_sequences.index.copy(),
            columns=pd.Index(kinase_columns, name="kinase"),
            dtype=float,
        )
        missing_residue_score_cells: dict[tuple[str, str], int] = defaultdict(int)

        site_by_id = {site.site_id: site for site in site_validations}
        for site_id in raw_scores.index.astype(str).tolist():
            site = site_by_id[site_id]
            if site.status != KINASE_LIBRARY_SITE_STATUS_VALID_UNSCORED_SITE:
                continue
            if site.sequence_window is None or site.residue_class is None:
                continue
            for lane in valid_lanes:
                if lane.residue_class != site.residue_class:
                    continue
                score = _score_sequence_window(
                    sequence_window=site.sequence_window,
                    score_table=lane.score_table,
                    positions=window_config.positions,
                )
                if score is None:
                    missing_residue_score_cells[(lane.kinase, lane.residue_class)] += 1
                    continue
                raw_scores.at[site_id, lane.kinase] = score

        site_diagnostics = _build_site_diagnostics(
            site_validations=site_validations,
            raw_scores=raw_scores,
        )
        kinase_diagnostics = _finish_kinase_diagnostics(
            kinase_diagnostics=kinase_diagnostics,
            site_validations=site_validations,
            raw_scores=raw_scores,
            missing_residue_score_cells=missing_residue_score_cells,
        )

        percentile_ranks: pd.DataFrame | None = None
        reference_ranks: pd.DataFrame | None = None
        distribution_lookup = _coerce_reference_distributions(reference_distributions)
        if reference_distributions is not None:
            percentile_ranks, reference_ranks = _score_reference_distributions(
                raw_scores=raw_scores,
                site_diagnostics=site_diagnostics,
                distribution_lookup=distribution_lookup,
                higher_is_better=bool(higher_is_better),
            )

        resolved_score_scale = (
            str(score_scale)
            if score_scale is not None
            else str(resource_score_scale or "kinase_library_raw_position_sum")
        )
        metadata = KinaseLibraryScoreScaleMetadata(
            score_scale=resolved_score_scale,
            raw_score_formula=_RAW_SCORE_FORMULA,
            higher_is_better=bool(higher_is_better),
            percentile_method=(
                _percentile_method_label(higher_is_better=bool(higher_is_better))
                if reference_distributions is not None
                else None
            ),
            rank_method=(
                _rank_method_label(higher_is_better=bool(higher_is_better))
                if reference_distributions is not None
                else None
            ),
            sequence_window=window_config.to_payload(),
            residue_classes=tuple(allowed_classes),
        )
        return KinaseLibraryMotifScoringResult(
            raw_scores=raw_scores,
            percentile_ranks=percentile_ranks,
            reference_ranks=reference_ranks,
            site_diagnostics=site_diagnostics,
            kinase_diagnostics=kinase_diagnostics,
            sequence_windows=sequence_windows,
            score_scale_metadata=metadata,
        )


def score_kinase_library_motifs(
    **kwargs: Any,
) -> KinaseLibraryMotifScoringResult:
    """Convenience wrapper for ``KinaseLibraryMotifScorer().run(...)``."""

    return KinaseLibraryMotifScorer().run(**kwargs)


def _coerce_residue_classes(values: Sequence[object] | None) -> tuple[str, ...]:
    if values is None:
        return ()
    if isinstance(values, (str, bytes)) or hasattr(values, "value"):
        values = (values,)
    classes: list[str] = []
    for value in values:
        residue_class = normalize_kinase_library_residue_class(value)
        if residue_class not in classes:
            classes.append(residue_class)
    return tuple(classes)


def _coerce_matrix_lanes(
    matrices: object,
    *,
    matrix_residue_classes: Mapping[str, object] | None,
    default_residue_class: str | None,
) -> tuple[list[_MatrixLane], object | None, str | None]:
    if hasattr(matrices, "matrices") and not isinstance(matrices, Mapping):
        matrix_resource = cast(Any, matrices)
        entries = list(matrix_resource.matrices)
        resource_window = (
            matrix_resource.sequence_window
            if hasattr(matrix_resource, "sequence_window")
            else None
        )
        resource_score_scale = (
            matrix_resource.score_scale
            if hasattr(matrix_resource, "score_scale")
            else None
        )
        return (
            [
                _coerce_matrix_entry(entry, source_order=index)
                for index, entry in enumerate(entries)
            ],
            resource_window,
            None if resource_score_scale is None else str(resource_score_scale),
        )

    lanes: list[_MatrixLane] = []
    if isinstance(matrices, Mapping):
        source_order = 0
        for kinase, value in matrices.items():
            kinase_id = str(kinase)
            if isinstance(value, pd.DataFrame):
                residue_class = _resolve_mapping_matrix_residue_class(
                    kinase_id,
                    matrix_residue_classes=matrix_residue_classes,
                    default_residue_class=default_residue_class,
                )
                lanes.append(
                    _coerce_matrix_entry(
                        {
                            "kinase": kinase_id,
                            "residue_class": residue_class,
                            "score_table": value,
                        },
                        source_order=source_order,
                    )
                )
                source_order += 1
                continue
            if isinstance(value, Mapping):
                for residue_class, score_table in value.items():
                    lanes.append(
                        _coerce_matrix_entry(
                            {
                                "kinase": kinase_id,
                                "residue_class": residue_class,
                                "score_table": score_table,
                            },
                            source_order=source_order,
                        )
                    )
                    source_order += 1
                continue
            lanes.append(_coerce_matrix_entry(value, source_order=source_order))
            source_order += 1
        return lanes, None, None

    if isinstance(matrices, Sequence) and not isinstance(matrices, (str, bytes)):
        return (
            [
                _coerce_matrix_entry(entry, source_order=index)
                for index, entry in enumerate(matrices)
            ],
            None,
            None,
        )
    return ([_coerce_matrix_entry(matrices, source_order=0)], None, None)


def _resolve_mapping_matrix_residue_class(
    kinase: str,
    *,
    matrix_residue_classes: Mapping[str, object] | None,
    default_residue_class: str | None,
) -> str:
    if matrix_residue_classes is not None and kinase in matrix_residue_classes:
        return normalize_kinase_library_residue_class(matrix_residue_classes[kinase])
    if default_residue_class is not None:
        return default_residue_class
    raise ValueError(
        "matrix_residue_classes is required for mapping-style DataFrame matrices "
        "unless exactly one residue_class is requested"
    )


def _coerce_matrix_entry(entry: object, *, source_order: int) -> _MatrixLane:
    if isinstance(entry, Mapping):
        kinase = entry.get("kinase")
        residue_class = entry.get("residue_class")
        score_table = entry.get("score_table")
        kinase_family = entry.get("kinase_family")
        kinase_group = entry.get("kinase_group")
    else:
        kinase = getattr(entry, "kinase", None)
        residue_class = getattr(entry, "residue_class", None)
        score_table = getattr(entry, "score_table", None)
        kinase_family = getattr(entry, "kinase_family", None)
        kinase_group = getattr(entry, "kinase_group", None)

    if kinase is None or residue_class is None or score_table is None:
        raise ValueError(
            "each Kinase Library matrix must provide kinase, residue_class, "
            "and score_table"
        )
    matrix = KinaseLibraryMotifMatrix(
        kinase=str(kinase),
        residue_class=residue_class,
        score_table=score_table,
        kinase_family=None if kinase_family is None else str(kinase_family),
        kinase_group=None if kinase_group is None else str(kinase_group),
    )
    return _MatrixLane(
        kinase=matrix.kinase,
        residue_class=matrix.residue_class,
        score_table=matrix.score_table,
        kinase_family=matrix.kinase_family,
        kinase_group=matrix.kinase_group,
        source_order=int(source_order),
    )


def _classes_from_lanes(matrix_lanes: Sequence[_MatrixLane]) -> tuple[str, ...]:
    classes: list[str] = []
    for lane in matrix_lanes:
        if lane.residue_class not in classes:
            classes.append(lane.residue_class)
    return tuple(classes)


def _resolve_window_config(
    *,
    window: KinaseLibraryWindowConfig | object | None,
    upstream_residues: int | None,
    downstream_residues: int | None,
    flank_size: int | None,
    sequence_semantics: SequenceSemantics | None,
    matrix_lanes: Sequence[_MatrixLane],
    resource_window: object | None,
) -> KinaseLibraryWindowConfig:
    if flank_size is not None:
        if upstream_residues is not None or downstream_residues is not None:
            raise ValueError(
                "flank_size cannot be combined with upstream_residues or "
                "downstream_residues"
            )
        flank = int(flank_size)
        upstream_residues = flank
        downstream_residues = flank

    base_window = window if window is not None else resource_window
    if base_window is not None:
        window_obj = cast(Any, base_window)
        base_sequence_semantics = (
            window_obj.sequence_semantics
            if hasattr(window_obj, "sequence_semantics")
            else SEQUENCE_SEMANTICS_CENTRED_WINDOW
        )
        resolved = KinaseLibraryWindowConfig(
            upstream_residues=int(window_obj.upstream_residues),
            downstream_residues=int(window_obj.downstream_residues),
            sequence_semantics=sequence_semantics or base_sequence_semantics,
        )
        if upstream_residues is not None or downstream_residues is not None:
            return KinaseLibraryWindowConfig(
                upstream_residues=(
                    resolved.upstream_residues
                    if upstream_residues is None
                    else int(upstream_residues)
                ),
                downstream_residues=(
                    resolved.downstream_residues
                    if downstream_residues is None
                    else int(downstream_residues)
                ),
                sequence_semantics=resolved.sequence_semantics,
            )
        return resolved

    if upstream_residues is None or downstream_residues is None:
        inferred = _infer_window_from_matrices(matrix_lanes)
        if inferred is not None:
            upstream_residues = inferred.upstream_residues
            downstream_residues = inferred.downstream_residues
        else:
            upstream_residues = DEFAULT_MOTIF_FLANK_SIZE
            downstream_residues = DEFAULT_MOTIF_FLANK_SIZE

    return KinaseLibraryWindowConfig(
        upstream_residues=int(upstream_residues),
        downstream_residues=int(downstream_residues),
        sequence_semantics=sequence_semantics or SEQUENCE_SEMANTICS_CENTRED_WINDOW,
    )


def _infer_window_from_matrices(
    matrix_lanes: Sequence[_MatrixLane],
) -> KinaseLibraryWindowConfig | None:
    position_sets: list[tuple[int, ...]] = []
    for lane in matrix_lanes:
        try:
            positions = tuple(int(column) for column in lane.score_table.columns)
        except (TypeError, ValueError):
            continue
        if 0 not in positions:
            continue
        sorted_positions = tuple(sorted(positions))
        expected = tuple(range(min(sorted_positions), max(sorted_positions) + 1))
        if sorted_positions != expected:
            continue
        position_sets.append(sorted_positions)
    if not position_sets:
        return None
    first = position_sets[0]
    if any(positions != first for positions in position_sets):
        return None
    return KinaseLibraryWindowConfig(
        upstream_residues=abs(min(first)),
        downstream_residues=max(first),
        sequence_semantics=SEQUENCE_SEMANTICS_CENTRED_WINDOW,
    )


def _coerce_sequence_series(
    seqs: Mapping[str, object] | Sequence[object] | pd.Series,
    *,
    site_index: Sequence[object] | None,
) -> pd.Series:
    if isinstance(seqs, pd.Series):
        series = seqs.copy()
    elif isinstance(seqs, Mapping):
        series = pd.Series(dict(seqs), dtype=object)
    else:
        seq_list = list(seqs)
        if site_index is None:
            series = pd.Series(seq_list, dtype=object)
        else:
            if len(seq_list) != len(site_index):
                raise ValueError("site_index must have same length as site_sequences")
            series = pd.Series(seq_list, index=list(site_index), dtype=object)
    series.index = pd.Index(
        [str(item) for item in series.index], name=series.index.name
    )
    if site_index is not None:
        ordered_index = [str(item) for item in site_index]
        missing = [site for site in ordered_index if site not in series.index]
        if missing:
            raise ValueError(
                f"site_sequences missing entries for {', '.join(map(str, missing[:5]))}"
            )
        series = series.loc[ordered_index]
    return series.map(_normalize_sequence_value)


def _coerce_site_identity_series(
    values: Mapping[str, object] | pd.Series | None,
    *,
    site_index: Sequence[object],
) -> pd.Series:
    labels = [str(site_id) for site_id in site_index]
    if values is None:
        return pd.Series(labels, index=pd.Index(labels), dtype=object)
    if isinstance(values, pd.Series):
        series = values.copy()
    else:
        series = pd.Series(dict(values), dtype=object)
    series.index = pd.Index([str(item) for item in series.index])
    missing = [site for site in labels if site not in series.index]
    if missing:
        raise ValueError(
            f"site_identities missing entries for {', '.join(map(str, missing[:5]))}"
        )
    return series.loc[labels].map(
        lambda value: None if _is_missing_scalar(value) else str(value)
    )


def _validate_site_sequence(
    *,
    site_id: str,
    sequence_input: object,
    site_identity: object,
    window_config: KinaseLibraryWindowConfig,
    allowed_residue_classes: Sequence[str],
) -> _SiteValidation:
    expected_residue = _parse_expected_site_residue(
        None if _is_missing_scalar(site_identity) else str(site_identity)
    )
    if _is_missing_scalar(sequence_input) or not isinstance(sequence_input, str):
        return _SiteValidation(
            site_id=site_id,
            sequence=None,
            sequence_window=None,
            sequence_length=None,
            observed_central_residue=None,
            expected_central_residue=expected_residue,
            residue_class=None,
            status=KINASE_LIBRARY_SITE_STATUS_MISSING_SEQUENCE,
            reason="site_sequence is missing or blank",
        )

    sequence = str(sequence_input).strip().upper()
    if sequence == "":
        return _SiteValidation(
            site_id=site_id,
            sequence=None,
            sequence_window=None,
            sequence_length=None,
            observed_central_residue=None,
            expected_central_residue=expected_residue,
            residue_class=None,
            status=KINASE_LIBRARY_SITE_STATUS_MISSING_SEQUENCE,
            reason="site_sequence is missing or blank",
        )

    sequence_length = len(sequence)
    window, central_index, length_reason = _extract_kinase_library_window(
        sequence=sequence,
        window_config=window_config,
    )
    if window is None or central_index is None:
        return _SiteValidation(
            site_id=site_id,
            sequence=sequence,
            sequence_window=None,
            sequence_length=sequence_length,
            observed_central_residue=None,
            expected_central_residue=expected_residue,
            residue_class=None,
            status=KINASE_LIBRARY_SITE_STATUS_UNSUPPORTED_SEQUENCE_LENGTH,
            reason=length_reason,
        )

    observed_residue = sequence[central_index]
    if expected_residue is not None and observed_residue != expected_residue:
        return _SiteValidation(
            site_id=site_id,
            sequence=sequence,
            sequence_window=window,
            sequence_length=sequence_length,
            observed_central_residue=observed_residue,
            expected_central_residue=expected_residue,
            residue_class=_PHOSPHO_RESIDUE_CLASSES.get(observed_residue),
            status=KINASE_LIBRARY_SITE_STATUS_WRONG_CENTRAL_RESIDUE,
            reason=(
                f"expected centre residue '{expected_residue}', observed "
                f"'{observed_residue}'"
            ),
        )

    residue_class = _PHOSPHO_RESIDUE_CLASSES.get(observed_residue)
    if residue_class is None:
        return _SiteValidation(
            site_id=site_id,
            sequence=sequence,
            sequence_window=window,
            sequence_length=sequence_length,
            observed_central_residue=observed_residue,
            expected_central_residue=expected_residue,
            residue_class=None,
            status=KINASE_LIBRARY_SITE_STATUS_WRONG_CENTRAL_RESIDUE,
            reason=(
                f"centre residue '{observed_residue}' is not supported; expected "
                "S, T, or Y"
            ),
        )

    unsupported = sorted(
        {
            character
            for character in window
            if character not in _SUPPORTED_AMINO_ACID_SET
        }
    )
    if unsupported:
        return _SiteValidation(
            site_id=site_id,
            sequence=sequence,
            sequence_window=window,
            sequence_length=sequence_length,
            observed_central_residue=observed_residue,
            expected_central_residue=expected_residue,
            residue_class=residue_class,
            status=KINASE_LIBRARY_SITE_STATUS_UNSUPPORTED_RESIDUE,
            reason=f"unsupported residue character(s): {', '.join(unsupported)}",
        )

    if residue_class not in allowed_residue_classes:
        return _SiteValidation(
            site_id=site_id,
            sequence=sequence,
            sequence_window=window,
            sequence_length=sequence_length,
            observed_central_residue=observed_residue,
            expected_central_residue=expected_residue,
            residue_class=residue_class,
            status=KINASE_LIBRARY_SITE_STATUS_WRONG_RESIDUE_CLASS,
            reason=(
                f"site residue class '{residue_class}' is outside configured "
                f"classes: {', '.join(allowed_residue_classes)}"
            ),
        )

    return _SiteValidation(
        site_id=site_id,
        sequence=sequence,
        sequence_window=window,
        sequence_length=sequence_length,
        observed_central_residue=observed_residue,
        expected_central_residue=expected_residue,
        residue_class=residue_class,
        status=KINASE_LIBRARY_SITE_STATUS_VALID_UNSCORED_SITE,
        reason=None,
    )


def _extract_kinase_library_window(
    *,
    sequence: str,
    window_config: KinaseLibraryWindowConfig,
) -> tuple[str | None, int | None, str | None]:
    if window_config.sequence_semantics == SEQUENCE_SEMANTICS_CENTRED_WINDOW:
        if len(sequence) != window_config.window_size:
            return (
                None,
                None,
                f"sequence length {len(sequence)} does not match expected "
                f"window length {window_config.window_size}",
            )
        return sequence, window_config.centre_index, None

    if window_config.sequence_semantics != SEQUENCE_SEMANTICS_CENTRED_SEQUENCE:
        raise ValueError(
            "sequence_semantics must be 'centred_window' or 'centred_sequence'"
        )
    if len(sequence) < window_config.window_size:
        return (
            None,
            None,
            f"sequence length {len(sequence)} is shorter than expected window "
            f"length {window_config.window_size}",
        )
    if len(sequence) % 2 == 0:
        return None, None, "centred sequences must have odd length"
    central_index = len(sequence) // 2
    start = central_index - window_config.upstream_residues
    stop = central_index + window_config.downstream_residues + 1
    if start < 0 or stop > len(sequence):
        return (
            None,
            None,
            "configured upstream/downstream residues exceed supplied sequence",
        )
    return sequence[start:stop], central_index, None


def _parse_expected_site_residue(site_identity: str | None) -> str | None:
    if site_identity is None:
        return None
    match = _SITE_IDENTITY_PATTERN.fullmatch(site_identity)
    if match is None:
        return None
    site_token = match.group("site").strip().upper()
    if not site_token:
        return None
    residue = site_token[0]
    return residue if residue in _SUPPORTED_AMINO_ACID_SET else None


def _validate_matrix_lanes(
    *,
    matrix_lanes: Sequence[_MatrixLane],
    window_config: KinaseLibraryWindowConfig,
    allowed_residue_classes: Sequence[str],
) -> tuple[list[_MatrixLane], pd.DataFrame]:
    rows: list[dict[str, object]] = []
    valid_lanes: list[_MatrixLane] = []
    seen_keys: set[tuple[str, str]] = set()
    for lane in matrix_lanes:
        row = _base_kinase_diagnostic_row(lane, window_config)
        if lane.residue_class not in allowed_residue_classes:
            row.update(
                {
                    "status": KINASE_LIBRARY_MATRIX_STATUS_FILTERED_RESIDUE_CLASS,
                    "reason": (
                        f"matrix residue class '{lane.residue_class}' is outside "
                        "configured residue classes"
                    ),
                }
            )
            rows.append(row)
            continue

        normalized, error = _normalize_score_table(lane.score_table)
        if error is not None or normalized is None:
            row.update(
                {
                    "status": KINASE_LIBRARY_MATRIX_STATUS_INVALID,
                    "reason": error,
                }
            )
            rows.append(row)
            continue

        observed_positions = tuple(int(position) for position in normalized.columns)
        expected_positions = window_config.positions
        if set(observed_positions) != set(expected_positions):
            missing = [
                position
                for position in expected_positions
                if position not in observed_positions
            ]
            unexpected = [
                position
                for position in observed_positions
                if position not in expected_positions
            ]
            row.update(
                {
                    "status": KINASE_LIBRARY_MATRIX_STATUS_UNSUPPORTED_WINDOW,
                    "reason": _format_position_mismatch_reason(
                        missing=missing,
                        unexpected=unexpected,
                    ),
                    "matrix_positions": "|".join(map(str, observed_positions)),
                }
            )
            rows.append(row)
            continue

        key = (lane.kinase, lane.residue_class)
        if key in seen_keys:
            row.update(
                {
                    "status": KINASE_LIBRARY_MATRIX_STATUS_DUPLICATE,
                    "reason": "duplicate kinase/residue_class matrix lane",
                    "matrix_positions": "|".join(map(str, observed_positions)),
                }
            )
            rows.append(row)
            continue
        seen_keys.add(key)

        normalized = normalized.loc[:, list(expected_positions)]
        valid_lane = _MatrixLane(
            kinase=lane.kinase,
            residue_class=lane.residue_class,
            score_table=normalized,
            kinase_family=lane.kinase_family,
            kinase_group=lane.kinase_group,
            source_order=lane.source_order,
        )
        row.update(
            {
                "status": KINASE_LIBRARY_MATRIX_STATUS_VALID,
                "reason": None,
                "matrix_positions": "|".join(map(str, expected_positions)),
                "matrix_amino_acid_count": int(normalized.shape[0]),
            }
        )
        rows.append(row)
        valid_lanes.append(valid_lane)

    diagnostics = pd.DataFrame(rows)
    if diagnostics.empty:
        diagnostics = pd.DataFrame(
            columns=[
                "kinase",
                "residue_class",
                "status",
                "reason",
                "kinase_family",
                "kinase_group",
                "matrix_positions",
                "expected_positions",
                "matrix_amino_acid_count",
                "sites_eligible_for_residue_class",
                "sites_scored",
                "missing_residue_score_cells",
            ]
        )
        diagnostics.index = pd.MultiIndex.from_arrays(
            [[], []],
            names=["kinase", "residue_class"],
        )
        return valid_lanes, diagnostics

    diagnostics.index = pd.MultiIndex.from_frame(
        diagnostics.loc[:, ["kinase", "residue_class"]],
        names=["kinase", "residue_class"],
    )
    return valid_lanes, diagnostics


def _base_kinase_diagnostic_row(
    lane: _MatrixLane,
    window_config: KinaseLibraryWindowConfig,
) -> dict[str, object]:
    return {
        "kinase": lane.kinase,
        "residue_class": lane.residue_class,
        "status": KINASE_LIBRARY_MATRIX_STATUS_INVALID,
        "reason": None,
        "kinase_family": lane.kinase_family,
        "kinase_group": lane.kinase_group,
        "matrix_positions": None,
        "expected_positions": "|".join(map(str, window_config.positions)),
        "matrix_amino_acid_count": 0,
        "sites_eligible_for_residue_class": 0,
        "sites_scored": 0,
        "missing_residue_score_cells": 0,
    }


def _normalize_score_table(
    score_table: pd.DataFrame,
) -> tuple[pd.DataFrame | None, str | None]:
    if not isinstance(score_table, pd.DataFrame):
        return None, "score_table must be a pandas DataFrame"
    if score_table.empty:
        return None, "score_table must be non-empty"
    amino_acids = [str(value).strip().upper() for value in score_table.index]
    if len(set(amino_acids)) != len(amino_acids):
        return None, "score_table contains duplicate amino-acid rows"
    invalid_amino_acids = [
        value
        for value in amino_acids
        if len(value) != 1 or value not in _SUPPORTED_AMINO_ACID_SET
    ]
    if invalid_amino_acids:
        return None, "score_table contains unsupported amino-acid rows"
    try:
        positions = [int(column) for column in score_table.columns]
    except (TypeError, ValueError):
        return None, "score_table columns must be integer relative positions"
    if len(set(positions)) != len(positions):
        return None, "score_table contains duplicate position columns"
    try:
        normalized = score_table.copy(deep=True).astype(float)
    except ValueError:
        return None, "score_table must contain numeric score values"
    if normalized.isna().to_numpy().any():
        return None, "score_table contains missing score values"
    values = normalized.to_numpy(dtype=float, copy=False)
    if not np.isfinite(values).all():
        return None, "score_table contains non-finite score values"
    normalized.index = pd.Index(amino_acids, name="amino_acid")
    normalized.columns = pd.Index(positions, name="position")
    return normalized, None


def _format_position_mismatch_reason(
    *,
    missing: Sequence[int],
    unexpected: Sequence[int],
) -> str:
    parts: list[str] = []
    if missing:
        parts.append(f"missing required positions: {', '.join(map(str, missing))}")
    if unexpected:
        parts.append(f"unexpected positions: {', '.join(map(str, unexpected))}")
    return "; ".join(parts) if parts else "matrix positions do not match window"


def _stable_kinase_columns(valid_lanes: Sequence[_MatrixLane]) -> list[str]:
    columns: list[str] = []
    for lane in sorted(valid_lanes, key=lambda item: item.source_order):
        if lane.kinase not in columns:
            columns.append(lane.kinase)
    return columns


def _score_sequence_window(
    *,
    sequence_window: str,
    score_table: pd.DataFrame,
    positions: Sequence[int],
) -> float | None:
    score = 0.0
    for amino_acid, position in zip(sequence_window, positions, strict=True):
        if amino_acid not in score_table.index:
            return None
        score += float(cast(Any, score_table.at[amino_acid, position]))
    if not math.isfinite(score):
        return None
    return score


def _build_site_diagnostics(
    *,
    site_validations: Sequence[_SiteValidation],
    raw_scores: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for site in site_validations:
        scored_kinase_count = 0
        status = site.status
        reason = site.reason
        if site.site_id in raw_scores.index:
            site_scores = cast(pd.Series, raw_scores.loc[site.site_id])
            scored_kinase_count = int(site_scores.notna().sum())
        if site.status == KINASE_LIBRARY_SITE_STATUS_VALID_UNSCORED_SITE:
            if scored_kinase_count > 0:
                status = KINASE_LIBRARY_SITE_STATUS_VALID_SCORED_SITE
                reason = None
            else:
                status = KINASE_LIBRARY_SITE_STATUS_VALID_UNSCORED_SITE
                reason = (
                    "site passed sequence validation but no kinase matrix scored it"
                )
        rows.append(
            {
                "site_id": site.site_id,
                "status": status,
                "reason": reason,
                "sequence": site.sequence,
                "sequence_window": site.sequence_window,
                "sequence_length": site.sequence_length,
                "observed_central_residue": site.observed_central_residue,
                "expected_central_residue": site.expected_central_residue,
                "residue_class": site.residue_class,
                "scored_kinase_count": scored_kinase_count,
                "excluded_from_scoring": status
                != KINASE_LIBRARY_SITE_STATUS_VALID_SCORED_SITE,
            }
        )
    diagnostics = pd.DataFrame(rows)
    diagnostics.index = pd.Index([row["site_id"] for row in rows], name="site_id")
    return diagnostics


def _finish_kinase_diagnostics(
    *,
    kinase_diagnostics: pd.DataFrame,
    site_validations: Sequence[_SiteValidation],
    raw_scores: pd.DataFrame,
    missing_residue_score_cells: Mapping[tuple[str, str], int],
) -> pd.DataFrame:
    diagnostics = kinase_diagnostics.copy(deep=True)
    eligible_by_class: dict[str, int] = defaultdict(int)
    for site in site_validations:
        if (
            site.status == KINASE_LIBRARY_SITE_STATUS_VALID_UNSCORED_SITE
            and site.residue_class is not None
        ):
            eligible_by_class[site.residue_class] += 1

    for index in diagnostics.index:
        kinase, residue_class = index
        diagnostics.at[index, "sites_eligible_for_residue_class"] = int(
            eligible_by_class[str(residue_class)]
        )
        diagnostics.at[index, "missing_residue_score_cells"] = int(
            missing_residue_score_cells.get((str(kinase), str(residue_class)), 0)
        )
        if (
            diagnostics.at[index, "status"] == KINASE_LIBRARY_MATRIX_STATUS_VALID
            and str(kinase) in raw_scores.columns
        ):
            class_site_ids = [
                site.site_id
                for site in site_validations
                if (
                    site.status == KINASE_LIBRARY_SITE_STATUS_VALID_UNSCORED_SITE
                    and site.residue_class == residue_class
                    and site.site_id in raw_scores.index
                )
            ]
            if class_site_ids:
                diagnostics.at[index, "sites_scored"] = int(
                    raw_scores.loc[class_site_ids, str(kinase)].notna().sum()
                )
            else:
                diagnostics.at[index, "sites_scored"] = 0
    return diagnostics


def _coerce_reference_distributions(
    reference_distributions: Mapping[object, object] | None,
) -> dict[tuple[str, str | None], np.ndarray]:
    if reference_distributions is None:
        return {}
    lookup: dict[tuple[str, str | None], np.ndarray] = {}
    for key, value in reference_distributions.items():
        if isinstance(key, tuple) and len(key) == 2:
            kinase = str(key[0])
            residue_class = normalize_kinase_library_residue_class(key[1])
            lookup[(kinase, residue_class)] = _coerce_distribution_values(value)
            continue
        if isinstance(value, Mapping) and not isinstance(value, pd.Series):
            for residue_class_key, distribution in value.items():
                lookup[
                    (
                        str(key),
                        normalize_kinase_library_residue_class(residue_class_key),
                    )
                ] = _coerce_distribution_values(distribution)
            continue
        lookup[(str(key), None)] = _coerce_distribution_values(value)
    return lookup


def _coerce_distribution_values(value: object) -> np.ndarray:
    if isinstance(value, pd.Series):
        values = value.to_numpy(dtype=float, copy=True)
    else:
        values = np.asarray(list(value), dtype=float)  # type: ignore[arg-type]
    return values[np.isfinite(values)]


def _score_reference_distributions(
    *,
    raw_scores: pd.DataFrame,
    site_diagnostics: pd.DataFrame,
    distribution_lookup: Mapping[tuple[str, str | None], np.ndarray],
    higher_is_better: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    percentiles = pd.DataFrame(
        np.nan,
        index=raw_scores.index.copy(),
        columns=raw_scores.columns.copy(),
        dtype=float,
    )
    ranks = percentiles.copy(deep=True)
    for site_id in raw_scores.index.astype(str):
        residue_class = site_diagnostics.at[site_id, "residue_class"]
        if residue_class is None or pd.isna(residue_class):
            continue
        for kinase in raw_scores.columns.astype(str):
            score = raw_scores.at[site_id, kinase]
            if pd.isna(score):
                continue
            distribution = _distribution_for(
                distribution_lookup=distribution_lookup,
                kinase=kinase,
                residue_class=str(residue_class),
            )
            if distribution.size == 0:
                continue
            score_value = float(cast(Any, score))
            if higher_is_better:
                percentiles.at[site_id, kinase] = (
                    float(np.count_nonzero(distribution <= score_value))
                    / float(distribution.size)
                    * 100.0
                )
                ranks.at[site_id, kinase] = float(
                    1 + np.count_nonzero(distribution > score_value)
                )
            else:
                percentiles.at[site_id, kinase] = (
                    float(np.count_nonzero(distribution >= score_value))
                    / float(distribution.size)
                    * 100.0
                )
                ranks.at[site_id, kinase] = float(
                    1 + np.count_nonzero(distribution < score_value)
                )
    return percentiles, ranks


def _distribution_for(
    *,
    distribution_lookup: Mapping[tuple[str, str | None], np.ndarray],
    kinase: str,
    residue_class: str,
) -> np.ndarray:
    class_specific = distribution_lookup.get((kinase, residue_class))
    if class_specific is not None:
        return class_specific
    return distribution_lookup.get((kinase, None), np.asarray([], dtype=float))


def _percentile_method_label(*, higher_is_better: bool) -> str:
    if higher_is_better:
        return _PERCENTILE_METHOD_HIGHER_IS_BETTER
    return _PERCENTILE_METHOD_LOWER_IS_BETTER


def _rank_method_label(*, higher_is_better: bool) -> str:
    if higher_is_better:
        return _RANK_METHOD_HIGHER_IS_BETTER
    return _RANK_METHOD_LOWER_IS_BETTER


__all__ = [
    "KinaseLibraryMotifScorer",
    "score_kinase_library_motifs",
]

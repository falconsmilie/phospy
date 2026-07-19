from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Literal

import numpy as np
import pandas as pd

from phospy.frames.ownership import (
    export_dataframe,
    export_optional_dataframe,
    export_series,
    own_dataframe,
    own_optional_dataframe,
    own_series,
)
from phospy.science.prediction.sequence_validation import (
    SUPPORTED_AMINO_ACIDS,
    SequenceValidationResult,
    SequenceValidationStatus,
)

AMINO_ACIDS: tuple[str, ...] = SUPPORTED_AMINO_ACIDS
DEFAULT_MOTIF_FLANK_SIZE = 7
SEQUENCE_SEMANTICS_CENTRED_WINDOW = "centred_window"
SEQUENCE_SEMANTICS_CENTRED_SEQUENCE = "centred_sequence"
SequenceSemantics = Literal["centred_window", "centred_sequence"]
KINASE_LIBRARY_RESIDUE_CLASS_SER_THR = "ser_thr"
KINASE_LIBRARY_RESIDUE_CLASS_TYR = "tyr"
KinaseLibraryResidueClassName = Literal["ser_thr", "tyr"]
KINASE_LIBRARY_RESIDUE_CLASSES: tuple[str, ...] = (
    KINASE_LIBRARY_RESIDUE_CLASS_SER_THR,
    KINASE_LIBRARY_RESIDUE_CLASS_TYR,
)
KINASE_LIBRARY_SITE_STATUS_MISSING_SEQUENCE = "missing_sequence"
KINASE_LIBRARY_SITE_STATUS_WRONG_CENTRAL_RESIDUE = "wrong_central_residue"
KINASE_LIBRARY_SITE_STATUS_WRONG_RESIDUE_CLASS = "wrong_residue_class"
KINASE_LIBRARY_SITE_STATUS_UNSUPPORTED_SEQUENCE_LENGTH = "unsupported_sequence_length"
KINASE_LIBRARY_SITE_STATUS_UNSUPPORTED_RESIDUE = "unsupported_residue_character"
KINASE_LIBRARY_SITE_STATUS_VALID_SCORED_SITE = "valid_scored_site"
KINASE_LIBRARY_SITE_STATUS_VALID_UNSCORED_SITE = "valid_unscored_site"
KinaseLibrarySiteStatus = Literal[
    "missing_sequence",
    "wrong_central_residue",
    "wrong_residue_class",
    "unsupported_sequence_length",
    "unsupported_residue_character",
    "valid_scored_site",
    "valid_unscored_site",
]
KINASE_LIBRARY_MATRIX_STATUS_VALID = "valid"
KINASE_LIBRARY_MATRIX_STATUS_FILTERED_RESIDUE_CLASS = "filtered_residue_class"
KINASE_LIBRARY_MATRIX_STATUS_UNSUPPORTED_WINDOW = "unsupported_matrix_window"
KINASE_LIBRARY_MATRIX_STATUS_DUPLICATE = "duplicate_matrix"
KINASE_LIBRARY_MATRIX_STATUS_INVALID = "invalid_matrix"
KinaseLibraryMatrixStatus = Literal[
    "valid",
    "filtered_residue_class",
    "unsupported_matrix_window",
    "duplicate_matrix",
    "invalid_matrix",
]
_MOTIF_LIBRARY_VALIDATION_ATTR = "motif_library_validation"
_ASCII_LOOKUP_SIZE = 256
_INVALID_AMINO_ACID_INDEX = -1
_AMINO_ACID_INDEX_LOOKUP = np.full(
    _ASCII_LOOKUP_SIZE,
    _INVALID_AMINO_ACID_INDEX,
    dtype=np.int16,
)
for _row_index, _amino_acid in enumerate(AMINO_ACIDS):
    _AMINO_ACID_INDEX_LOOKUP[ord(_amino_acid)] = _row_index


@dataclass(frozen=True, slots=True, init=False)
class MotifScoringResult:
    """Motif score matrices and window metadata for one scoring run."""

    sequence_validation: SequenceValidationResult
    library_validation: MotifLibraryValidationResult | None = None
    _motif_scores: pd.DataFrame = field(init=False, repr=False)
    _motif_sizes: pd.Series = field(init=False, repr=False)
    _sequence_windows: pd.Series = field(init=False, repr=False)

    def __init__(
        self,
        motif_scores: pd.DataFrame,
        motif_sizes: pd.Series,
        sequence_windows: pd.Series,
        sequence_validation: SequenceValidationResult,
        library_validation: MotifLibraryValidationResult | None = None,
    ) -> None:
        object.__setattr__(
            self,
            "_motif_scores",
            own_dataframe(
                motif_scores,
                field_name="motif_scoring_result.motif_scores",
            ),
        )
        object.__setattr__(
            self,
            "_motif_sizes",
            own_series(
                motif_sizes,
                field_name="motif_scoring_result.motif_sizes",
            ),
        )
        object.__setattr__(
            self,
            "_sequence_windows",
            own_series(
                sequence_windows,
                field_name="motif_scoring_result.sequence_windows",
            ),
        )
        object.__setattr__(self, "sequence_validation", sequence_validation)
        object.__setattr__(self, "library_validation", library_validation)

    @property
    def motif_scores(self) -> pd.DataFrame:
        return export_dataframe(self._motif_scores)

    @property
    def motif_sizes(self) -> pd.Series:
        return export_series(self._motif_sizes)

    @property
    def sequence_windows(self) -> pd.Series:
        return export_series(self._sequence_windows)


@dataclass(frozen=True, slots=True)
class KinaseLibraryWindowConfig:
    """Sequence-window contract for Kinase Library-style motif scoring."""

    upstream_residues: int = DEFAULT_MOTIF_FLANK_SIZE
    downstream_residues: int = DEFAULT_MOTIF_FLANK_SIZE
    sequence_semantics: SequenceSemantics = SEQUENCE_SEMANTICS_CENTRED_WINDOW

    def __post_init__(self) -> None:
        upstream = int(self.upstream_residues)
        downstream = int(self.downstream_residues)
        if upstream < 0:
            raise ValueError("upstream_residues must be >= 0")
        if downstream < 0:
            raise ValueError("downstream_residues must be >= 0")
        if self.sequence_semantics not in {
            SEQUENCE_SEMANTICS_CENTRED_WINDOW,
            SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
        }:
            raise ValueError(
                "sequence_semantics must be 'centred_window' or 'centred_sequence'"
            )
        object.__setattr__(self, "upstream_residues", upstream)
        object.__setattr__(self, "downstream_residues", downstream)

    @property
    def window_size(self) -> int:
        return int(self.upstream_residues) + 1 + int(self.downstream_residues)

    @property
    def centre_index(self) -> int:
        return int(self.upstream_residues)

    @property
    def positions(self) -> tuple[int, ...]:
        return tuple(
            range(-int(self.upstream_residues), int(self.downstream_residues) + 1)
        )

    def to_payload(self) -> dict[str, object]:
        return {
            "upstream_residues": int(self.upstream_residues),
            "downstream_residues": int(self.downstream_residues),
            "window_size": int(self.window_size),
            "central_residue_required": True,
            "sequence_semantics": self.sequence_semantics,
        }


@dataclass(frozen=True, slots=True)
class KinaseLibraryMotifMatrix:
    """Lightweight Kinase Library-style matrix accepted by the pure scorer."""

    kinase: str
    residue_class: KinaseLibraryResidueClassName | str
    score_table: pd.DataFrame
    kinase_family: str | None = None
    kinase_group: str | None = None

    def __post_init__(self) -> None:
        kinase = str(self.kinase).strip()
        if kinase == "":
            raise ValueError("kinase must be a non-empty string")
        residue_class = normalize_kinase_library_residue_class(self.residue_class)
        if not isinstance(self.score_table, pd.DataFrame):
            raise TypeError("score_table must be a pandas DataFrame")
        object.__setattr__(self, "kinase", kinase)
        object.__setattr__(self, "residue_class", residue_class)
        object.__setattr__(self, "score_table", self.score_table.copy(deep=True))
        object.__setattr__(
            self,
            "kinase_family",
            None if self.kinase_family is None else str(self.kinase_family),
        )
        object.__setattr__(
            self,
            "kinase_group",
            None if self.kinase_group is None else str(self.kinase_group),
        )


@dataclass(frozen=True, slots=True)
class KinaseLibraryScoreScaleMetadata:
    """Machine-readable score-scale metadata for one motif-scoring run."""

    score_scale: str
    raw_score_formula: str
    higher_is_better: bool
    percentile_method: str | None
    rank_method: str | None
    sequence_window: Mapping[str, object]
    residue_classes: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "score_scale", str(self.score_scale))
        object.__setattr__(self, "raw_score_formula", str(self.raw_score_formula))
        object.__setattr__(self, "higher_is_better", bool(self.higher_is_better))
        object.__setattr__(
            self,
            "percentile_method",
            None if self.percentile_method is None else str(self.percentile_method),
        )
        object.__setattr__(
            self,
            "rank_method",
            None if self.rank_method is None else str(self.rank_method),
        )
        object.__setattr__(
            self,
            "sequence_window",
            MappingProxyType(
                {str(key): value for key, value in self.sequence_window.items()}
            ),
        )
        object.__setattr__(
            self,
            "residue_classes",
            tuple(
                normalize_kinase_library_residue_class(item)
                for item in self.residue_classes
            ),
        )

    def to_payload(self) -> dict[str, object]:
        return {
            "score_scale": self.score_scale,
            "raw_score_formula": self.raw_score_formula,
            "higher_is_better": self.higher_is_better,
            "percentile_method": self.percentile_method,
            "rank_method": self.rank_method,
            "sequence_window": dict(self.sequence_window),
            "residue_classes": list(self.residue_classes),
        }


@dataclass(frozen=True, slots=True, init=False)
class KinaseLibraryMotifScoringResult:
    """Outputs from the pure Kinase Library-style motif scoring engine."""

    score_scale_metadata: KinaseLibraryScoreScaleMetadata
    _raw_scores: pd.DataFrame = field(init=False, repr=False)
    _percentile_ranks: pd.DataFrame | None = field(init=False, repr=False)
    _reference_ranks: pd.DataFrame | None = field(init=False, repr=False)
    _site_diagnostics: pd.DataFrame = field(init=False, repr=False)
    _kinase_diagnostics: pd.DataFrame = field(init=False, repr=False)
    _sequence_windows: pd.Series = field(init=False, repr=False)

    def __init__(
        self,
        raw_scores: pd.DataFrame,
        percentile_ranks: pd.DataFrame | None,
        reference_ranks: pd.DataFrame | None,
        site_diagnostics: pd.DataFrame,
        kinase_diagnostics: pd.DataFrame,
        sequence_windows: pd.Series,
        score_scale_metadata: KinaseLibraryScoreScaleMetadata,
    ) -> None:
        object.__setattr__(
            self,
            "_raw_scores",
            own_dataframe(
                raw_scores,
                field_name="kinase_library_motif_scoring_result.raw_scores",
            ),
        )
        object.__setattr__(
            self,
            "_percentile_ranks",
            own_optional_dataframe(
                percentile_ranks,
                field_name=("kinase_library_motif_scoring_result.percentile_ranks"),
            ),
        )
        object.__setattr__(
            self,
            "_reference_ranks",
            own_optional_dataframe(
                reference_ranks,
                field_name="kinase_library_motif_scoring_result.reference_ranks",
            ),
        )
        object.__setattr__(
            self,
            "_site_diagnostics",
            own_dataframe(
                site_diagnostics,
                field_name="kinase_library_motif_scoring_result.site_diagnostics",
            ),
        )
        object.__setattr__(
            self,
            "_kinase_diagnostics",
            own_dataframe(
                kinase_diagnostics,
                field_name="kinase_library_motif_scoring_result.kinase_diagnostics",
            ),
        )
        object.__setattr__(
            self,
            "_sequence_windows",
            own_series(
                sequence_windows,
                field_name="kinase_library_motif_scoring_result.sequence_windows",
            ),
        )
        object.__setattr__(self, "score_scale_metadata", score_scale_metadata)

    @property
    def raw_scores(self) -> pd.DataFrame:
        return export_dataframe(self._raw_scores)

    @property
    def percentile_ranks(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._percentile_ranks)

    @property
    def reference_ranks(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._reference_ranks)

    @property
    def site_diagnostics(self) -> pd.DataFrame:
        return export_dataframe(self._site_diagnostics)

    @property
    def kinase_diagnostics(self) -> pd.DataFrame:
        return export_dataframe(self._kinase_diagnostics)

    @property
    def sequence_windows(self) -> pd.Series:
        return export_series(self._sequence_windows)

    def raw_scores_dataframe(self) -> pd.DataFrame:
        """Return raw motif scores isolated from this result."""

        return export_dataframe(self._raw_scores)

    def percentile_ranks_dataframe(self) -> pd.DataFrame | None:
        """Return optional percentile ranks isolated from this result."""

        return export_optional_dataframe(self._percentile_ranks)

    def reference_ranks_dataframe(self) -> pd.DataFrame | None:
        """Return optional reference ranks isolated from this result."""

        return export_optional_dataframe(self._reference_ranks)

    def site_diagnostics_dataframe(self) -> pd.DataFrame:
        """Return site diagnostics isolated from this result."""

        return export_dataframe(self._site_diagnostics)

    def kinase_diagnostics_dataframe(self) -> pd.DataFrame:
        """Return kinase diagnostics isolated from this result."""

        return export_dataframe(self._kinase_diagnostics)

    def sequence_windows_series(self) -> pd.Series:
        """Return sequence windows isolated from this result."""

        return export_series(self._sequence_windows)

    def score_scale_payload(self) -> dict[str, object]:
        return self.score_scale_metadata.to_payload()


@dataclass(frozen=True, slots=True)
class MotifLibraryValidationRow:
    """Per-reference motif-library validation outcome."""

    reference_id: str
    site_id: str | None
    kinase: str
    sequence: str | None
    status: SequenceValidationStatus
    reason: str | None
    expected_centre_residue: str | None
    observed_centre_residue: str | None
    sequence_length: int | None


@dataclass(frozen=True, slots=True)
class MotifLibraryValidationResult:
    """Structured diagnostics for motif-library/reference sequence validation."""

    total_reference_sequences: int
    accepted_reference_sequences: int
    excluded_reference_sequences: int
    missing_sequences: int
    short_sequences: int
    off_centre_sequences: int
    site_residue_mismatches: int
    invalid_site_ids: int
    non_phospho_centre_residues: int
    unsupported_residue_characters: int
    sequences_excluded_from_motif_profile_construction: int
    expected_window_size: int
    supported_amino_acids: tuple[str, ...]
    accepted_window_length_policy: str
    unsupported_residue_policy: str
    excluded_reference_ids: tuple[str, ...]
    rows: tuple[MotifLibraryValidationRow, ...]

    def summary(self) -> dict[str, object]:
        """Return compact diagnostics for run-level reporting/provenance."""

        return {
            "reference_sequences_provided": self.total_reference_sequences,
            "reference_sequences_accepted": self.accepted_reference_sequences,
            "reference_sequences_excluded": self.excluded_reference_sequences,
            "excluded_missing_sequence": self.missing_sequences,
            "excluded_short_window": self.short_sequences,
            "excluded_unsupported_residue": self.unsupported_residue_characters,
            "excluded_off_centre_residue": self.off_centre_sequences,
            "excluded_site_residue_mismatch": self.site_residue_mismatches,
            "excluded_invalid_site_id": self.invalid_site_ids,
            "excluded_non_phospho_centre_residue": (self.non_phospho_centre_residues),
            "sequences_excluded_from_motif_profile_construction": (
                self.sequences_excluded_from_motif_profile_construction
            ),
            "expected_window_size": self.expected_window_size,
            "accepted_window_length_policy": self.accepted_window_length_policy,
            "unsupported_residue_policy": self.unsupported_residue_policy,
            "supported_amino_acids": list(self.supported_amino_acids),
        }

    def to_frame(self) -> pd.DataFrame:
        """Return row-level library validation provenance as a compact table."""

        columns = [
            "reference_id",
            "site_id",
            "kinase",
            "sequence",
            "status",
            "reason",
            "observed_centre_residue",
            "expected_centre_residue",
            "sequence_length",
        ]
        return pd.DataFrame(
            [
                {
                    "reference_id": row.reference_id,
                    "site_id": row.site_id,
                    "kinase": row.kinase,
                    "sequence": row.sequence,
                    "status": row.status,
                    "reason": row.reason,
                    "observed_centre_residue": row.observed_centre_residue,
                    "expected_centre_residue": row.expected_centre_residue,
                    "sequence_length": row.sequence_length,
                }
                for row in self.rows
            ],
            columns=pd.Index(columns),
        )


@dataclass(frozen=True, slots=True)
class _LibraryCandidate:
    reference_id: str
    site_id: str | None
    kinase: str
    sequence_input: object


@dataclass(frozen=True, slots=True)
class ExplicitMotifSequence:
    """Structured explicit motif-library sequence metadata."""

    reference_id: str
    site_id: str | None
    kinase: str
    sequence: object


def normalize_kinase_library_residue_class(value: object) -> str:
    """Return the stable scorer residue-class token for supported inputs."""

    text = getattr(value, "value", value)
    normalized = str(text).strip().lower().replace("-", "_").replace("/", "_")
    normalized = "_".join(part for part in normalized.split("_") if part)
    if normalized in {"ser_thr", "st", "s_t", "serine_threonine"}:
        return KINASE_LIBRARY_RESIDUE_CLASS_SER_THR
    if normalized in {"tyr", "y", "tyrosine"}:
        return KINASE_LIBRARY_RESIDUE_CLASS_TYR
    raise ValueError("residue_class must be 'ser_thr' or 'tyr'")


def _normalize_sequence_value(value: object) -> object:
    if value is None:
        return np.nan
    if bool(pd.Series((value,), dtype="object").isna().iat[0]):
        return np.nan
    if not isinstance(value, str):
        return value
    sequence = value.strip().upper()
    if sequence == "":
        return np.nan
    return sequence


def _is_missing_scalar(value: object) -> bool:
    if value is None:
        return True
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _extract_sequence_window(value: object, flank_size: int | None) -> object:
    if _is_missing_scalar(value):
        return np.nan
    sequence = str(value).upper()
    if sequence == "":
        return np.nan
    if flank_size is None:
        return sequence
    window_size = (2 * flank_size) + 1
    if len(sequence) <= window_size:
        return sequence
    mid = len(sequence) // 2
    start = mid - flank_size
    stop = mid + flank_size + 1
    return sequence[start:stop]

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Protocol, runtime_checkable

import numpy as np
import pandas as pd

from .constants import (
    BUNDLED_REFERENCE_ALIASES,
    BUNDLED_REFERENCE_AUTO,
    BUNDLED_REFERENCE_DEFAULTS,
    BUNDLED_REFERENCE_PROVIDER_NAME,
    BUNDLED_REFERENCE_SOURCE,
    BUNDLED_REFERENCE_SPECIES_ALIASES,
    BUNDLED_REFERENCE_VERSION,
)
from .types import KinaseMotifSequenceMap, KinaseSubstrateMap
from .validation.collections import normalize_sequence_mapping
from .validation.errors import (
    InputCompatibilityError,
    TableSchemaError,
)
from .validation.scalars import validate_non_negative_int, validate_positive_int

AMINO_ACIDS: tuple[str, ...] = (
    "A",
    "R",
    "N",
    "D",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
)

_ASCII_LOOKUP_SIZE = 256
_INVALID_AMINO_ACID_INDEX = -1
_AMINO_ACID_INDEX_LOOKUP = np.full(
    _ASCII_LOOKUP_SIZE,
    _INVALID_AMINO_ACID_INDEX,
    dtype=np.int16,
)
for _row_index, _amino_acid in enumerate(AMINO_ACIDS):
    _AMINO_ACID_INDEX_LOOKUP[ord(_amino_acid)] = _row_index

_ALLOWED_SEQUENCE_CHARACTERS: frozenset[str] = frozenset(AMINO_ACIDS) | frozenset({"_"})


@dataclass(frozen=True, slots=True)
class ReferenceBundleSourceMetadata:
    """Source metadata for one kinase reference bundle."""

    source: str
    reference: str
    version: str | None = None

    def __post_init__(self) -> None:
        if not self.source.strip():
            msg = "ReferenceBundle source_metadata.source must not be empty"
            raise InputCompatibilityError(msg)
        if not self.reference.strip():
            msg = "ReferenceBundle source_metadata.reference must not be empty"
            raise InputCompatibilityError(msg)
        if self.version is not None and not self.version.strip():
            msg = "ReferenceBundle source_metadata.version must not be empty when provided"
            raise InputCompatibilityError(msg)


@dataclass(frozen=True, slots=True)
class ReferenceBundleProvenance:
    """Provenance describing how a kinase reference bundle was resolved."""

    provider: str
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.provider.strip():
            msg = "ReferenceBundle provenance.provider must not be empty"
            raise InputCompatibilityError(msg)
        normalized_notes = tuple(str(note) for note in self.notes)
        if any(not note.strip() for note in normalized_notes):
            msg = "ReferenceBundle provenance.notes must not contain empty entries"
            raise InputCompatibilityError(msg)
        object.__setattr__(self, "notes", normalized_notes)


@dataclass(frozen=True, slots=True, init=False)
class ReferenceBundle:
    """Typed kinase-prior contract between reference resolution and workflow setup."""

    substrate_map: dict[str, tuple[str, ...]]
    motif_sequences: dict[str, tuple[str, ...]]
    species: str
    source_metadata: ReferenceBundleSourceMetadata
    provenance: ReferenceBundleProvenance

    def __init__(
        self,
        *,
        substrate_map: KinaseSubstrateMap,
        motif_sequences: KinaseMotifSequenceMap,
        species: str,
        source_metadata: ReferenceBundleSourceMetadata,
        provenance: ReferenceBundleProvenance,
    ) -> None:
        normalized_substrate_map = normalize_sequence_mapping(
            substrate_map,
            field_name="substrate_map",
            empty_message="ReferenceBundle substrate_map must not be empty",
        )
        normalized_motif_sequences = normalize_sequence_mapping(
            motif_sequences,
            field_name="motif_sequences",
            empty_message="ReferenceBundle motif_sequences must not be empty",
        )
        resolved_species = str(species).strip()
        if not resolved_species:
            msg = "ReferenceBundle species must not be empty"
            raise InputCompatibilityError(msg)

        _validate_reference_mapping_values(
            normalized_substrate_map,
            field_name="substrate_map",
        )
        _validate_reference_mapping_values(
            normalized_motif_sequences,
            field_name="motif_sequences",
        )

        substrate_kinases = set(normalized_substrate_map)
        motif_kinases = set(normalized_motif_sequences)
        if substrate_kinases != motif_kinases:
            missing_in_motifs = sorted(substrate_kinases - motif_kinases)
            missing_in_substrates = sorted(motif_kinases - substrate_kinases)
            parts: list[str] = []
            if missing_in_motifs:
                parts.append(
                    "missing from motif_sequences: " + ", ".join(missing_in_motifs)
                )
            if missing_in_substrates:
                parts.append(
                    "missing from substrate_map: " + ", ".join(missing_in_substrates)
                )
            msg = "ReferenceBundle kinase sets must match exactly"
            if parts:
                msg = f"{msg} ({'; '.join(parts)})"
            raise InputCompatibilityError(msg)

        build_validated_motif_library(
            motif_sequences=normalized_motif_sequences,
            context="ReferenceBundle motif_sequences",
        )

        object.__setattr__(self, "substrate_map", dict(normalized_substrate_map))
        object.__setattr__(self, "motif_sequences", dict(normalized_motif_sequences))
        object.__setattr__(self, "species", resolved_species)
        object.__setattr__(self, "source_metadata", source_metadata)
        object.__setattr__(self, "provenance", provenance)


@runtime_checkable
class ReferenceProvider(Protocol):
    """Protocol for resolving kinase prior inputs into a ReferenceBundle."""

    def resolve(
        self,
        *,
        species: str,
        reference: str = "auto",
    ) -> ReferenceBundle: ...


@dataclass(frozen=True, slots=True)
class BundledReferenceProvider:
    """Resolve packaged kinase priors for the currently supported species lane."""

    source: str = BUNDLED_REFERENCE_SOURCE
    version: str = BUNDLED_REFERENCE_VERSION

    def resolve(
        self,
        *,
        species: str,
        reference: str = BUNDLED_REFERENCE_AUTO,
    ) -> ReferenceBundle:
        resolved_species = _normalize_bundled_species(species)
        resolved_reference = _normalize_bundled_reference(
            species=resolved_species,
            reference=reference,
        )
        substrate_map = self._load_substrate_map(
            species=resolved_species,
            reference=resolved_reference,
        )
        site_sequences = self._load_site_sequences(
            species=resolved_species,
            reference=resolved_reference,
        )
        motif_sequences = _build_reference_motif_sequences(
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            species=resolved_species,
            reference=resolved_reference,
        )
        return ReferenceBundle(
            substrate_map=substrate_map,
            motif_sequences=motif_sequences,
            species=resolved_species,
            source_metadata=ReferenceBundleSourceMetadata(
                source=self.source,
                reference=resolved_reference,
                version=self.version,
            ),
            provenance=ReferenceBundleProvenance(
                provider=BUNDLED_REFERENCE_PROVIDER_NAME,
                notes=(
                    f"resolved species={resolved_species}",
                    f"resolved reference={resolved_reference}",
                ),
            ),
        )

    @classmethod
    def supported_species(cls) -> tuple[str, ...]:
        return tuple(BUNDLED_REFERENCE_DEFAULTS)

    @classmethod
    def supported_references_for_species(cls, species: str) -> tuple[str, ...]:
        resolved_species = _normalize_bundled_species(species)
        canonical_references = {
            resolved_reference
            for resolved_reference in BUNDLED_REFERENCE_ALIASES[
                resolved_species
            ].values()
            if resolved_reference != BUNDLED_REFERENCE_AUTO
        }
        return tuple(sorted(canonical_references))

    def _load_substrate_map(
        self,
        *,
        species: str,
        reference: str,
    ) -> dict[str, tuple[str, ...]]:
        return _load_grouped_mapping_file(
            _bundled_reference_resource_path(
                species=species,
                reference=reference,
                filename="substrate_map.csv",
            ),
            group_column="kinase",
            value_column="site_id",
        )

    def _load_site_sequences(
        self,
        *,
        species: str,
        reference: str,
    ) -> dict[str, str]:
        return _load_string_mapping_file(
            _bundled_reference_resource_path(
                species=species,
                reference=reference,
                filename="site_sequences.csv",
            ),
            key_column="site_id",
            value_column="centralized_sequence",
        )


@dataclass(slots=True)
class MotifScoringResult:
    """Motif scoring tables and sequence windows for one scoring run."""

    motif_scores: pd.DataFrame
    motif_sizes: pd.Series
    sequence_windows: pd.Series


@dataclass(slots=True)
class ValidatedMotifLibrary:
    """Trusted validated motif library bundle used by motif scoring setup.

    The contained frequency matrices and size series remain mutable pandas data,
    so this wrapper is not presented as an immutable value object.
    """

    motif_frequency_matrices: dict[str, pd.DataFrame]
    motif_sizes: pd.Series


class KinaseMotifScorer:
    """Score phosphosite sequences against per-kinase motif frequency matrices."""

    def __init__(
        self,
        motif_frequency_matrices: Mapping[str, pd.DataFrame],
        motif_sizes: pd.Series,
        flank_size: int = 7,
    ) -> None:
        validate_non_negative_int(flank_size, name="flank_size")
        if not motif_frequency_matrices:
            msg = "motif_frequency_matrices must not be empty"
            raise InputCompatibilityError(msg)

        self.motif_frequency_matrices = {
            kinase: _coerce_frequency_matrix(matrix)
            for kinase, matrix in motif_frequency_matrices.items()
        }
        motif_widths = {
            _motif_matrix_width(matrix)
            for matrix in self.motif_frequency_matrices.values()
        }
        if len(motif_widths) != 1:
            msg = "All motif frequency matrices must have the same window width"
            raise InputCompatibilityError(msg)
        self._motif_width = motif_widths.pop()
        self._motif_frequency_values = {
            kinase: matrix.to_numpy(dtype=float, copy=False)
            for kinase, matrix in self.motif_frequency_matrices.items()
        }
        self.motif_sizes = motif_sizes.astype(float).copy()
        self.flank_size = flank_size

        missing = [
            kinase
            for kinase in self.motif_frequency_matrices
            if kinase not in self.motif_sizes.index
        ]
        if missing:
            msg = f"motif_sizes is missing entries for: {', '.join(missing)}"
            raise InputCompatibilityError(msg)

    @classmethod
    def from_substrate_sequences(
        cls,
        motif_sequences: Mapping[str, Sequence[str]],
        flank_size: int = 7,
    ) -> KinaseMotifScorer:
        validated_library = build_validated_motif_library(
            motif_sequences=motif_sequences,
            flank_size=flank_size,
        )
        return cls(
            motif_frequency_matrices=validated_library.motif_frequency_matrices,
            motif_sizes=validated_library.motif_sizes,
            flank_size=flank_size,
        )

    def score_sequences(
        self,
        seqs: Mapping[str, str] | Sequence[str] | pd.Series,
        site_index: Sequence[str] | None = None,
        min_motif_size: int = 1,
    ) -> MotifScoringResult:
        validate_positive_int(min_motif_size, name="min_motif_size")
        windows = _coerce_sequence_series(
            seqs=seqs,
            site_index=site_index,
            flank_size=self.flank_size,
        )

        kinases = [
            kinase
            for kinase in self.motif_frequency_matrices
            if float(self.motif_sizes.loc[kinase]) >= float(min_motif_size)
        ]
        motif_scores = pd.DataFrame(
            np.nan,
            index=windows.index.copy(),
            columns=kinases,
            dtype=float,
        )
        if kinases:
            encoded_windows = _encode_sequence_positions(
                windows,
                width=self._motif_width,
            )
            for kinase in kinases:
                motif_scores.loc[:, kinase] = _score_encoded_sequences(
                    encoded_sequences=encoded_windows,
                    frequency_values=self._motif_frequency_values[kinase],
                )
        motif_scores = minmax_scale_columns(motif_scores)

        motif_sizes = self.motif_sizes.loc[kinases].astype(float).copy()
        motif_sizes.index.name = "kinase"
        return MotifScoringResult(
            motif_scores=motif_scores,
            motif_sizes=motif_sizes,
            sequence_windows=windows,
        )


def create_frequency_matrix(
    substrates_seq: Sequence[str] | pd.Series,
    flank_size: int = 7,
) -> pd.DataFrame:
    """Create an amino-acid frequency matrix from substrate sequences."""

    validate_non_negative_int(flank_size, name="flank_size")
    windows = _coerce_sequence_series(substrates_seq, flank_size=flank_size)
    return _build_frequency_matrix_from_windows(
        windows,
        context="substrates_seq",
    )


def build_validated_motif_library(
    motif_sequences: Mapping[str, Sequence[str]],
    *,
    flank_size: int = 7,
    context: str = "motif_sequences",
) -> ValidatedMotifLibrary:
    validate_non_negative_int(flank_size, name="flank_size")

    matrices: dict[str, pd.DataFrame] = {}
    sizes: dict[str, float] = {}
    widths: set[int] = set()

    for kinase, sequences in motif_sequences.items():
        windows = _coerce_sequence_series(sequences, flank_size=flank_size)
        try:
            matrices[kinase] = _build_frequency_matrix_from_windows(
                windows,
                context=f"{context} for kinase {kinase}",
            )
        except TableSchemaError as error:
            if "same window length" in str(error):
                msg = f"{context} for kinase {kinase} must use a consistent sequence width"
                raise InputCompatibilityError(msg) from error
            raise
        sizes[kinase] = float(len(windows))
        widths.add(matrices[kinase].shape[1])

    if len(widths) > 1:
        msg = f"{context} must use the same sequence width across kinases"
        raise InputCompatibilityError(msg)

    return ValidatedMotifLibrary(
        motif_frequency_matrices=matrices,
        motif_sizes=pd.Series(sizes, dtype=float),
    )


def frequency_scoring(
    sequence_list: Sequence[str] | pd.Series,
    frequency_mat: pd.DataFrame,
) -> pd.Series:
    """Score phosphosite windows against a single motif frequency matrix."""

    frequency_mat = _coerce_frequency_matrix(frequency_mat)
    sequences = _coerce_sequence_series(sequence_list, flank_size=None)
    frequency_values = frequency_mat.to_numpy(dtype=float, copy=False)
    encoded_sequences = _encode_sequence_positions(
        sequences,
        width=frequency_values.shape[1],
    )
    score_values = _score_encoded_sequences(
        encoded_sequences=encoded_sequences,
        frequency_values=frequency_values,
    )
    return pd.Series(score_values, index=sequences.index.copy(), dtype=float)


def score_phosphosite_motifs(
    seqs: Mapping[str, str] | Sequence[str] | pd.Series,
    motif_frequency_matrices: Mapping[str, pd.DataFrame],
    motif_sizes: pd.Series,
    site_index: Sequence[str] | None = None,
    min_motif_size: int = 1,
    flank_size: int = 7,
) -> MotifScoringResult:
    scorer = KinaseMotifScorer(
        motif_frequency_matrices=motif_frequency_matrices,
        motif_sizes=motif_sizes,
        flank_size=flank_size,
    )
    return scorer.score_sequences(
        seqs=seqs,
        site_index=site_index,
        min_motif_size=min_motif_size,
    )


def minmax_scale_columns(mat: pd.DataFrame) -> pd.DataFrame:
    """Apply PhosR-style column-wise min-max scaling."""

    scaled = mat.astype(float).copy()
    for column in scaled.columns:
        values = scaled.loc[:, column]
        min_value = float(values.min())
        max_value = float(values.max())
        denominator = max_value - min_value
        if denominator == 0.0:
            scaled.loc[:, column] = np.nan
        else:
            scaled.loc[:, column] = (values - min_value) / denominator
    return scaled


def _build_frequency_matrix_from_windows(
    windows: pd.Series,
    *,
    context: str,
) -> pd.DataFrame:
    if windows.empty:
        msg = f"{context} must contain at least one sequence"
        raise TableSchemaError(msg)

    width = len(str(windows.iloc[0]))
    if any(len(str(window)) != width for window in windows):
        msg = "All sequences must have the same window length"
        raise TableSchemaError(msg)

    frequency_values = np.zeros((len(AMINO_ACIDS), width), dtype=float)
    encoded_windows = _encode_sequence_positions(windows, width=width)
    valid_rows, valid_cols = np.nonzero(encoded_windows >= 0)
    if valid_rows.size > 0:
        np.add.at(
            frequency_values,
            (encoded_windows[valid_rows, valid_cols], valid_cols),
            1.0,
        )

    frequency_values /= float(len(windows))
    return pd.DataFrame(
        frequency_values,
        index=list(AMINO_ACIDS),
        columns=[f"p{i}" for i in range(1, width + 1)],
        dtype=float,
    )


def _encode_sequence_positions(
    sequences: Sequence[object] | pd.Series,
    width: int,
) -> np.ndarray:
    """Encode sequence characters into amino-acid row indices for vectorised scoring."""
    encoded = np.full(
        (len(sequences), width),
        _INVALID_AMINO_ACID_INDEX,
        dtype=np.int16,
    )
    if width == 0:
        return encoded

    for row_index, sequence in enumerate(sequences):
        if sequence is None or pd.isna(sequence):
            continue

        text = str(sequence)[:width]
        if text == "":
            continue

        code_points = np.fromiter(
            (ord(character) for character in text),
            dtype=np.int32,
            count=len(text),
        )
        valid_ascii = code_points < _ASCII_LOOKUP_SIZE
        if not np.any(valid_ascii):
            continue

        valid_positions = np.flatnonzero(valid_ascii)
        encoded[row_index, valid_positions] = _AMINO_ACID_INDEX_LOOKUP[
            code_points[valid_positions]
        ]
    return encoded


def _score_encoded_sequences(
    encoded_sequences: np.ndarray,
    frequency_values: np.ndarray,
) -> np.ndarray:
    """Score encoded sequence windows against a frequency matrix using NumPy indexing."""
    if encoded_sequences.size == 0:
        return np.zeros(encoded_sequences.shape[0], dtype=float)

    valid_mask = encoded_sequences >= 0
    safe_indices = np.where(valid_mask, encoded_sequences, 0)
    position_indices = np.arange(encoded_sequences.shape[1])
    position_scores = frequency_values[safe_indices, position_indices]
    position_scores = np.where(valid_mask, position_scores, 0.0)
    return position_scores.sum(axis=1, dtype=float)


def _load_grouped_mapping_file(
    path: Path,
    *,
    group_column: str,
    value_column: str,
) -> dict[str, tuple[str, ...]]:
    frame = pd.read_csv(path)
    grouped: dict[str, list[str]] = {}
    for group, value in frame.loc[:, [group_column, value_column]].itertuples(
        index=False
    ):
        grouped.setdefault(str(group).strip(), []).append(str(value).strip())
    return {key: tuple(values) for key, values in grouped.items()}


def _load_string_mapping_file(
    path: Path,
    *,
    key_column: str,
    value_column: str,
) -> dict[str, str]:
    frame = pd.read_csv(path)
    return {
        str(key).strip(): str(value).strip()
        for key, value in frame.loc[:, [key_column, value_column]].itertuples(
            index=False
        )
    }


def _build_reference_motif_sequences(
    *,
    substrate_map: Mapping[str, Sequence[str]],
    site_sequences: Mapping[str, str],
    species: str,
    reference: str,
) -> dict[str, tuple[str, ...]]:
    motif_sequences: dict[str, tuple[str, ...]] = {}
    missing_sites: set[str] = set()
    for kinase, site_ids in substrate_map.items():
        sequences: list[str] = []
        for site_id in site_ids:
            sequence = site_sequences.get(str(site_id))
            if sequence is None:
                missing_sites.add(str(site_id))
                continue
            sequences.append(str(sequence))
        motif_sequences[str(kinase)] = tuple(sequences)
    if missing_sites:
        missing_text = ", ".join(sorted(missing_sites))
        msg = (
            "BundledReferenceProvider reference data is incomplete for "
            f"species '{species}' and reference '{reference}'; missing site sequences for: {missing_text}"
        )
        raise InputCompatibilityError(msg)
    return motif_sequences


def _normalize_bundled_species(species: str) -> str:
    normalized = str(species).strip().lower()
    resolved_species = BUNDLED_REFERENCE_SPECIES_ALIASES.get(normalized)
    if resolved_species is None:
        supported = ", ".join(sorted(BUNDLED_REFERENCE_DEFAULTS))
        msg = (
            f"Unsupported bundled reference species '{species}'. "
            f"Supported species: {supported}"
        )
        raise InputCompatibilityError(msg)
    return resolved_species


def _normalize_bundled_reference(*, species: str, reference: str) -> str:
    resolved_reference = BUNDLED_REFERENCE_ALIASES[species].get(
        str(reference).strip().lower()
    )
    if resolved_reference is None:
        supported = ", ".join(
            BundledReferenceProvider.supported_references_for_species(species)
        )
        msg = (
            f"Unsupported bundled reference '{reference}' for species '{species}'. "
            f"Supported references: {supported}"
        )
        raise InputCompatibilityError(msg)
    return resolved_reference


def _bundled_reference_resource_path(
    *,
    species: str,
    reference: str,
    filename: str,
) -> Path:
    resource = resources.files("phospy").joinpath(
        "data",
        "reference_bundles",
        species,
        reference,
        filename,
    )
    if not resource.is_file():
        msg = (
            "BundledReferenceProvider could not find packaged reference data for "
            f"species '{species}' and reference '{reference}' ({filename})"
        )
        raise InputCompatibilityError(msg)
    with resources.as_file(resource) as resolved_path:
        return resolved_path


def _validate_reference_mapping_values(
    mapping: dict[str, tuple[str, ...]],
    *,
    field_name: str,
) -> None:
    empty_kinases = sorted(kinase for kinase, values in mapping.items() if not values)
    if empty_kinases:
        msg = f"ReferenceBundle {field_name} entries must not be empty: {', '.join(empty_kinases)}"
        raise InputCompatibilityError(msg)


def _motif_matrix_width(frequency_mat: pd.DataFrame) -> int:
    return int(frequency_mat.shape[1])


def _coerce_frequency_matrix(frequency_mat: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(frequency_mat, pd.DataFrame):
        msg = "frequency_mat must be a pandas DataFrame"
        raise TypeError(msg)

    matrix = frequency_mat.astype(float).copy()
    matrix = matrix.reindex(index=list(AMINO_ACIDS), fill_value=0.0)
    matrix.columns = [f"p{i}" for i in range(1, matrix.shape[1] + 1)]
    return matrix


def _coerce_sequence_series(
    seqs: Mapping[str, str] | Sequence[str] | pd.Series,
    site_index: Sequence[str] | None = None,
    flank_size: int | None = 7,
) -> pd.Series:
    if isinstance(seqs, pd.Series):
        series = seqs
    elif isinstance(seqs, Mapping):
        series = pd.Series(dict(seqs), dtype=object)
    else:
        seq_list = list(seqs)
        if site_index is None:
            series = pd.Series(seq_list, dtype=object)
        else:
            if len(seq_list) != len(site_index):
                msg = "site_index must have the same length as seqs"
                raise TableSchemaError(msg)
            series = pd.Series(seq_list, index=list(site_index), dtype=object)

    if site_index is not None:
        missing = [site for site in site_index if site not in series.index]
        if missing:
            msg = f"seqs is missing entries for: {', '.join(missing)}"
            raise TableSchemaError(msg)
        series = series.loc[list(site_index)]

    return series.map(lambda value: _extract_sequence_window(value, flank_size))


def _extract_sequence_window(value: object, flank_size: int | None) -> str:
    if value is None or pd.isna(value):
        return np.nan  # type: ignore[return-value]

    sequence = str(value).upper()
    if sequence == "":
        return "_"

    if flank_size is None:
        return _validate_sequence_window(sequence)

    window_size = (2 * flank_size) + 1
    if len(sequence) <= window_size:
        return _validate_sequence_window(sequence)

    mid = len(sequence) // 2
    start = mid - flank_size
    stop = mid + flank_size + 1
    return _validate_sequence_window(sequence[start:stop])


def _validate_sequence_window(sequence: str) -> str:
    invalid_characters = sorted(
        {
            character
            for character in sequence
            if character not in _ALLOWED_SEQUENCE_CHARACTERS
        }
    )
    if invalid_characters:
        invalid_text = ", ".join(repr(character) for character in invalid_characters)
        msg = (
            "sequence contains invalid amino-acid characters: "
            f"{invalid_text}; allowed characters are the 20 standard amino acids and '_'"
        )
        raise TableSchemaError(msg)
    return sequence

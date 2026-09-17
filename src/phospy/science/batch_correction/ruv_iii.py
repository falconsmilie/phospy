"""Replicate-aware RUV-III numerical correction.

The kernel in this module is intentionally separate from the established
``sps_ruv_style`` executor.  It follows the finite-``k`` branch of
``ruv::RUVIII`` with samples as observations and phosphosites as features.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from numbers import Integral
from typing import cast

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.sites.validation import require_site_key_index

RUV_III_ALGORITHM_ID = "replicate_aware_ruv_iii_v1"
RUV_III_METHOD = "ruv_iii"


@dataclass(frozen=True, slots=True)
class RuvIIIReplicateStructure:
    """Validated replicate-set assignments in an explicit sample order.

    Every sample belongs to exactly one replicate set and every set contains at
    least two samples.  The explicit order prevents replicate relationships
    from being inferred from matrix-column adjacency.
    """

    sample_order: tuple[str, ...]
    set_by_sample: tuple[str, ...]

    def __post_init__(self) -> None:
        sample_order = _require_label_tuple(
            self.sample_order,
            field_name="RUV-III replicate sample_order",
        )
        set_by_sample = _require_label_tuple(
            self.set_by_sample,
            field_name="RUV-III replicate set_by_sample",
            require_unique=False,
        )
        if len(sample_order) != len(set_by_sample):
            raise PhosPyInputError(
                "RUV-III replicate structure requires one replicate-set "
                "assignment for every sample"
            )
        if len(sample_order) < 2:
            raise PhosPyInputError(
                "RUV-III replicate structure requires at least two samples"
            )

        counts: dict[str, int] = {}
        for set_id in set_by_sample:
            counts[set_id] = counts.get(set_id, 0) + 1
        unsuitable = tuple(set_id for set_id, count in counts.items() if count < 2)
        if unsuitable:
            raise PhosPyInputError(
                "RUV-III replicate sets must each contain at least two samples; "
                "singleton set(s): " + ", ".join(repr(value) for value in unsuitable)
            )

        object.__setattr__(self, "sample_order", sample_order)
        object.__setattr__(self, "set_by_sample", set_by_sample)

    @classmethod
    def from_assignments(
        cls,
        *,
        sample_order: Sequence[str],
        replicate_by_sample: Mapping[str, str],
    ) -> RuvIIIReplicateStructure:
        """Build validated replicate sets from explicit name-based assignments."""

        ordered_samples = tuple(sample_order)
        _require_label_tuple(
            ordered_samples,
            field_name="RUV-III replicate sample_order",
        )
        assignment_keys = tuple(replicate_by_sample.keys())
        missing = tuple(
            sample for sample in ordered_samples if sample not in replicate_by_sample
        )
        unexpected = tuple(
            sample for sample in assignment_keys if sample not in ordered_samples
        )
        if missing or unexpected:
            details: list[str] = []
            if missing:
                details.append(
                    "missing assignments=" + ", ".join(repr(value) for value in missing)
                )
            if unexpected:
                details.append(
                    "unexpected samples="
                    + ", ".join(repr(value) for value in unexpected)
                )
            raise PhosPyInputError(
                "RUV-III replicate assignments must exactly cover sample_order; "
                + "; ".join(details)
            )
        return cls(
            sample_order=ordered_samples,
            set_by_sample=tuple(
                replicate_by_sample[sample] for sample in ordered_samples
            ),
        )

    @property
    def set_ids(self) -> tuple[str, ...]:
        """Return replicate-set identifiers in deterministic first-seen order."""

        return tuple(dict.fromkeys(self.set_by_sample))

    @property
    def set_count(self) -> int:
        """Return the number of replicate sets."""

        return len(self.set_ids)

    def to_mapping_matrix(self) -> np.ndarray:
        """Return the samples-by-replicate-sets indicator matrix."""

        set_position = {
            set_id: position for position, set_id in enumerate(self.set_ids)
        }
        mapping = np.zeros((len(self.sample_order), self.set_count), dtype="float64")
        for sample_position, set_id in enumerate(self.set_by_sample):
            mapping[sample_position, set_position[set_id]] = 1.0
        return mapping


@dataclass(frozen=True, slots=True)
class RuvIIIDiagnostics:
    """Immutable, comparable diagnostics for one RUV-III fit."""

    method: str
    algorithm_id: str
    requested_k: int
    effective_k: int
    sample_count: int
    feature_count: int
    control_count: int
    replicate_set_count: int
    input_matrix_rank: int
    replicate_mapping_rank: int
    replicate_residual_degrees_of_freedom: int
    replicate_residual_rank: int
    control_loading_rank: int
    estimated_unwanted_factor_dimensions: tuple[int, int]
    residual_singular_values: tuple[float, ...]

    def to_payload(self) -> dict[str, object]:
        """Return deterministic JSON-compatible diagnostics."""

        return {
            "method": self.method,
            "algorithm_id": self.algorithm_id,
            "requested_k": self.requested_k,
            "effective_k": self.effective_k,
            "sample_count": self.sample_count,
            "feature_count": self.feature_count,
            "control_count": self.control_count,
            "replicate_set_count": self.replicate_set_count,
            "input_matrix_rank": self.input_matrix_rank,
            "replicate_mapping_rank": self.replicate_mapping_rank,
            "replicate_residual_degrees_of_freedom": (
                self.replicate_residual_degrees_of_freedom
            ),
            "replicate_residual_rank": self.replicate_residual_rank,
            "control_loading_rank": self.control_loading_rank,
            "estimated_unwanted_factor_dimensions": list(
                self.estimated_unwanted_factor_dimensions
            ),
            "residual_singular_values": list(self.residual_singular_values),
        }


@dataclass(frozen=True, slots=True, eq=False)
class RuvIIIResult:
    """Corrected phosphosite matrix and compact factor-estimation output."""

    corrected_matrix: pd.DataFrame
    estimated_unwanted_factors: pd.DataFrame
    diagnostics: RuvIIIDiagnostics

    @property
    def corrected(self) -> pd.DataFrame:
        """Return the corrected feature-by-sample matrix."""

        return self.corrected_matrix


class RuvIIIKernel:
    """Apply replicate-aware RUV-III correction to a finite matrix."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        control_site_keys: Sequence[str],
        replicate_structure: RuvIIIReplicateStructure,
        k: int,
    ) -> RuvIIIResult:
        """Fit and apply the finite-``k`` RUV-III model."""

        values, site_keys, sample_ids = _validated_matrix(phospho)
        requested_k = _require_k(k)
        if tuple(replicate_structure.sample_order) != sample_ids:
            raise PhosPyInputError(
                "RUV-III replicate sample_order must exactly match phospho columns"
            )
        control_positions = _control_positions(
            site_keys=site_keys,
            control_site_keys=control_site_keys,
        )
        mapping = replicate_structure.to_mapping_matrix()
        mapping_rank = _matrix_rank(mapping)
        if mapping_rank != replicate_structure.set_count:
            raise PhosPyInputError("RUV-III replicate mapping matrix is rank-deficient")
        residual_degrees_of_freedom = int(values.shape[0]) - mapping_rank
        if residual_degrees_of_freedom < 1:
            raise PhosPyInputError(
                "RUV-III replicate structure leaves no within-set sample degrees "
                "of freedom"
            )
        if requested_k > residual_degrees_of_freedom:
            raise PhosPyInputError(
                f"RUV-III k={requested_k} exceeds the available replicate-residual "
                f"degrees of freedom={residual_degrees_of_freedom}"
            )
        if requested_k > len(control_positions):
            raise PhosPyInputError(
                f"RUV-III k={requested_k} requires at least {requested_k} negative "
                f"controls; received control_count={len(control_positions)}"
            )

        try:
            with np.errstate(over="ignore", invalid="ignore"):
                replicate_means = np.linalg.solve(
                    mapping.T @ mapping,
                    mapping.T @ values,
                )
                residual = values - mapping @ replicate_means
            _require_finite_intermediate(
                residual,
                name="replicate-residual matrix",
            )
            left_vectors, residual_singular_values, _ = np.linalg.svd(
                residual,
                full_matrices=False,
            )
            _require_finite_intermediate(
                residual_singular_values,
                name="replicate-residual singular values",
            )
        except np.linalg.LinAlgError as error:
            raise PhosPyInputError(
                "RUV-III numerical decomposition failed while estimating "
                "replicate-residual variation; verify matrix scaling and rank"
            ) from error
        residual_rank = _rank_from_singular_values(
            matrix_shape=residual.shape,
            singular_values=residual_singular_values,
        )
        if requested_k > residual_rank:
            raise PhosPyInputError(
                f"RUV-III k={requested_k} is not estimable from replicate-residual "
                f"matrix rank={residual_rank}"
            )
        if _singular_value_tie_at_cutoff(
            singular_values=residual_singular_values,
            matrix_shape=residual.shape,
            k=requested_k,
        ):
            raise PhosPyInputError(
                f"RUV-III k={requested_k} splits a numerically tied singular-value "
                "block at the requested factor cutoff; choose k that includes or "
                "excludes the complete tied block"
            )

        if requested_k == 0:
            corrected_values = values.copy()
            factors = np.empty((int(values.shape[0]), 0), dtype="float64")
            control_loading_rank = 0
        else:
            selected_vectors = _canonicalize_vector_basis(
                left_vectors[:, :requested_k],
                singular_values=residual_singular_values[:requested_k],
                sample_ids=sample_ids,
                matrix_shape=residual.shape,
            )
            with np.errstate(over="ignore", invalid="ignore"):
                loadings = selected_vectors.T @ values
            _require_finite_intermediate(loadings, name="unwanted-factor loadings")
            control_loadings = loadings[:, list(control_positions)]
            try:
                factor_solution, _, recovered_rank, _ = np.linalg.lstsq(
                    control_loadings.T,
                    values[:, list(control_positions)].T,
                    rcond=None,
                )
            except np.linalg.LinAlgError as error:
                raise PhosPyInputError(
                    "RUV-III least-squares factor recovery failed; verify negative "
                    "controls, matrix scaling, and requested k"
                ) from error
            control_loading_rank = int(recovered_rank)
            if control_loading_rank < requested_k:
                raise PhosPyInputError(
                    "RUV-III negative-control loading matrix is rank-deficient for "
                    f"k={requested_k}; control_loading_rank={control_loading_rank}"
                )
            factors = factor_solution.T
            _require_finite_intermediate(factors, name="estimated unwanted factors")
            with np.errstate(over="ignore", invalid="ignore"):
                corrected_values = values - factors @ loadings

        if not np.isfinite(corrected_values).all() or not np.isfinite(factors).all():
            raise PhosPyInputError(
                "RUV-III numerical estimation produced non-finite output"
            )

        corrected_matrix = pd.DataFrame(
            corrected_values.T,
            index=phospho.index.copy(),
            columns=phospho.columns.copy(),
        )
        estimated_factors = pd.DataFrame(
            factors,
            index=phospho.columns.copy(),
            columns=pd.Index(
                tuple(
                    f"unwanted_factor_{position + 1}" for position in range(requested_k)
                ),
                name="factor",
            ),
        )
        diagnostics = RuvIIIDiagnostics(
            method=RUV_III_METHOD,
            algorithm_id=RUV_III_ALGORITHM_ID,
            requested_k=requested_k,
            effective_k=requested_k,
            sample_count=int(values.shape[0]),
            feature_count=int(values.shape[1]),
            control_count=len(control_positions),
            replicate_set_count=replicate_structure.set_count,
            input_matrix_rank=_matrix_rank(values),
            replicate_mapping_rank=mapping_rank,
            replicate_residual_degrees_of_freedom=residual_degrees_of_freedom,
            replicate_residual_rank=residual_rank,
            control_loading_rank=control_loading_rank,
            estimated_unwanted_factor_dimensions=(int(factors.shape[0]), requested_k),
            residual_singular_values=tuple(
                float(value) for value in residual_singular_values[:requested_k]
            ),
        )
        return RuvIIIResult(
            corrected_matrix=corrected_matrix,
            estimated_unwanted_factors=estimated_factors,
            diagnostics=diagnostics,
        )


def run_ruv_iii(
    phospho: pd.DataFrame,
    *,
    control_site_keys: Sequence[str],
    replicate_structure: RuvIIIReplicateStructure,
    k: int,
) -> RuvIIIResult:
    """Run the standalone replicate-aware RUV-III scientific kernel."""

    return RuvIIIKernel().run(
        phospho=phospho,
        control_site_keys=control_site_keys,
        replicate_structure=replicate_structure,
        k=k,
    )


def _validated_matrix(
    phospho: pd.DataFrame,
) -> tuple[np.ndarray, tuple[str, ...], tuple[str, ...]]:
    if not isinstance(cast(object, phospho), pd.DataFrame):
        raise PhosPyInputError("RUV-III requires phospho to be a pandas DataFrame")
    if phospho.shape[0] < 1:
        raise PhosPyInputError("RUV-III requires at least one phosphosite feature")
    if phospho.shape[1] < 2:
        raise PhosPyInputError("RUV-III requires at least two samples")
    if phospho.index.name != "site_key":
        raise PhosPyInputError("RUV-III requires phospho.index.name='site_key'")
    require_site_key_index(
        phospho.index,
        field_name="RUV-III phospho.index",
        error_type=PhosPyInputError,
    )
    site_keys = tuple(str(value) for value in phospho.index.tolist())
    sample_ids = _require_label_tuple(
        tuple(phospho.columns.tolist()),
        field_name="RUV-III phospho sample columns",
    )
    complex_columns = tuple(
        str(column)
        for column in phospho.columns
        if pd.api.types.is_complex_dtype(phospho.loc[:, column])
    )
    if complex_columns:
        raise PhosPyInputError(
            "RUV-III requires real-valued phospho intensities; complex-valued "
            "column(s): " + ", ".join(repr(value) for value in complex_columns)
        )
    non_numeric = tuple(
        str(column)
        for column in phospho.columns
        if (
            not pd.api.types.is_numeric_dtype(phospho.loc[:, column])
            or pd.api.types.is_bool_dtype(phospho.loc[:, column])
        )
    )
    if non_numeric:
        raise PhosPyInputError(
            "RUV-III requires numeric, non-boolean phospho columns; invalid "
            "column(s): " + ", ".join(repr(value) for value in non_numeric)
        )
    values = phospho.to_numpy(dtype="float64", copy=True)
    if not np.isfinite(values).all():
        row_position, column_position = np.argwhere(~np.isfinite(values))[0]
        raise PhosPyInputError(
            "RUV-III requires a complete finite phosphosite matrix; found a "
            f"non-finite value at site_key={site_keys[int(row_position)]!r}, "
            f"sample_id={sample_ids[int(column_position)]!r}. The low-level "
            "kernel does not impute missing values."
        )
    return values.T, site_keys, sample_ids


def _control_positions(
    *,
    site_keys: tuple[str, ...],
    control_site_keys: Sequence[str],
) -> tuple[int, ...]:
    if isinstance(cast(object, control_site_keys), (str, bytes)):
        raise PhosPyInputError(
            "RUV-III control_site_keys must be a sequence of governed site_key values"
        )
    controls = _require_label_tuple(
        tuple(control_site_keys),
        field_name="RUV-III control_site_keys",
    )
    if not controls:
        raise PhosPyInputError("RUV-III requires at least one negative-control site")
    position_by_key = {
        site_key: position for position, site_key in enumerate(site_keys)
    }
    missing = tuple(
        site_key for site_key in controls if site_key not in position_by_key
    )
    if missing:
        raise PhosPyInputError(
            "RUV-III control_site_keys are absent from the phospho site_key index: "
            + ", ".join(repr(value) for value in missing)
        )
    return tuple(position_by_key[site_key] for site_key in controls)


def _require_k(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise PhosPyInputError(
            "RUV-III k must be an integer greater than or equal to 0"
        )
    resolved = int(value)
    if resolved < 0:
        raise PhosPyInputError("RUV-III k must be greater than or equal to 0")
    return resolved


def _require_label_tuple(
    values: tuple[object, ...],
    *,
    field_name: str,
    require_unique: bool = True,
) -> tuple[str, ...]:
    labels: list[str] = []
    for position, value in enumerate(values):
        if not isinstance(value, str) or value.strip() == "":
            raise PhosPyInputError(
                f"{field_name}[{position}] must be a non-empty string"
            )
        labels.append(value)
    if require_unique and len(set(labels)) != len(labels):
        duplicates = tuple(
            label for label in dict.fromkeys(labels) if labels.count(label) > 1
        )
        raise PhosPyInputError(
            f"{field_name} must be unique; duplicate label(s): "
            + ", ".join(repr(value) for value in duplicates)
        )
    return tuple(labels)


def _matrix_rank(matrix: np.ndarray) -> int:
    try:
        return int(np.linalg.matrix_rank(matrix))
    except np.linalg.LinAlgError as error:
        raise PhosPyInputError(
            "RUV-III numerical rank estimation failed; verify matrix scaling and "
            "finite input values"
        ) from error


def _rank_from_singular_values(
    *,
    matrix_shape: tuple[int, int],
    singular_values: np.ndarray,
) -> int:
    if singular_values.size == 0:
        return 0
    tolerance = (
        max(matrix_shape) * np.finfo(np.float64).eps * float(np.max(singular_values))
    )
    return int(np.sum(singular_values > tolerance))


def _singular_value_tie_at_cutoff(
    *,
    singular_values: np.ndarray,
    matrix_shape: tuple[int, int],
    k: int,
) -> bool:
    if k == 0 or k >= int(singular_values.size):
        return False
    return _singular_values_are_tied(
        left=float(singular_values[k - 1]),
        right=float(singular_values[k]),
        singular_values=singular_values,
        matrix_shape=matrix_shape,
    )


def _singular_values_are_tied(
    *,
    left: float,
    right: float,
    singular_values: np.ndarray,
    matrix_shape: tuple[int, int],
) -> bool:
    scale = max(abs(left), abs(right))
    rank_tolerance = (
        max(matrix_shape) * np.finfo(np.float64).eps * float(np.max(singular_values))
    )
    tie_tolerance = max(rank_tolerance, 1e-12 * scale)
    return abs(left - right) <= tie_tolerance


def _require_finite_intermediate(matrix: np.ndarray, *, name: str) -> None:
    if not np.isfinite(matrix).all():
        raise PhosPyInputError(
            f"RUV-III numerical estimation produced a non-finite {name}; "
            "rescale the finite input matrix and retry"
        )


def _canonicalize_vector_basis(
    vectors: np.ndarray,
    *,
    singular_values: np.ndarray,
    sample_ids: tuple[str, ...],
    matrix_shape: tuple[int, int],
) -> np.ndarray:
    canonical = vectors.copy()
    block_start = 0
    component_count = int(singular_values.size)
    while block_start < component_count:
        block_end = block_start + 1
        while block_end < component_count and _singular_values_are_tied(
            left=float(singular_values[block_end - 1]),
            right=float(singular_values[block_end]),
            singular_values=singular_values,
            matrix_shape=matrix_shape,
        ):
            block_end += 1
        if block_end - block_start > 1:
            canonical[:, block_start:block_end] = _canonical_basis_for_subspace(
                vectors[:, block_start:block_end],
                sample_ids=sample_ids,
            )
        block_start = block_end
    return _canonicalize_vector_signs(canonical, sample_ids=sample_ids)


def _canonical_basis_for_subspace(
    vectors: np.ndarray,
    *,
    sample_ids: tuple[str, ...],
) -> np.ndarray:
    dimension = int(vectors.shape[1])
    basis = np.zeros_like(vectors)
    basis_count = 0
    tolerance = max(vectors.shape) * np.finfo(np.float64).eps * 10.0
    lexical_positions = sorted(range(len(sample_ids)), key=sample_ids.__getitem__)
    for sample_position in lexical_positions:
        candidate = vectors @ vectors[sample_position, :]
        for basis_position in range(basis_count):
            existing = basis[:, basis_position]
            candidate -= existing * float(existing @ candidate)
        norm = float(np.linalg.norm(candidate))
        if norm <= tolerance:
            continue
        basis[:, basis_count] = candidate / norm
        basis_count += 1
        if basis_count == dimension:
            return basis
    raise PhosPyInputError(
        "RUV-III could not construct a deterministic basis for a tied "
        "singular-value subspace; choose a different k"
    )


def _canonicalize_vector_signs(
    vectors: np.ndarray,
    *,
    sample_ids: tuple[str, ...],
) -> np.ndarray:
    canonical = vectors.copy()
    lexical_positions = tuple(
        sorted(range(len(sample_ids)), key=sample_ids.__getitem__)
    )
    for column_position in range(int(canonical.shape[1])):
        vector = canonical[:, column_position]
        maximum = float(np.max(np.abs(vector)))
        candidates = {
            position
            for position, value in enumerate(vector)
            if np.isclose(abs(float(value)), maximum, rtol=1e-12, atol=0.0)
        }
        anchor = next(
            position for position in lexical_positions if position in candidates
        )
        if float(vector[anchor]) < 0.0:
            canonical[:, column_position] *= -1.0
    return canonical


__all__ = [
    "RUV_III_ALGORITHM_ID",
    "RUV_III_METHOD",
    "RuvIIIDiagnostics",
    "RuvIIIKernel",
    "RuvIIIReplicateStructure",
    "RuvIIIResult",
    "run_ruv_iii",
]

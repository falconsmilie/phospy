"""Kinase Library resource validator."""

from __future__ import annotations

import math
from typing import cast

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.provenance.models import KinaseLibraryResourceProvenance
from phospy.science.references.kinase_library import (
    KinaseLibraryMatrix,
    KinaseLibraryResource,
)
from phospy.science.references.models import SequenceWindowDefinition


class KinaseLibraryResourceValidator:
    """Validate the stable Kinase Library-style resource contract."""

    def run(self, resource: KinaseLibraryResource) -> None:
        if not isinstance(resource, KinaseLibraryResource):
            raise ReferenceValidationError(
                "kinase_library resource must be KinaseLibraryResource"
            )
        self._validate_sequence_window(resource.sequence_window)
        self._validate_metadata(resource)
        self._validate_matrices(
            matrices=resource.matrices,
            sequence_window=resource.sequence_window,
        )
        self._validate_provenance(resource.provenance)

    def _validate_metadata(self, resource: KinaseLibraryResource) -> None:
        for field_name, value in (
            ("source_name", resource.source_name),
            ("source_version", resource.source_version),
            ("score_scale", resource.score_scale),
            ("license", resource.license),
        ):
            if not isinstance(value, str) or not value.strip():
                raise ReferenceValidationError(
                    f"kinase_library.{field_name} must be a non-empty string"
                )
        if not resource.organisms:
            raise ReferenceValidationError("kinase_library.organisms must not be empty")
        for organism in resource.organisms:
            if not isinstance(organism, str) or not organism.strip():
                raise ReferenceValidationError(
                    "kinase_library.organisms must contain non-empty strings"
                )

    @staticmethod
    def _validate_sequence_window(
        sequence_window: SequenceWindowDefinition,
    ) -> None:
        if not isinstance(cast(object, sequence_window), SequenceWindowDefinition):
            raise ReferenceValidationError(
                "kinase_library.sequence_window must be SequenceWindowDefinition"
            )
        if sequence_window.upstream_residues < 0:
            raise ReferenceValidationError(
                "kinase_library.sequence_window.upstream_residues must be >= 0"
            )
        if sequence_window.downstream_residues < 0:
            raise ReferenceValidationError(
                "kinase_library.sequence_window.downstream_residues must be >= 0"
            )
        if not isinstance(sequence_window.central_residue_required, bool):
            raise ReferenceValidationError(
                "kinase_library.sequence_window.central_residue_required must be bool"
            )

    def _validate_matrices(
        self,
        *,
        matrices: tuple[KinaseLibraryMatrix, ...],
        sequence_window: SequenceWindowDefinition,
    ) -> None:
        if not matrices:
            raise ReferenceValidationError("kinase_library.matrices must not be empty")
        expected_positions = tuple(
            range(
                -int(sequence_window.upstream_residues),
                int(sequence_window.downstream_residues) + 1,
            )
        )
        seen_keys: set[tuple[str, str]] = set()
        for matrix in matrices:
            if not isinstance(matrix, KinaseLibraryMatrix):
                raise ReferenceValidationError(
                    "kinase_library.matrices must contain KinaseLibraryMatrix values"
                )
            key = (matrix.kinase, matrix.residue_class.value)
            if key in seen_keys:
                raise ReferenceValidationError(
                    "kinase_library.matrices contains duplicate kinase/residue_class "
                    f"entry: {matrix.kinase}/{matrix.residue_class.value}"
                )
            seen_keys.add(key)
            self._validate_score_table(
                matrix.score_table,
                expected_positions=expected_positions,
                context=f"kinase_library[{matrix.kinase}:{matrix.residue_class.value}]",
            )

    @staticmethod
    def _validate_score_table(
        score_table: pd.DataFrame,
        *,
        expected_positions: tuple[int, ...],
        context: str,
    ) -> None:
        if score_table.empty:
            raise ReferenceValidationError(f"{context}.score_table must be non-empty")
        observed_positions = tuple(int(position) for position in score_table.columns)
        missing_positions = [
            position
            for position in expected_positions
            if position not in observed_positions
        ]
        unexpected_positions = [
            position
            for position in observed_positions
            if position not in expected_positions
        ]
        if missing_positions:
            raise ReferenceValidationError(
                f"{context}.score_table is missing required positions: "
                f"{_format_positions(missing_positions)}"
            )
        if unexpected_positions:
            raise ReferenceValidationError(
                f"{context}.score_table contains positions outside sequence_window: "
                f"{_format_positions(unexpected_positions)}"
            )
        if observed_positions != expected_positions:
            raise ReferenceValidationError(
                f"{context}.score_table positions must be ordered as "
                f"{_format_positions(list(expected_positions))}"
            )
        if score_table.isna().to_numpy().any():
            raise ReferenceValidationError(
                f"{context}.score_table contains missing score values"
            )
        values = score_table.to_numpy(dtype=float, copy=False)
        for value in values.ravel().tolist():
            if not math.isfinite(float(value)):
                raise ReferenceValidationError(
                    f"{context}.score_table contains non-finite score values"
                )

    @staticmethod
    def _validate_provenance(
        provenance: KinaseLibraryResourceProvenance,
    ) -> None:
        if not isinstance(cast(object, provenance), KinaseLibraryResourceProvenance):
            raise ReferenceValidationError(
                "kinase_library.provenance must be KinaseLibraryResourceProvenance"
            )
        if provenance.source_type != "local":
            raise ReferenceValidationError(
                "kinase_library.provenance.source_type must be 'local'"
            )
        if not provenance.source_files:
            raise ReferenceValidationError(
                "kinase_library.provenance.source_files must not be empty"
            )
        if not provenance.table_fingerprints:
            raise ReferenceValidationError(
                "kinase_library.provenance.table_fingerprints must not be empty"
            )


def _format_positions(positions: list[int]) -> str:
    return ", ".join(str(position) for position in positions)


__all__ = ["KinaseLibraryResourceValidator"]

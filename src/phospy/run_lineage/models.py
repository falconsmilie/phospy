"""Cross-capability run-lineage provenance model."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import cast

from phospy.errors.input import PhosPyInputError
from phospy.preprocessing.provenance.models import (
    PreprocessingStageProvenance,
    _required_preprocessing_stage_tuple,  # pyright: ignore[reportPrivateUsage] - run lineage reuses the owner validator to preserve the accepted validation order and error text
)
from phospy.provenance.immutability import freeze_json_mapping, thaw_json_mapping
from phospy.provenance.models._shared import JsonValue
from phospy.provenance.models.environment import EnvironmentProvenance
from phospy.provenance.models.tables import (
    TableFingerprint,
    _required_table_fingerprint_tuple,  # pyright: ignore[reportPrivateUsage] - run lineage reuses the primitive-owner validator to preserve accepted validation behavior
)
from phospy.provenance.scientific_policy_models import ScientificPolicyRecord
from phospy.references.provenance.context import ReferenceContext
from phospy.references.provenance.models import (
    ReferenceProvenance,
    _require_run_provenance_reference_context_organism_coherence,  # pyright: ignore[reportPrivateUsage] - run lineage owns the composite check while References owns its exact coherence validator
)


@dataclass(frozen=True, slots=True)
class RunProvenance:
    """Machine-readable run provenance payload.

    Workflow-specific audit details that are not table-transforming stages, such
    as protein-aware preparation summaries, belong in `workflow_parameters`.
    `reference_context` records the biological reference context of the input
    dataset for downstream compatibility checks. `reference` remains the
    workflow reference resource provenance when a workflow consumes an explicit
    reference bundle.
    """

    environment: EnvironmentProvenance
    input_tables: tuple[TableFingerprint, ...]
    preprocessing_stages: tuple[PreprocessingStageProvenance, ...]
    reference: ReferenceProvenance | None
    workflow_name: str | None
    workflow_parameters: Mapping[str, object]
    random_state: int | None
    random_seed_policy: str | None
    output_tables: tuple[TableFingerprint, ...]
    scientific_policies: tuple[ScientificPolicyRecord, ...] = ()
    reference_context: ReferenceContext | None = None
    _frozen_workflow_parameters: Mapping[str, JsonValue] = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __getattribute__(self, name: str) -> object:
        if name == "workflow_parameters":
            try:
                frozen = object.__getattribute__(self, "_frozen_workflow_parameters")
            except AttributeError:
                return object.__getattribute__(self, name)
            return thaw_json_mapping(
                frozen,
                field_name="run_provenance.workflow_parameters",
            )
        return object.__getattribute__(self, name)

    def __post_init__(self) -> None:
        if not isinstance(  # pyright: ignore[reportUnnecessaryIsInstance] - runtime construction can bypass static annotations and must retain the accepted concrete-type error
            self.environment,
            EnvironmentProvenance,
        ):
            raise PhosPyInputError(
                "run_provenance.environment must be EnvironmentProvenance"
            )
        object.__setattr__(
            self,
            "input_tables",
            _required_table_fingerprint_tuple(
                self.input_tables,
                field_name="run_provenance.input_tables",
            ),
        )
        object.__setattr__(
            self,
            "preprocessing_stages",
            _required_preprocessing_stage_tuple(self.preprocessing_stages),
        )
        if self.reference is not None and not isinstance(  # pyright: ignore[reportUnnecessaryIsInstance] - runtime construction must retain the accepted concrete nested-type validation
            self.reference,
            ReferenceProvenance,
        ):
            raise PhosPyInputError(
                "run_provenance.reference must be ReferenceProvenance or None"
            )
        raw_workflow_parameters = object.__getattribute__(self, "workflow_parameters")
        frozen_workflow_parameters = freeze_json_mapping(
            raw_workflow_parameters,
            field_name="run_provenance.workflow_parameters",
        )
        object.__setattr__(
            self,
            "_frozen_workflow_parameters",
            frozen_workflow_parameters,
        )
        object.__setattr__(
            self,
            "workflow_parameters",
            frozen_workflow_parameters,
        )
        object.__setattr__(
            self,
            "output_tables",
            _required_table_fingerprint_tuple(
                self.output_tables,
                field_name="run_provenance.output_tables",
            ),
        )
        object.__setattr__(
            self,
            "scientific_policies",
            _required_scientific_policy_tuple(self.scientific_policies),
        )
        _require_run_provenance_reference_context_organism_coherence(
            reference=self.reference,
            reference_context=self.reference_context,
        )


def _required_scientific_policy_tuple(
    values: object,
) -> tuple[ScientificPolicyRecord, ...]:
    field_name = "run_provenance.scientific_policies"
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    policies: list[ScientificPolicyRecord] = []
    for policy in cast(Sequence[object], values):
        if not isinstance(policy, ScientificPolicyRecord):
            raise PhosPyInputError(
                f"{field_name} must contain only ScientificPolicyRecord values"
            )
        policies.append(policy)
    return tuple(policies)

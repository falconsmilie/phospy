from __future__ import annotations

from collections.abc import Mapping
from typing import cast

import pytest

from phospy.core.errors.input import PhosPyInputError
from phospy.core.errors.validation import ReferenceValidationError
from phospy.preprocessing.provenance.models import PreprocessingStageProvenance
from phospy.provenance.models.environment import EnvironmentProvenance
from phospy.provenance.models.tables import TableFingerprint
from phospy.references.provenance.context import ReferenceContext
from phospy.references.provenance.models import ReferenceProvenance
from phospy.run_lineage.models import RunProvenance
from phospy.taxonomy.organisms import Organism


def test_run_lineage_owner_contracts_preserve_sequence_order() -> None:
    first_fingerprint = _fingerprint("first")
    second_fingerprint = _fingerprint("second")
    first_stage = _stage("first")
    second_stage = _stage("second")

    provenance = _run_provenance(
        input_tables=[first_fingerprint, second_fingerprint],
        preprocessing_stages=[first_stage, second_stage],
        reference=_reference(Organism.RAT),
        output_tables=[second_fingerprint, first_fingerprint],
        reference_context=_context(Organism.RAT),
    )

    assert provenance.input_tables == (first_fingerprint, second_fingerprint)
    assert provenance.preprocessing_stages == (first_stage, second_stage)
    assert provenance.output_tables == (second_fingerprint, first_fingerprint)


def test_run_lineage_rejects_invalid_fingerprint_members_with_contract_message() -> (
    None
):
    with pytest.raises(PhosPyInputError) as exc_info:
        _run_provenance(input_tables=[object()])

    assert type(exc_info.value) is PhosPyInputError
    assert str(exc_info.value) == (
        "run_provenance.input_tables must contain only TableFingerprint values"
    )


def test_run_lineage_rejects_invalid_stage_members_with_contract_message() -> None:
    with pytest.raises(PhosPyInputError) as exc_info:
        _run_provenance(preprocessing_stages=[object()])

    assert type(exc_info.value) is PhosPyInputError
    assert str(exc_info.value) == (
        "run_provenance.preprocessing_stages must contain only "
        "PreprocessingStageProvenance values"
    )


def test_run_lineage_multiple_invalid_fields_preserve_validation_order() -> None:
    with pytest.raises(PhosPyInputError) as input_exc:
        _run_provenance(
            input_tables="invalid",
            preprocessing_stages="invalid",
            reference=object(),
            output_tables="invalid",
        )
    assert str(input_exc.value) == "run_provenance.input_tables must be a sequence"

    with pytest.raises(PhosPyInputError) as stage_exc:
        _run_provenance(
            preprocessing_stages="invalid",
            reference=object(),
            output_tables="invalid",
        )
    assert str(stage_exc.value) == (
        "run_provenance.preprocessing_stages must be a sequence"
    )

    with pytest.raises(PhosPyInputError) as reference_exc:
        _run_provenance(reference=object(), output_tables="invalid")
    assert str(reference_exc.value) == (
        "run_provenance.reference must be ReferenceProvenance or None"
    )

    with pytest.raises(PhosPyInputError) as output_exc:
        _run_provenance(
            output_tables="invalid",
            reference=_reference(Organism.RAT),
            reference_context=_context(Organism.HUMAN),
        )
    assert str(output_exc.value) == "run_provenance.output_tables must be a sequence"


def test_run_lineage_reference_coherence_preserves_reference_error_contract() -> None:
    with pytest.raises(ReferenceValidationError) as exc_info:
        _run_provenance(
            reference=_reference(Organism.RAT),
            reference_context=_context(Organism.HUMAN),
        )

    assert type(exc_info.value) is ReferenceValidationError
    assert str(exc_info.value) == (
        "Run provenance reference-context organism mismatch: "
        "run_provenance.reference.organism='rat'; "
        "run_provenance.reference_context.organism='human' "
        "resolved_to='human'"
    )


def _run_provenance(
    *,
    input_tables: object = (),
    preprocessing_stages: object = (),
    reference: object = None,
    workflow_parameters: object = None,
    output_tables: object = (),
    reference_context: object = None,
) -> RunProvenance:
    parameters = {} if workflow_parameters is None else workflow_parameters
    return RunProvenance(
        environment=EnvironmentProvenance(
            package_name="phospy",
            package_version="test",
            python_version="test",
            dependency_versions={},
        ),
        input_tables=cast(tuple[TableFingerprint, ...], input_tables),
        preprocessing_stages=cast(
            tuple[PreprocessingStageProvenance, ...],
            preprocessing_stages,
        ),
        reference=cast(ReferenceProvenance | None, reference),
        workflow_name="unit-test",
        workflow_parameters=cast(Mapping[str, object], parameters),
        random_state=None,
        random_seed_policy=None,
        output_tables=cast(tuple[TableFingerprint, ...], output_tables),
        reference_context=cast(ReferenceContext | None, reference_context),
    )


def _fingerprint(name: str) -> TableFingerprint:
    return TableFingerprint(
        name=name,
        rows=1,
        columns=1,
        index_name="site_key",
        column_names=("sample",),
        dtypes=("float64",),
        exact_hash_algorithm="sha256",
        exact_hash_value=f"{name}-exact",
        tolerance_hash_algorithm="sha256",
        tolerance_hash_value=f"{name}-tolerance",
    )


def _stage(name: str) -> PreprocessingStageProvenance:
    return PreprocessingStageProvenance(
        stage=name,
        operation="identity",
        parameters={},
        input_shape=(1, 1),
        output_shape=(1, 1),
        input_hash=f"{name}-input",
        output_hash=f"{name}-output",
        dropped_row_ids=(),
        dropped_row_count=0,
    )


def _reference(organism: Organism) -> ReferenceProvenance:
    return ReferenceProvenance(
        source_type="unit-test",
        organism=organism,
        bundle_id=None,
        table_fingerprints=(),
    )


def _context(organism: Organism) -> ReferenceContext:
    return ReferenceContext(
        organism=organism,
        protein_namespace="gene_symbol",
        source_name="unit-test",
        source_version="v1",
        proteome_version=None,
        reference_table_sha256="a" * 64,
    )

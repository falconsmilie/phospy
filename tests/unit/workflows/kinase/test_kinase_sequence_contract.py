from __future__ import annotations

from contextlib import nullcontext

import pandas as pd
import pytest

from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api.configs.kinase import (
    KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
)
from phospy.errors import WorkflowValidationError
from phospy.provenance.hashing import fingerprint_table
from phospy.provenance.models import KinaseLibraryResourceProvenance
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.prediction.motif_scoring.models import AMINO_ACIDS
from phospy.science.references.kinase_library import (
    KinaseLibraryMatrix,
    KinaseLibraryResidueClass,
    KinaseLibraryResource,
)
from phospy.science.references.models import SequenceWindowDefinition
from phospy.workflows.kinase.validator import (
    KinaseWorkflowValidator,
    _selected_explicit_reference_sequence_context,
)
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)
from tests.support.unsafe_dataset_states import (
    unsafe_drop_dataset_site_metadata_columns,
    unsafe_set_dataset_site_metadata_columns,
)

_PANDAS_COPY_ON_WRITE_CAN_BE_DISABLED = (
    int(pd.__version__.split(".", maxsplit=1)[0]) < 3
)


def _window(residue: str, *, flank: int = 7) -> str:
    return ("A" * flank) + residue + ("A" * flank)


def _dataset(
    *,
    sequences: tuple[str, str] = (_window("S"), _window("T")),
) -> AnalysisReadyPhosphoDataset:
    display_ids = ["GENE1;S10;", "GENE2;T20;"]
    site_index = site_key_index_from_display_ids(display_ids)
    return trusted_analysis_ready_dataset_from_tables(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0], "sample_b": [1.5, 2.5]},
            index=site_index.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_index.astype(str).tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_index),
                "gene_symbol": ["GENE1", "GENE2"],
                "site": ["S10", "T20"],
                "protein_id": ["GENE1", "GENE2"],
                "site_sequence": list(sequences),
            },
            index=site_index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references(
    *,
    sequences: tuple[str, str] = (_window("S"), _window("T")),
) -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["KLIB_ST", "KLIB_ST"],
                "substrate_site": ["GENE1;S10;", "GENE2;T20;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": list(sequences)},
            index=pd.Index(["GENE1;S10;", "GENE2;T20;"], name="site_id"),
        ),
    )


def _kinase_library_resource() -> KinaseLibraryResource:
    positions = tuple(range(-7, 8))
    score_table = pd.DataFrame(
        0.0,
        index=pd.Index(AMINO_ACIDS, name="amino_acid"),
        columns=pd.Index(positions, name="position"),
    )
    score_table.loc["S", 0] = 2.0
    score_table.loc["T", 0] = 1.0
    matrix = KinaseLibraryMatrix(
        kinase="KLIB_ST",
        residue_class=KinaseLibraryResidueClass.SER_THR,
        score_table=score_table,
    )
    sequence_window = SequenceWindowDefinition(
        upstream_residues=7,
        downstream_residues=7,
        central_residue_required=True,
    )
    provenance = KinaseLibraryResourceProvenance(
        source_type="local",
        source_name="synthetic_kinase_library",
        source_version="test",
        license="test-only",
        score_scale="synthetic_raw_position_sum",
        organisms=(Organism.RAT.value,),
        sequence_window=sequence_window.to_payload(),
        source_files={"kinase_library": {"path": "synthetic"}},
        table_fingerprints=(
            fingerprint_table(
                score_table,
                name="references.kinase_library.score_table.klib_st.ser_thr",
            ),
        ),
        manifest={
            "resource_type": "kinase_library",
            "source_name": "synthetic_kinase_library",
            "source_version": "test",
            "score_scale": "synthetic_raw_position_sum",
            "organisms": (Organism.RAT.value,),
            "sequence_window": sequence_window.to_payload(),
        },
    )
    return KinaseLibraryResource(
        matrices=(matrix,),
        source_name="synthetic_kinase_library",
        source_version="test",
        score_scale="synthetic_raw_position_sum",
        sequence_window=sequence_window,
        organisms=(Organism.RAT.value,),
        license="test-only",
        provenance=provenance,
    )


def _request(
    *,
    dataset: AnalysisReadyPhosphoDataset | None = None,
    references: ReferenceBundle | None = None,
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=_dataset() if dataset is None else dataset,
        references=_references() if references is None else references,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
            mode="deterministic_ranking",
        ),
        activity_config=None,
        kinase_library_resource=_kinase_library_resource(),
    )


def test_valid_centered_15aa_window_passes_for_kinase_library_mode() -> None:
    request = _request()

    validated = KinaseWorkflowValidator().run(request)

    assert validated is request


def test_missing_site_sequence_fails_before_kinase_motif_scoring() -> None:
    dataset = _dataset()
    unsafe_drop_dataset_site_metadata_columns(dataset, "site_sequence")

    with pytest.raises(WorkflowValidationError) as exc_info:
        KinaseWorkflowValidator().run(_request(dataset=dataset))

    assert "site_sequence" in str(exc_info.value)
    assert "kinase workflow request" in str(exc_info.value)


def test_empty_site_sequence_fails_before_kinase_motif_scoring() -> None:
    dataset = _dataset()
    unsafe_set_dataset_site_metadata_columns(
        dataset,
        {"site_sequence": ["", _window("T")]},
    )

    with pytest.raises(WorkflowValidationError) as exc_info:
        KinaseWorkflowValidator().run(_request(dataset=dataset))

    message = str(exc_info.value)
    assert "missing or blank" in message
    assert str(dataset.site_metadata.index[0]) in message


@pytest.mark.parametrize(
    ("bad_sequence", "expected"),
    [
        (_window("S", flank=15), "expected_length=15"),
        ("AAAAAAAAAAAAAAA", "expected_center_residues=S/T/Y"),
        ("AAAAAAA*AAAAAAA", "unsupported_characters='*'"),
        (_window("s"), "lowercase characters are not allowed"),
        ("__AAAAASAAAAAAA", "terminal padding is not allowed"),
    ],
)
def test_kinase_library_mode_rejects_invalid_selected_sequence_context(
    bad_sequence: str,
    expected: str,
) -> None:
    references = _references(sequences=(bad_sequence, _window("T")))
    site_key = str(_dataset().site_metadata.index[0])

    with pytest.raises(WorkflowValidationError) as exc_info:
        KinaseWorkflowValidator().run(_request(references=references))

    message = str(exc_info.value)
    assert "kinase_library_contextual_motif" in message
    assert "expected_length=15" in message
    assert "expected_center_index=7" in message
    assert site_key in message
    assert expected in message


def test_base_dataset_sequence_can_pass_construction_but_fail_kinase_fixed_window() -> (
    None
):
    long_centered = _window("S", flank=15)
    dataset = _dataset(sequences=(long_centered, _window("T")))
    references = _references(sequences=(long_centered, _window("T")))

    assert dataset.site_metadata.loc[
        dataset.site_metadata.index[0], "site_sequence"
    ] == (long_centered)

    with pytest.raises(WorkflowValidationError) as exc_info:
        KinaseWorkflowValidator().run(_request(dataset=dataset, references=references))

    message = str(exc_info.value)
    assert "workflow-specific sequence context contract failed" in message
    assert "expected_length=15" in message
    assert "observed_length=31" in message


@pytest.mark.parametrize(
    "copy_on_write",
    (
        [
            pytest.param(False, id="copy-on-write-disabled"),
            pytest.param(True, id="copy-on-write-enabled"),
        ]
        if _PANDAS_COPY_ON_WRITE_CAN_BE_DISABLED
        else [pytest.param(None, id="copy-on-write-always-enabled")]
    ),
)
def test_selected_explicit_reference_sequence_context_does_not_mutate_input_frame(
    copy_on_write: bool | None,
) -> None:
    dataset = _dataset(sequences=(_window("S"), _window("T")))
    site_metadata = dataset.site_metadata
    original_site_metadata = site_metadata.copy(deep=True)
    references = _references(sequences=(_window("T"), _window("T")))
    first_site_key = site_metadata.index[0]
    copy_mode = (
        nullcontext()
        if copy_on_write is None
        else pd.option_context("mode.copy_on_write", copy_on_write)
    )

    with copy_mode:
        selected, source_by_site = _selected_explicit_reference_sequence_context(
            dataset=dataset,
            site_metadata=site_metadata,
            references=references,
            conflict_policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
        )

    pd.testing.assert_frame_equal(site_metadata, original_site_metadata)
    assert selected.loc[first_site_key, "site_sequence"] == _window("T")
    assert site_metadata.loc[first_site_key, "site_sequence"] == _window("S")
    assert source_by_site[str(first_site_key)] == "reference"

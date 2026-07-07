from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
import pytest

from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.api.configs.kinase import (
    KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
    KINASE_SCORING_MODES,
    KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY,
    normalize_kinase_scoring_mode,
)
from phospy.errors import WorkflowBoundaryError, WorkflowValidationError
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
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.kinase.scoring_mode_contracts import (
    KINASE_ANALYSIS_READY_SITE_METADATA_COLUMNS,
    KINASE_SCORING_MODE_INPUT_CONTRACTS,
    kinase_scoring_mode_input_contract,
)
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
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

_DEFAULT_KINASE_LIBRARY_RESOURCE = object()


@dataclass(frozen=True, slots=True)
class ExpectedKinaseScoringModeContract:
    scoring_mode: str
    requires_site_sequence: bool
    requires_centered_sequence_context: bool
    requires_localisation_probability: bool
    requires_substrate_reference_overlap: bool
    requires_kinase_library_resource: bool
    requires_profile_construction: bool


_EXPECTED_REQUIRED_DATASET_COLUMNS = KINASE_ANALYSIS_READY_SITE_METADATA_COLUMNS
_EXPECTED_CONTRACTS = (
    ExpectedKinaseScoringModeContract(
        scoring_mode=KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        requires_site_sequence=True,
        requires_centered_sequence_context=True,
        requires_localisation_probability=False,
        requires_substrate_reference_overlap=True,
        requires_kinase_library_resource=False,
        requires_profile_construction=True,
    ),
    ExpectedKinaseScoringModeContract(
        scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
        requires_site_sequence=True,
        requires_centered_sequence_context=True,
        requires_localisation_probability=False,
        requires_substrate_reference_overlap=True,
        requires_kinase_library_resource=True,
        requires_profile_construction=True,
    ),
    ExpectedKinaseScoringModeContract(
        scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
        requires_site_sequence=True,
        requires_centered_sequence_context=True,
        requires_localisation_probability=False,
        requires_substrate_reference_overlap=False,
        requires_kinase_library_resource=True,
        requires_profile_construction=False,
    ),
    ExpectedKinaseScoringModeContract(
        scoring_mode=KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
        requires_site_sequence=True,
        requires_centered_sequence_context=True,
        requires_localisation_probability=False,
        requires_substrate_reference_overlap=True,
        requires_kinase_library_resource=True,
        requires_profile_construction=True,
    ),
)


@pytest.mark.parametrize(
    "expected",
    _EXPECTED_CONTRACTS,
    ids=lambda expected: expected.scoring_mode,
)
def test_supported_kinase_scoring_mode_input_contract_is_explicit(
    expected: ExpectedKinaseScoringModeContract,
) -> None:
    contract = kinase_scoring_mode_input_contract(expected.scoring_mode)

    assert contract.scoring_mode == expected.scoring_mode
    assert contract.required_dataset_columns == _EXPECTED_REQUIRED_DATASET_COLUMNS
    assert contract.requires_site_sequence is expected.requires_site_sequence
    assert (
        contract.requires_centered_sequence_context
        is expected.requires_centered_sequence_context
    )
    assert (
        contract.requires_localisation_probability
        is expected.requires_localisation_probability
    )
    assert (
        contract.requires_substrate_reference_overlap
        is expected.requires_substrate_reference_overlap
    )
    assert (
        contract.requires_kinase_library_resource
        is expected.requires_kinase_library_resource
    )
    assert contract.requires_profile_construction is (
        expected.requires_profile_construction
    )


def test_every_canonical_kinase_scoring_mode_has_input_contract() -> None:
    canonical_config_modes = {
        normalize_kinase_scoring_mode(scoring_mode)
        for scoring_mode in KINASE_SCORING_MODES
    }
    contract_modes = set(KINASE_SCORING_MODE_INPUT_CONTRACTS)

    assert canonical_config_modes - contract_modes == set(), (
        "add missing modes to KINASE_SCORING_MODE_INPUT_CONTRACTS"
    )
    assert contract_modes - canonical_config_modes == set()


def test_deprecated_kinase_library_motif_alias_uses_contextual_contract() -> None:
    contract = kinase_scoring_mode_input_contract(
        KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF
    )

    assert contract.scoring_mode == KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF
    assert contract.requires_substrate_reference_overlap is True
    assert contract.requires_profile_construction is True


def test_public_kinase_library_resource_mode_constant_matches_mode_contracts() -> None:
    resource_modes = {
        scoring_mode
        for scoring_mode in KINASE_SCORING_MODES
        if kinase_scoring_mode_input_contract(
            scoring_mode
        ).requires_kinase_library_resource
    }

    assert resource_modes == set(KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY)


@pytest.mark.parametrize(
    "scoring_mode",
    [
        KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
        KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
        KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    ],
)
def test_kinase_validator_requires_site_sequence_for_every_scoring_mode(
    scoring_mode: str,
) -> None:
    dataset = _dataset()
    unsafe_drop_dataset_site_metadata_columns(dataset, "site_sequence")

    with pytest.raises(WorkflowValidationError, match="site_sequence"):
        KinaseWorkflowValidator().run(
            _request(dataset=dataset, scoring_mode=scoring_mode)
        )


@pytest.mark.parametrize(
    "scoring_mode",
    [
        KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
        KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
        KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    ],
)
def test_kinase_validator_requires_centered_sequence_context_for_every_scoring_mode(
    scoring_mode: str,
) -> None:
    dataset = _dataset()
    unsafe_set_dataset_site_metadata_columns(
        dataset,
        {"site_sequence": ["AAAAAAAAAAAAAAA", _window("T")]},
    )

    with pytest.raises(WorkflowValidationError, match="center residue"):
        KinaseWorkflowValidator().run(
            _request(dataset=dataset, scoring_mode=scoring_mode)
        )


@pytest.mark.parametrize(
    "scoring_mode",
    [
        KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
        KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
        KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    ],
)
def test_scoring_modes_do_not_require_localisation_probability_by_default(
    scoring_mode: str,
) -> None:
    request = _request(scoring_mode=scoring_mode)

    validated = KinaseWorkflowValidator().run(request)

    assert validated is request


@pytest.mark.parametrize(
    "scoring_mode",
    [
        KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
        KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
        KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    ],
)
def test_kinase_library_modes_require_resource_before_execution(
    scoring_mode: str,
) -> None:
    with pytest.raises(WorkflowValidationError, match="kinase_library_resource"):
        KinaseWorkflowValidator().run(
            _request(scoring_mode=scoring_mode, kinase_library_resource=None)
        )


@pytest.mark.parametrize(
    "scoring_mode",
    [
        KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
        KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    ],
)
def test_profile_modes_require_substrate_reference_overlap(
    scoring_mode: str,
) -> None:
    request = _request(
        references=_references(mapped_display_ids=("UNMATCHED;S99;",)),
        scoring_mode=scoring_mode,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflowInterpreter().run(request)

    assert exc_info.value.seam == "kinase.interpreter.reference_coverage"


def test_kinase_library_motif_only_does_not_require_substrate_reference_overlap() -> (
    None
):
    request = _request(
        references=_references(mapped_display_ids=("UNMATCHED;S99;",)),
        scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
    )

    resolved = KinaseWorkflowInterpreter().run(request)

    assert resolved.kinase_substrate_map.empty
    assert resolved.execution_config.scoring_mode == (
        KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY
    )


def _window(residue: str) -> str:
    return ("A" * 7) + residue + ("A" * 7)


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["GENE1;S10;", "GENE2;T20;"]
    site_index = site_key_index_from_display_ids(display_ids)
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [1.5, 2.5]},
        index=site_index.copy(),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": ["GENE1", "GENE2"],
            "site": ["S10", "T20"],
            "protein_id": ["GENE1", "GENE2"],
            "site_sequence": [_window("S"), _window("T")],
        },
        index=site_index.copy(),
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references(
    *,
    mapped_display_ids: tuple[str, ...] = ("GENE1;S10;", "GENE2;T20;"),
) -> ReferenceBundle:
    sequence_display_ids = tuple(
        dict.fromkeys(("GENE1;S10;", "GENE2;T20;", *mapped_display_ids))
    )
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["KPROFILE" for _ in mapped_display_ids],
                "substrate_site": list(mapped_display_ids),
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    _window(_residue_from_display_id(display_id))
                    for display_id in sequence_display_ids
                ]
            },
            index=pd.Index(list(sequence_display_ids), name="site_id"),
        ),
    )


def _residue_from_display_id(display_id: str) -> str:
    return display_id.split(";", 2)[1][0]


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
    scoring_mode: str,
    kinase_library_resource: object = _DEFAULT_KINASE_LIBRARY_RESOURCE,
) -> KinaseWorkflowRequest:
    resource = kinase_library_resource
    if resource is _DEFAULT_KINASE_LIBRARY_RESOURCE:
        resource = (
            None
            if scoring_mode == KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED
            else _kinase_library_resource()
        )
    if resource is not None and not isinstance(resource, KinaseLibraryResource):
        raise TypeError(
            "test fixture kinase_library_resource must be a resource or None"
        )
    return KinaseWorkflowRequest(
        dataset=dataset or _dataset(),
        references=references or _references(),
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            scoring_mode=scoring_mode,
        ),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
            mode="deterministic_ranking",
        ),
        activity_config=None,
        kinase_library_resource=resource,
    )

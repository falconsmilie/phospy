from __future__ import annotations

import pandas as pd
import pytest

from phospy.api.configs import KinasePredictionConfig, KinaseScoringConfig
from phospy.api.requests import KinaseWorkflowRequest
from phospy.errors import (
    PhosPyInputError,
    PhosPyWorkflowError,
    ReferenceResolutionError,
)
from phospy.science.datasets.builders.reader import DatasetInputReader
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.references.models import Organism, ReferencePreset
from phospy.science.references.resolution import ReferenceResolver
from phospy.validation.references.compatibility import ReferenceCompatibilityValidator
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import site_key_index_from_display_ids


class _BrokenSourceValidator:
    def run(self, value: object, *, field_name: str) -> object:
        return 123


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_id = "MAPK14;Y182;"
    index = site_key_index_from_display_ids([display_id])
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=index),
        site_metadata=pd.DataFrame(
            {
                "site_key": index.astype(str).tolist(),
                "display_id": [display_id],
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def test_dataset_input_reader_rejects_validator_contract_breach() -> None:
    reader = DatasetInputReader(source_validator=_BrokenSourceValidator())  # type: ignore[arg-type]
    with pytest.raises(
        PhosPyInputError,
        match="source validator produced unsupported value type int",
    ):
        reader.run(pd.DataFrame({"sample_a": [1.0]}), field_name="phospho")


def test_reference_resolver_rejects_non_reference_input_types() -> None:
    with pytest.raises(
        ReferenceResolutionError,
        match="reference input must be a ReferencePreset or ReferenceBundle; got str",
    ):
        ReferenceResolver().run("rat", dataset_organism=Organism.RAT)  # type: ignore[arg-type]


def test_reference_compatibility_auto_resolution_stays_explicit_without_assert(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    validator = ReferenceCompatibilityValidator()

    def _skip_run(
        reference_input: object, *, dataset_organism: Organism | None
    ) -> None:
        return None

    monkeypatch.setattr(validator, "run", _skip_run)
    with pytest.raises(
        ReferenceResolutionError,
        match="ReferencePreset.AUTO requires dataset.organism",
    ):
        validator.resolve_preset_organism(
            preset=ReferencePreset.AUTO,
            dataset_organism=None,
        )


def test_kinase_interpreter_rejects_invalid_overlap_summary_state() -> None:
    interpreter = KinaseWorkflowInterpreter()
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
    )
    with pytest.raises(
        PhosPyWorkflowError,
        match="overlap_counts\\['per_kinase_quantified'\\] to be a pandas Series",
    ):
        # Intentional private-seam guard: this validates a defensive contract on the
        # internal overlap summary shape, which cannot be reached via public requests.
        interpreter._validate_eligible_kinases(
            overlap_counts={"per_kinase_quantified": []},  # type: ignore[arg-type]
            request=request,
        )

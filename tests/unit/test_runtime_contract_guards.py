from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.api.configs import KinasePredictionConfig, KinaseScoringConfig
from phospy.api.requests import KinaseWorkflowRequest
from phospy.errors import (
    PhosPyInputError,
    PhosPyWorkflowError,
    ReferenceCompatibilityError,
    ReferenceResolutionError,
)
from phospy.science.datasets.builders.reader import DatasetInputReader
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.science.references.resolution import ReferenceResolver
from phospy.validation.references.compatibility import ReferenceCompatibilityValidator
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.kinase.resolved_validator import (
    ResolvedKinaseEligibilityValidator,
    ResolvedKinaseInputs,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)


class _BrokenSourceValidator:
    def run(self, value: object, *, field_name: str) -> object:
        return 123


class _DatasetPathReader:
    def run(self, path: Path, *, field_name: str) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "path": [str(path)],
                "field_name": [field_name],
            }
        )


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_id = "MAPK14;Y182;"
    index = site_key_index_from_display_ids([display_id])
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=index),
        site_metadata=pd.DataFrame(
            {
                "site_key": index.astype(str).tolist(),
                "display_id": [display_id],
                **site_key_context_columns(index),
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


def test_dataset_input_reader_requires_injected_path_reader_for_paths() -> None:
    reader = DatasetInputReader()

    with pytest.raises(
        PhosPyInputError,
        match="path inputs require an injected DatasetPathTableReader",
    ):
        reader.run("phospho.csv", field_name="phospho")


def test_dataset_input_reader_delegates_path_sources_to_injected_reader() -> None:
    reader = DatasetInputReader(path_reader=_DatasetPathReader())

    resolved = reader.run("phospho.csv", field_name="phospho")

    assert resolved.to_dict(orient="records") == [
        {"path": "phospho.csv", "field_name": "phospho"}
    ]


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


def test_reference_compatibility_auto_canonicalizes_dataset_organism_alias() -> None:
    resolved = ReferenceCompatibilityValidator().resolve_preset_organism(
        preset=ReferencePreset.AUTO,
        dataset_organism="Rattus norvegicus",
    )

    assert resolved is Organism.RAT


def test_reference_compatibility_explicit_bundle_uses_canonical_dataset_organism() -> (
    None
):
    bundle = ReferenceBundle(
        organism=Organism.HUMAN,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["K"], "substrate_site": ["MAPK14;Y182;"]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["AAAAAYAAAAA"]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
    )
    validator = ReferenceCompatibilityValidator()

    validator.run(bundle, dataset_organism="Homo sapiens")
    with pytest.raises(
        ReferenceCompatibilityError,
        match="references\\.organism must match dataset\\.organism",
    ):
        validator.run(bundle, dataset_organism="rat")


def test_resolved_kinase_validator_rejects_invalid_overlap_summary_state() -> None:
    validator = ResolvedKinaseEligibilityValidator()
    dataset = _dataset()
    request = KinaseWorkflowRequest(
        dataset=dataset,
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
        validator._validate_eligible_kinases(
            overlap_counts={"per_kinase_quantified": []},  # type: ignore[arg-type]
            resolved_inputs=ResolvedKinaseInputs(
                dataset=dataset,
                dataset_phospho=dataset.phospho,
                references=ReferenceBundle(
                    organism=Organism.RAT,
                    kinase_substrate_map=pd.DataFrame(
                        {"kinase": ["K"], "substrate_site": ["MAPK14;Y182;"]}
                    ),
                    site_sequences=pd.DataFrame(
                        {"site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"]},
                        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
                    ),
                ),
                kinase_substrate_map=pd.DataFrame(),
                site_sequences=pd.DataFrame(),
                site_identity_map=pd.DataFrame(),
                scoring_site_index=pd.Index([], name=dataset.phospho.index.name),
                activity_phospho_matrix=dataset.phospho.iloc[0:0],
                execution_config=KinaseWorkflowInterpreter._resolve_execution_config(
                    request
                ),
                reference_site_count=0,
            ),
        )

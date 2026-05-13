from __future__ import annotations

from pathlib import Path

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SampleDesignRecord,
    SignalomeClusteringConfig,
    SignalomeConfig,
    SignalomeOutputConfig,
    SignalomePerformanceConfig,
    SignalomeScientificConfig,
    SignalomeValidationConfig,
    SignalomeWorkflowRequest,
    TechnicalReplicatePolicy,
)
from phospy.api.results import KinasePredictionResult, KinaseScoringResult

ROOT = Path(__file__).resolve().parents[2]
API_DOCS_DIR = ROOT / "docs" / "api"
DIFFERENTIAL_DOC = API_DOCS_DIR / "differential-workflow.md"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _tiny_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    phospho = pd.DataFrame(
        {
            "A_1": [100.0, 70.0],
            "A_2": [101.0, 72.0],
            "B_1": [120.0, 80.0],
            "B_2": [118.0, 82.0],
        },
        index=["TSC2;S939;", "GSK3B;S9;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["TSC2", "GSK3B"],
            "site": ["S939", "S9"],
            "site_sequence": [
                "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
                "ATMSGRPRTTSFAESCKPVQQPSAFGQAAAL",
            ],
            "protein_id": ["TSC2", "GSK3B"],
            "localisation_confidence": [0.95, 0.92],
        },
        index=phospho.index.copy(),
    )
    return phospho, site_metadata


def _build_dataset():
    phospho, site_metadata = _tiny_inputs()
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )


def _minimal_kinase_result() -> KinaseWorkflowResult:
    dataset = _build_dataset()
    site_ids = dataset.phospho.index.astype(str).tolist()
    index = pd.Index(site_ids, name="site_id")
    kinases = pd.Index(["K1"], name="kinase")
    score_matrix = pd.DataFrame([[0.9], [0.8]], index=index, columns=kinases)
    prediction_matrix = pd.DataFrame([[0.7], [0.6]], index=index, columns=kinases)
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K1"],
                "substrate_site": site_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31, "C" * 31]},
            index=index.copy(),
        ),
    )
    return KinaseWorkflowResult(
        dataset=dataset,
        references=references,
        scoring_result=KinaseScoringResult(profile_scores=score_matrix),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


def test_api_docs_public_imports_are_valid() -> None:
    namespace: dict[str, object] = {}

    exec(
        """from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
    KinaseWorkflow,
    SignalomeWorkflow,
)""",
        namespace,
    )
    exec(
        """from phospy.api import (
    DatasetBuildRequest,
    ExperimentalDesign,
    Contrast,
    SampleDesignRecord,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflowRequest,
    UnsupportedInputFormatError,
    WorkflowValidationError,
)""",
        namespace,
    )

    assert "AnalysisReadyDatasetBuilder" in namespace
    assert "DifferentialAnalysisWorkflow" in namespace
    assert "DatasetBuildRequest" in namespace
    assert "DifferentialAnalysisRequest" in namespace
    assert "SignalomeWorkflowRequest" in namespace


def test_api_docs_dataset_build_request_example_is_constructible() -> None:
    phospho, site_metadata = _tiny_inputs()
    preprocessing = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=1.0,
        ),
        normalisation=DatasetNormalisationConfig(policy="median_center"),
        missing_data=DatasetMissingDataConfig(policy="forbid"),
    )

    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        preprocessing_config=preprocessing,
    )

    assert request.organism is Organism.RAT
    assert request.preprocessing_config.normalisation.policy == "median_center"


def test_api_docs_differential_request_example_is_constructible() -> None:
    dataset = _build_dataset()
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
            ),
        )
    )
    contrasts = (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )
    request = DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=contrasts,
        config=DifferentialAnalysisConfig(
            technical_replicate_policy=TechnicalReplicatePolicy.MEAN,
            minimum_condition_replicates=2,
        ),
    )

    assert request.design.samples[0].sample_id == "A_1"
    assert request.contrasts[0].name == "B_vs_A"
    assert request.config.minimum_condition_replicates == 2


def test_api_docs_kinase_request_example_is_constructible() -> None:
    dataset = _build_dataset()
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            include_diagnostic_scoring_tables=False,
            profile_missing_value_strategy="strict",
        ),
        prediction_config=KinasePredictionConfig(
            mode="deterministic_ranking",
            top_k=30,
            deterministic_max_selected_kinases=10,
        ),
        activity_config=KinaseActivityConfig(
            enabled=True,
            method="simplified_weighted_substrate_activity",
            threshold=0.6,
            min_substrates=3,
            top_n_substrates=20,
        ),
        site_sequence_conflict_policy="prefer_reference",
    )

    assert request.references is ReferencePreset.AUTO
    assert request.prediction_config.mode == "deterministic_ranking"
    assert request.site_sequence_conflict_policy == "prefer_reference"


def test_api_docs_signalome_request_example_is_constructible() -> None:
    kinase_result = _minimal_kinase_result()
    config = SignalomeConfig(
        scientific=SignalomeScientificConfig(
            substrate_support_cutoff=0.5,
            assignment_policy="cutoff_binary",
        ),
        clustering=SignalomeClusteringConfig(
            module_count=None,
            module_selection_primary_correlation_threshold=0.5,
            module_selection_fallback_correlation_threshold=0.1,
            module_selection_max_clusters=10,
            candidate_scoring_policy="full",
            clustering_engine="scipy_hierarchical",
        ),
        validation=SignalomeValidationConfig(
            score_preconditioning_policy="error_on_drop",
            allow_mixed_total_protein_quantitative_meaning=False,
        ),
        output=SignalomeOutputConfig(
            network_correlation_threshold=0.5,
            network_policy="signed",
        ),
        performance=SignalomePerformanceConfig(
            max_exact_tree_sites=2000,
            max_full_candidate_scoring_sites=2000,
        ),
    )
    request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=config,
    )

    assert request.kinase_result.dataset.organism is Organism.RAT
    assert request.config.output.network_policy == "signed"


def test_api_docs_differential_import_route_uses_supported_public_path() -> None:
    source = _read(DIFFERENTIAL_DOC)

    assert "from phospy import DifferentialAnalysisWorkflow" in source
    assert "from phospy.api import (" in source
    assert "DifferentialAnalysisRequest," in source
    assert "`from phospy import DifferentialAnalysis` and" in source
    assert "`from phospy.api import DifferentialAnalysis` are not supported" in source

from __future__ import annotations

from pathlib import Path

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    Contrast,
    DatasetBatchCorrectionConfig,
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetLocalisationConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    EnrichmentConfig,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    ExperimentalDesign,
    GeneSetCollection,
    IntensityScaleKind,
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
)
from phospy.api.results import KinasePredictionResult, KinaseScoringResult

ROOT = Path(__file__).resolve().parents[2]
README = ROOT / "README.md"
API_DOCS_DIR = ROOT / "docs" / "api"
WORKFLOW_DOCS_DIR = API_DOCS_DIR / "workflows"
DATASET_BUILD_DOC = API_DOCS_DIR / "dataset-build-workflow.md"
DIFFERENTIAL_DOC = WORKFLOW_DOCS_DIR / "differential-analysis.md"
ENRICHMENT_DOC = WORKFLOW_DOCS_DIR / "enrichment.md"
KINASE_DOC = WORKFLOW_DOCS_DIR / "kinase.md"
SIGNALOME_DOC = WORKFLOW_DOCS_DIR / "signalome.md"


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
                "ATMSGRPRTTSFAESSKPVQQPSAFGQAAAL",
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
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                localisation=DatasetLocalisationConfig(
                    mode="require_threshold",
                    confidence_column="localisation_confidence",
                    min_confidence=0.75,
                ),
            ),
        )
    )


def _minimal_kinase_result() -> KinaseWorkflowResult:
    dataset = _build_dataset()
    site_keys = dataset.phospho.index.astype(str).tolist()
    display_ids = dataset.site_metadata.loc[:, "display_id"].astype(str).tolist()
    index = pd.Index(site_keys, name="site_key")
    kinases = pd.Index(["K1"], name="kinase")
    score_matrix = pd.DataFrame([[0.9], [0.8]], index=index, columns=kinases)
    prediction_matrix = pd.DataFrame([[0.7], [0.6]], index=index, columns=kinases)
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K1"],
                "substrate_site": display_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31, "C" * 31]},
            index=pd.Index(display_ids, name="site_id"),
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
    exec(
        """from phospy.api import (
    EnrichmentConfig,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
)""",
        namespace,
    )

    assert "AnalysisReadyDatasetBuilder" in namespace
    assert "DifferentialAnalysisWorkflow" in namespace
    assert "EnrichmentWorkflow" in namespace
    assert "DatasetBuildRequest" in namespace
    assert "DifferentialAnalysisRequest" in namespace
    assert "EnrichmentWorkflowRequest" in namespace
    assert "SignalomeWorkflowRequest" in namespace


def test_each_public_workflow_has_dedicated_api_page_with_contract_classes() -> None:
    expected = {
        DIFFERENTIAL_DOC: (
            "DifferentialAnalysisWorkflow",
            "DifferentialAnalysisRequest",
            "DifferentialAnalysisConfig",
            "DifferentialAnalysisResult",
        ),
        ENRICHMENT_DOC: (
            "EnrichmentWorkflow",
            "EnrichmentWorkflowRequest",
            "EnrichmentConfig",
            "EnrichmentWorkflowResult",
        ),
        KINASE_DOC: (
            "KinaseWorkflow",
            "KinaseWorkflowRequest",
            "KinaseScoringConfig",
            "KinasePredictionConfig",
            "KinaseActivityConfig",
            "KinaseWorkflowResult",
        ),
        SIGNALOME_DOC: (
            "SignalomeWorkflow",
            "SignalomeWorkflowRequest",
            "SignalomeConfig",
            "SignalomeWorkflowResult",
        ),
    }

    for path, class_names in expected.items():
        assert path.exists(), f"missing workflow API page: {path}"
        source = _read(path)
        for class_name in class_names:
            assert class_name in source, f"{class_name} missing from {path}"


def test_api_docs_enrichment_example_uses_public_api_and_runs_offline() -> None:
    source = _read(ENRICHMENT_DOC)
    normalized = " ".join(source.split())

    assert "## Minimal example" in source
    assert "EnrichmentWorkflow" in source
    assert "GeneSetCollection(" in source
    assert "background_universe=" in source
    assert "offline/no-online-resource policy" in source
    assert "GO, KEGG, Reactome" in source
    assert "Enrichr, gseapy, clusterProfiler" in normalized
    assert "does not implement GSEA, ssGSEA, or PTM-SEA" in normalized

    collection = GeneSetCollection(
        sets={
            "kinase_response": ("AKT1", "MAPK1", "MTOR"),
            "cell_cycle": ("CDK1", "CDK2", "MAPK1"),
        },
        identifier_kind="gene_symbol",
        term_names={
            "kinase_response": "Kinase response",
            "cell_cycle": "Cell cycle",
        },
        source_name="example in-memory gene sets",
        source_version="2026-06",
    )
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind="gene_symbol",
        set_collection=collection,
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(
            method="over_representation",
            multiple_testing_correction="benjamini_hochberg",
        ),
    )

    result = EnrichmentWorkflow().run(request)

    assert tuple(result.table["term_id"]) == ("kinase_response", "cell_cycle")
    assert result.provenance is not None
    assert result.provenance.workflow_parameters["background_universe_size"] == 5


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


def test_api_docs_batch_correction_example_is_constructible() -> None:
    source = _read(DATASET_BUILD_DOC)

    assert "DatasetBatchCorrectionConfig(" in source
    assert 'method="linear_residualize_batch"' in source
    assert "dataset.preprocessing_report.batch_correction" in source
    assert "confounding-check status" in source
    assert "not ComBat, not RUV" in source
    assert "not mixed-effects modelling" in source

    phospho = pd.DataFrame(
        {
            "sample_1": [10.0, 2.0],
            "sample_2": [15.0, 7.0],
            "sample_3": [14.0, 1.0],
            "sample_4": [19.0, 6.0],
        },
        index=["MAPK14;Y182;", "AKT1;T308;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            ],
            "display_id": ["MAPK14;Y182;", "AKT1;T308;"],
            "organism": ["rat", "rat"],
            "protein_namespace": ["protein_id", "protein_id"],
            "protein_identifier": ["MAPK14", "AKT1"],
            "protein_id": ["MAPK14", "AKT1"],
            "localisation_confidence": [0.95, 0.92],
        },
        index=phospho.index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {
            "batch": ["run_1", "run_2", "run_1", "run_2"],
            "condition": ["control", "control", "treated", "treated"],
        },
        index=phospho.columns.copy(),
    )
    preprocessing = DatasetPreprocessingConfig(
        batch_correction=DatasetBatchCorrectionConfig(
            method="linear_residualize_batch",
            batch_column="batch",
            condition_column="condition",
            preserve_condition_effects=True,
        )
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            organism=Organism.RAT,
            input_intensity_scale=IntensityScaleKind.LOG2,
            preprocessing_config=preprocessing,
        )
    )

    assert dataset.preprocessing_report is not None
    report = dataset.preprocessing_report.batch_correction
    assert report is not None
    assert report.status == "applied"
    assert report.method == "linear_residualize_batch"
    assert report.confounding_check_status == "passed"
    assert report.batch_levels == ("run_1", "run_2")
    assert report.condition_levels == ("control", "treated")


def test_api_docs_protein_aware_preparation_boundary_is_documented() -> None:
    source = _read(DATASET_BUILD_DOC)
    normalized = " ".join(source.split())

    assert "## Protein-Aware Preparation Parameters" in source
    assert 'policy="prepare_model_inputs"' in source
    assert "does not change the phosphosite matrix" in normalized
    assert "does not subtract total protein" in normalized
    assert "does not normalise intensities" in normalized
    assert "does not run differential analysis" in normalized
    assert "PhosPy does not claim MSstatsPTM-style inference" in source
    assert "dataset.protein_aware_preparation" in source
    assert "report.site_eligibility_dataframe()" in source


def test_api_docs_differential_request_example_is_constructible() -> None:
    dataset = _build_dataset()
    assert dataset.intensity_scale_state.kind.value == "log2"
    assert dataset.intensity_scale_state.is_established

    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="control",
                biological_replicate_id="control_r1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="control",
                biological_replicate_id="control_r2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="treatment",
                biological_replicate_id="treatment_r1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="treatment",
                biological_replicate_id="treatment_r2",
            ),
        )
    )
    contrasts = (
        Contrast(
            name="treatment_vs_control",
            numerator_condition="treatment",
            denominator_condition="control",
        ),
    )
    request = DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=contrasts,
    )

    assert request.design.samples[0].sample_id == "A_1"
    assert request.contrasts[0].name == "treatment_vs_control"
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


def test_readme_and_differential_docs_keep_scientific_scope_contracts() -> None:
    readme_source = _read(README)
    readme_normalized = " ".join(readme_source.split())
    differential_source = _read(DIFFERENTIAL_DOC)

    assert "minimum_condition_replicates=1" not in readme_source
    assert "minimum_condition_replicates=1" not in differential_source
    assert "KinaseWorkflow().run(" in readme_source
    assert "KinaseWorkflowRequest(" in readme_source
    assert "site_sequence" in readme_source
    assert "ReferencePreset.AUTO" in readme_source
    assert "`linear_residualize_batch`" in readme_source
    assert "rejects confounded batch/condition metadata" in readme_source
    assert "not ComBat, RUV, limma `removeBatchEffect` parity" in readme_source
    assert "mixed-effects modelling" in readme_normalized
    assert "control_rep1" in differential_source
    assert "control_rep2" in differential_source
    assert "treatment_rep1" in differential_source
    assert "treatment_rep2" in differential_source
    assert 'policy="log2"' in differential_source


def test_public_workflow_docs_make_localisation_policy_explicit() -> None:
    readme_source = _read(README)
    dataset_source = _read(DATASET_BUILD_DOC)
    differential_source = _read(DIFFERENTIAL_DOC)
    kinase_source = _read(KINASE_DOC)
    signalome_source = _read(SIGNALOME_DOC)

    assert "DatasetLocalisationConfig(" in readme_source
    assert 'confidence_column="localisation_confidence"' in readme_source
    assert "min_confidence=0.75" in readme_source

    assert "## Localisation-Confidence Parameters" in dataset_source
    assert "fails dataset build when" in dataset_source
    assert 'mode="require_threshold"' in dataset_source
    assert 'confidence_column="localisation_confidence"' in dataset_source
    assert "min_confidence=0.75" in dataset_source

    assert "DatasetLocalisationConfig(" in differential_source
    assert "localisation_confidence" in differential_source
    assert "low-confidence phosphosite assignments fail fast" in differential_source

    assert "## Localisation Prerequisite" in kinase_source
    assert "DatasetLocalisationConfig(" in kinase_source
    assert "dataset build fails" in kinase_source

    assert "## Localisation Prerequisite" in signalome_source
    assert "DatasetLocalisationConfig(" in signalome_source
    assert "dataset build fails" in signalome_source


def test_kinase_docs_explain_reference_display_ambiguity_policy() -> None:
    source = _read(KINASE_DOC)

    assert "reference_display_ambiguity_policy" in source
    assert '"error"' in source
    assert '"allow_with_diagnostics"' in source
    assert "matched `site_key` values" in source
    assert "does not collapse duplicate display labels" in source

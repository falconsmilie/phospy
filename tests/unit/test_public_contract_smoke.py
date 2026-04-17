from __future__ import annotations

import pandas as pd

import phospy
import phospy_legacy
from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DatasetBuildRequest,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflow,
    SignalomeWorkflowRequest,
    SimpleKinaseWorkflow,
    SimpleKinaseWorkflowRequest,
)


def _dataset() -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=["MAPK14;Y182;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
        },
        index=phospho.index,
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
    )


def test_package_imports_from_new_tree() -> None:
    assert "phospy\\__init__.py" in phospy.__file__.replace("/", "\\")
    assert "legacy" in (phospy_legacy.__doc__ or "").lower()


def test_public_shells_import_and_instantiate() -> None:
    dataset = _dataset()
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["MAP2K6"], "substrate_site": ["MAPK14;Y182;"]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
    )
    request = SimpleKinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(),
        prediction_config=KinasePredictionConfig(),
        activity_config=KinaseActivityConfig(),
    )
    assert isinstance(request, SimpleKinaseWorkflowRequest)
    assert isinstance(SignalomeConfig(), SignalomeConfig)
    assert ReferencePreset.AUTO.value == "auto"


def test_builder_and_workflows_expose_run_and_execute_shells() -> None:
    dataset = _dataset()

    builder = AnalysisReadyDatasetBuilder()
    assert callable(getattr(builder, "run", None))
    built = builder.run(
        DatasetBuildRequest(
            phospho=dataset.phospho,
            site_metadata=dataset.site_metadata,
            organism=dataset.organism,
        )
    )
    assert isinstance(built, AnalysisReadyPhosphoDataset)

    kinase_workflow = SimpleKinaseWorkflow()
    assert callable(getattr(kinase_workflow, "run", None))
    kinase_result = kinase_workflow.run(
        SimpleKinaseWorkflowRequest(
            dataset=built,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=1),
            prediction_config=KinasePredictionConfig(top_k=1, ensemble_size=3),
            activity_config=KinaseActivityConfig(enabled=False),
        )
    )
    assert isinstance(kinase_result.prediction_result.pred_mat, pd.DataFrame)
    assert not hasattr(kinase_result, "pred_mat")

    signalome_workflow = SignalomeWorkflow()
    assert callable(getattr(signalome_workflow, "run", None))
    signalome_request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=SignalomeConfig(),
    )
    assert isinstance(signalome_request, SignalomeWorkflowRequest)

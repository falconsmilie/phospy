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
    phospho = pd.DataFrame({"sample_a": [1.0]}, index=["GENEA;S1;"])
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["GENEA"],
            "site": ["S1"],
            "site_sequence": ["AAAAAAA"],
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
            {"kinase": ["KINASE_A"], "substrate_site": ["GENEA;S1;"]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["AAAAAAA"]},
            index=pd.Index(["GENEA;S1;"], name="site_id"),
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
        SimpleKinaseWorkflowRequest(dataset=built, references=ReferencePreset.AUTO)
    )
    assert isinstance(kinase_result.prediction_result.pred_mat, pd.DataFrame)
    assert not hasattr(kinase_result, "pred_mat")

    signalome_workflow = SignalomeWorkflow()
    assert callable(getattr(signalome_workflow, "run", None))
    signalome_result = signalome_workflow.run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(),
        )
    )
    assert signalome_result.kinase_result is kinase_result

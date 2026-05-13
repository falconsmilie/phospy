from __future__ import annotations

import inspect

import phospy.workflows.kinase.interpreter as kinase_interpreter_module
import phospy.workflows.signalome.interpreter as signalome_interpreter_module
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.references.models import ReferenceBundle


def test_dataset_boundary_constructor_has_no_canonicalization_or_default_repair() -> (
    None
):
    post_init_source = inspect.getsource(AnalysisReadyPhosphoDataset.__post_init__)
    from_owned_source = inspect.getsource(AnalysisReadyPhosphoDataset._from_owned)
    assert "canonicalize_site_" not in post_init_source
    assert "IntensityScaleState.raw" not in from_owned_source


def test_reference_boundary_constructor_has_no_canonicalization_or_dedup_repair() -> (
    None
):
    post_init_source = inspect.getsource(ReferenceBundle.__post_init__)
    assert "canonicalize_site_" not in post_init_source
    assert "drop_duplicates" not in post_init_source
    assert ".str.strip" not in post_init_source


def test_workflow_interpreters_have_no_site_id_repair_calls() -> None:
    kinase_source = inspect.getsource(kinase_interpreter_module)
    signalome_source = inspect.getsource(signalome_interpreter_module)
    assert "canonicalize_site_" not in kinase_source
    assert "canonicalize_site_" not in signalome_source

from __future__ import annotations

import inspect

from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.references.models import ReferenceBundle


def test_dataset_boundary_constructor_has_no_canonicalization_or_default_repair() -> (
    None
):
    post_init_source = inspect.getsource(AnalysisReadyPhosphoDataset.__post_init__)
    from_owned_source = inspect.getsource(AnalysisReadyPhosphoDataset._from_owned)
    assert "canonicalize_site_" not in post_init_source
    assert "TransformationState.raw" not in from_owned_source


def test_reference_boundary_constructor_has_no_canonicalization_or_dedup_repair() -> (
    None
):
    post_init_source = inspect.getsource(ReferenceBundle.__post_init__)
    assert "canonicalize_site_" not in post_init_source
    assert "drop_duplicates" not in post_init_source
    assert ".str.strip" not in post_init_source

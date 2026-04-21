from __future__ import annotations

from pathlib import Path

import phospy
import phospy.api as public_api

ROOT = Path(__file__).resolve().parents[2]

USER_FACING_IMPORT_FILES = (
    ROOT / "README.md",
    ROOT / "docs" / "api.md",
    ROOT / "examples" / "dataset_builder_demo.py",
    ROOT / "examples" / "kinase_workflow_demo.py",
    ROOT / "examples" / "signalome_workflow_demo.py",
)


def test_top_level_facade_re_exports_canonical_api_namespace() -> None:
    assert set(public_api.__all__).issubset(set(phospy.__all__))
    for exported in public_api.__all__:
        assert getattr(phospy, exported) is getattr(public_api, exported)


def test_readme_and_api_guide_document_import_contract() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    api_guide = (ROOT / "docs" / "api.md").read_text(encoding="utf-8")

    canonical_phrase = (
        "`phospy.api` is the canonical namespace where public API types are defined"
    )
    primary_route_phrase = "Top-level `phospy` is the primary supported import route"

    assert canonical_phrase in readme
    assert canonical_phrase in api_guide
    assert primary_route_phrase in readme
    assert primary_route_phrase in api_guide


def test_user_facing_guides_and_examples_use_top_level_import_route() -> None:
    for file_path in USER_FACING_IMPORT_FILES:
        source = file_path.read_text(encoding="utf-8")
        assert "from phospy.api import" not in source
        assert "import phospy.api" not in source

from __future__ import annotations

import importlib
import tomllib
from pathlib import Path

import pytest

from phospy.api._compat import compatibility_exports

ROOT = Path(__file__).resolve().parents[3]


@pytest.mark.parametrize(
    ("statement", "name", "replacement_module"),
    (
        (
            "from phospy.api import KinaseScoringConfig",
            "KinaseScoringConfig",
            "phospy.advanced",
        ),
        (
            "from phospy.api.configs import KinaseScoringConfig",
            "KinaseScoringConfig",
            "phospy.advanced.configs",
        ),
        (
            "from phospy.api.configs.kinase import KinaseScoringConfig",
            "KinaseScoringConfig",
            "phospy.advanced.configs",
        ),
        (
            "from phospy.api.requests import SignalomeConfig",
            "SignalomeConfig",
            "phospy.advanced",
        ),
        (
            "from phospy.api.results import KinaseEligibilityReport",
            "KinaseEligibilityReport",
            "phospy.advanced.results",
        ),
        (
            "from phospy.api.configs import ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL",
            "ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL",
            "phospy.contracts.configs",
        ),
        (
            "from phospy.api.configs.kinase import KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED",
            "KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED",
            "phospy.contracts.configs.kinase",
        ),
    ),
)
def test_compatibility_imports_warn_with_actionable_replacement_path(
    statement: str,
    name: str,
    replacement_module: str,
) -> None:
    namespace: dict[str, object] = {}
    replacement = importlib.import_module(replacement_module)

    with pytest.warns(DeprecationWarning) as records:
        exec(statement, namespace)

    assert namespace[name] is getattr(replacement, name)
    warning_text = str(records[0].message)
    assert f"use `from {replacement_module} import {name}`" in warning_text
    assert "introduced in PhosPy 1.6.0" in warning_text
    assert "planned for removal in PhosPy 2.0.0" in warning_text


def test_every_compatibility_export_has_policy_metadata_and_live_owner() -> None:
    exports = compatibility_exports()
    keys = [(export.old_module, export.name) for export in exports]
    current_version = _version_tuple(
        tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"][
            "version"
        ]
    )

    assert exports
    assert len(keys) == len(set(keys))

    for export in exports:
        assert export.old_module.startswith("phospy.api")
        assert export.owner_module
        assert export.replacement_module
        assert export.introduced_version == "1.6.0"
        assert export.planned_removal_version == "2.0.0"
        assert _version_tuple(export.planned_removal_version) > current_version
        assert export.stability in {"advanced", "unsupported"}

        owner = importlib.import_module(export.owner_module)
        assert export.name in getattr(owner, "__all__", ()) or hasattr(
            owner,
            export.name,
        )


def test_advanced_compatibility_exports_point_to_advanced_namespace() -> None:
    for export in compatibility_exports():
        if export.stability != "advanced":
            continue
        assert export.replacement_module in {
            "phospy.advanced",
            "phospy.advanced.configs",
            "phospy.advanced.results",
        }


def test_unregistered_cross_submodule_compatibility_route_fails_closed() -> None:
    with pytest.raises(ImportError):
        exec("from phospy.api.configs.localisation import KinaseScoringConfig", {})


def _version_tuple(version: str) -> tuple[int, int, int]:
    major, minor, patch = version.split(".")
    return int(major), int(minor), int(patch)

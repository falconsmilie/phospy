from __future__ import annotations

import importlib
import tomllib
from pathlib import Path

import pytest

from phospy._deprecations import PhosPyDeprecationWarning, retained_deprecations
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

    with pytest.warns(PhosPyDeprecationWarning) as records:
        exec(statement, namespace)

    assert namespace[name] is getattr(replacement, name)
    warning_text = str(records[0].message)
    assert f"use `from {replacement_module} import {name}`" in warning_text
    assert "introduced in PhosPy 1.6.0" in warning_text
    assert "planned for removal in PhosPy 2.0.0" in warning_text


def test_every_compatibility_export_has_policy_metadata_and_live_owner() -> None:
    exports = compatibility_exports()
    keys = [(export.old_module, export.name) for export in exports]
    current_version = _current_project_version()

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


def test_every_retained_deprecation_has_complete_unique_metadata_and_live_replacement() -> (
    None
):
    records = retained_deprecations()
    identifiers = [record.identifier for record in records]
    deprecated_targets = [record.deprecated for record in records]
    current_version = _current_project_version()

    assert records
    assert len(identifiers) == len(set(identifiers))
    assert len(deprecated_targets) == len(set(deprecated_targets))

    for record in records:
        assert record.identifier
        assert record.kind
        assert record.owner_module
        assert record.deprecated
        assert record.replacement
        assert record.introduced_version
        assert record.planned_removal_version
        assert _version_tuple(record.planned_removal_version) > current_version
        assert record.stability in {"stable", "advanced", "unsupported", "internal"}
        assert record.replacement_module
        assert record.replacement_name
        assert record.source_uses

        importlib.import_module(record.owner_module)
        replacement_owner = importlib.import_module(record.replacement_module)
        assert record.replacement_name in getattr(
            replacement_owner,
            "__all__",
            (),
        ) or hasattr(replacement_owner, record.replacement_name)

        for source_use in record.source_uses:
            assert source_use.kind
            assert source_use.token
            if source_use.kind == "import-route":
                assert source_use.module
            if source_use.check_consumer_sources:
                assert not source_use.unchecked_reason
            else:
                assert source_use.unchecked_reason


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


def _current_project_version() -> tuple[int, int, int]:
    config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    return _version_tuple(config["project"]["version"])

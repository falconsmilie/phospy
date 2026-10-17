from __future__ import annotations

import importlib
import json
from pathlib import Path

from documentation_import_policy import (
    documented_import_problems,
    user_facing_markdown_paths,
)

import phospy
import phospy.advanced
import phospy.api

ROOT = Path(__file__).resolve().parents[2]
PUBLIC_BASELINE = ROOT / "architecture/public-api-baseline.json"


def _inventories() -> dict[str, frozenset[str]]:
    return {
        "phospy": frozenset(phospy.__all__),
        "phospy.api": frozenset(phospy.api.__all__),
        "phospy.advanced": frozenset(phospy.advanced.__all__),
    }


def test_every_current_user_facing_markdown_import_is_supported() -> None:
    paths = user_facing_markdown_paths(ROOT)
    assert paths
    assert all(path.is_file() for path in paths)
    problems: list[str] = []
    documented = []
    for path in paths:
        path_problems, path_imports = documented_import_problems(
            path=path.relative_to(ROOT).as_posix(),
            markdown=path.read_text(encoding="utf-8"),
            inventories=_inventories(),
        )
        problems.extend(path_problems)
        documented.extend(path_imports)

    assert problems == []
    assert documented


def test_documented_aggregate_names_resolve_to_exact_owner_identities() -> None:
    baseline = json.loads(PUBLIC_BASELINE.read_text(encoding="utf-8"))["surfaces"]
    definitions = {
        module: {item["name"]: item for item in surface["definitions"]}
        for module, surface in baseline.items()
    }
    documented = []
    for path in user_facing_markdown_paths(ROOT):
        problems, path_imports = documented_import_problems(
            path=path.relative_to(ROOT).as_posix(),
            markdown=path.read_text(encoding="utf-8"),
            inventories=_inventories(),
        )
        assert problems == []
        documented.extend(path_imports)

    for item in documented:
        facade = importlib.import_module(item.module)
        for name in item.names:
            record = definitions[item.module][name]
            value = getattr(facade, name)
            if record["kind"] == "callable_or_type":
                module_name, qualname = record["definition"].split(":", maxsplit=1)
                owner_value = importlib.import_module(module_name)
                for part in qualname.split("."):
                    owner_value = getattr(owner_value, part)
                assert value is owner_value
            else:
                assert repr(value) == record["value_repr"]


def test_documentation_import_policy_rejects_bypass_variants() -> None:
    variants = {
        "owner": "from phospy.datasets.models import AnalysisReadyPhosphoDataset",
        "retired": "from phospy.tables import DatasetTable",
        "wildcard": "from phospy.api import *",
        "alias": "from phospy.api import Organism as Species",
        "dynamic": 'import importlib\nimportlib.import_module("phospy.api")',
        "malformed": "from phospy.api import (Organism",
    }
    for name, source in variants.items():
        problems, _imports = documented_import_problems(
            path=f"{name}.md",
            markdown=f"```python\n{source}\n```\n",
            inventories=_inventories(),
        )
        assert problems, name

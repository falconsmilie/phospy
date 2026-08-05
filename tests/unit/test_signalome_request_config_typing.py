from __future__ import annotations

import ast
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import get_type_hints

import pytest

import phospy.advanced as advanced_api
import phospy.advanced.configs as advanced_configs
import phospy.api.requests as request_api
from phospy.api.requests import SignalomeWorkflowRequest
from phospy.contracts.configs import SignalomeConfig as ContractSignalomeConfig
from phospy.contracts.results.kinase import KinaseWorkflowResult

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PHOSPY_ROOT = SRC_ROOT / "phospy"
API_ROOT = PHOSPY_ROOT / "api"


def test_signalome_request_annotations_are_defined_without_runtime_patching() -> None:
    assert SignalomeWorkflowRequest.__annotations__ == {
        "kinase_result": "KinaseWorkflowResult",
        "config": "SignalomeConfig",
    }

    hints = get_type_hints(SignalomeWorkflowRequest)
    init_hints = get_type_hints(SignalomeWorkflowRequest.__init__)

    assert hints["kinase_result"] is KinaseWorkflowResult
    assert hints["config"] is ContractSignalomeConfig
    assert init_hints["kinase_result"] is KinaseWorkflowResult
    assert init_hints["config"] is ContractSignalomeConfig


def test_signalome_request_documentation_signature_matches_type_hints() -> None:
    import inspect

    signature = inspect.signature(SignalomeWorkflowRequest)
    config_parameter = signature.parameters["config"]

    assert config_parameter.annotation == "SignalomeConfig"
    assert get_type_hints(SignalomeWorkflowRequest)["config"] is ContractSignalomeConfig


def test_signalome_request_type_hints_do_not_depend_on_advanced_import_order() -> None:
    code = "\n".join(
        [
            "import sys",
            "from typing import get_type_hints",
            "from phospy.api.requests import SignalomeWorkflowRequest",
            "from phospy.contracts.configs import SignalomeConfig",
            "assert 'phospy.advanced.configs' not in sys.modules",
            "assert SignalomeWorkflowRequest.__annotations__['config'] == 'SignalomeConfig'",
            "assert get_type_hints(SignalomeWorkflowRequest)['config'] is SignalomeConfig",
            "assert 'phospy.advanced.configs' not in sys.modules",
        ]
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        [
            str(SRC_ROOT),
            environment.get("PYTHONPATH", ""),
        ]
    )

    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr


def test_signalome_config_reexports_preserve_one_class_identity() -> None:
    assert advanced_api.SignalomeConfig is ContractSignalomeConfig
    assert advanced_configs.SignalomeConfig is ContractSignalomeConfig

    with pytest.warns(
        DeprecationWarning,
        match=r"from phospy\.advanced import SignalomeConfig",
    ):
        from phospy.api.requests import SignalomeConfig as RequestCompatConfig

    with pytest.warns(
        DeprecationWarning,
        match=r"from phospy\.advanced\.configs import SignalomeConfig",
    ):
        from phospy.api.configs import SignalomeConfig as CompatConfig

    assert RequestCompatConfig is ContractSignalomeConfig
    assert CompatConfig is ContractSignalomeConfig


def test_stable_request_api_does_not_own_or_import_signalome_config() -> None:
    assert "SignalomeConfig" not in request_api.__dict__
    assert "SignalomeConfig" not in request_api.__all__


def test_stable_api_modules_do_not_static_import_advanced_modules() -> None:
    offenders: list[str] = []
    for path in sorted(API_ROOT.rglob("*.py")):
        if path.name == "_compat.py":
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for imported_module, line_number in _static_imports(tree):
            if imported_module == "phospy.advanced" or imported_module.startswith(
                "phospy.advanced."
            ):
                relative = path.relative_to(PROJECT_ROOT)
                offenders.append(f"{relative}:{line_number} -> {imported_module}")

    assert offenders == []


def test_source_does_not_mutate_runtime_annotations() -> None:
    offenders: list[str] = []
    for path in sorted(PHOSPY_ROOT.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if _mutates_annotations(node):
                relative = path.relative_to(PROJECT_ROOT)
                offenders.append(f"{relative}:{node.lineno}")

    assert offenders == []


def test_pyright_accepts_signalome_request_construction_with_advanced_reexport(
    tmp_path: Path,
) -> None:
    pytest.importorskip("pyright")

    typing_smoke = tmp_path / "signalome_request_construction.py"
    typing_smoke.write_text(
        "\n".join(
            [
                "from __future__ import annotations",
                "",
                "from typing import assert_type",
                "",
                "from phospy.advanced.configs import SignalomeConfig",
                "from phospy.api.requests import SignalomeWorkflowRequest",
                "from phospy.contracts.results.kinase import KinaseWorkflowResult",
                "",
                "def build(result: KinaseWorkflowResult) -> SignalomeWorkflowRequest:",
                "    explicit = SignalomeWorkflowRequest(",
                "        kinase_result=result,",
                "        config=SignalomeConfig.production(),",
                "    )",
                "    assert_type(explicit.config, SignalomeConfig)",
                "    defaulted = SignalomeWorkflowRequest(kinase_result=result)",
                "    assert_type(defaulted.config, SignalomeConfig)",
                "    return explicit",
                "",
            ]
        ),
        encoding="utf-8",
    )
    (tmp_path / "pyrightconfig.json").write_text(
        json.dumps(
            {
                "include": [typing_smoke.name],
                "extraPaths": [str(SRC_ROOT)],
                "pythonVersion": "3.11",
                "typeCheckingMode": "standard",
                "reportMissingTypeStubs": "none",
            }
        ),
        encoding="utf-8",
    )

    completed = subprocess.run(
        [sys.executable, "-m", "pyright", "--project", str(tmp_path)],
        cwd=PROJECT_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr


def _static_imports(tree: ast.AST) -> tuple[tuple[str, int], ...]:
    imports: list[tuple[str, int]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend((alias.name, node.lineno) for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imports.append((node.module, node.lineno))
    return tuple(imports)


def _mutates_annotations(node: ast.AST) -> bool:
    if isinstance(node, ast.Assign | ast.AnnAssign | ast.AugAssign):
        return any(_targets_annotations(target) for target in _assignment_targets(node))
    return False


def _assignment_targets(
    node: ast.Assign | ast.AnnAssign | ast.AugAssign,
) -> tuple[ast.AST, ...]:
    if isinstance(node, ast.Assign):
        return tuple(node.targets)
    return (node.target,)


def _targets_annotations(target: ast.AST) -> bool:
    if isinstance(target, ast.Attribute):
        return target.attr == "__annotations__"
    if isinstance(target, ast.Subscript):
        return _targets_annotations(target.value)
    if isinstance(target, ast.Tuple | ast.List):
        return any(_targets_annotations(element) for element in target.elts)
    return False

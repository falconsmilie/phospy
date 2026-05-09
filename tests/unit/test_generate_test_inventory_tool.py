from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOL_PATH = REPO_ROOT / "tools" / "testing" / "generate_test_inventory.py"
MANUAL_COLUMNS = (
    "primary_category",
    "recommended_action",
    "protected_risk_or_contract",
    "reviewer",
    "followup_ticket",
    "notes",
)


def _load_tool_module():
    module_name = "test_generate_test_inventory_tool_runtime"
    spec = importlib.util.spec_from_file_location(module_name, TOOL_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(module_name, None)
    return module


def _write_existing_inventory(
    csv_path: Path,
    *,
    encoding: str,
    path_value: str,
) -> dict[str, str]:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    manual_values = {
        "primary_category": "protected",
        "recommended_action": "keep",
        "protected_risk_or_contract": "parity infrastructure",
        "reviewer": "codex",
        "followup_ticket": "TST-AUDIT-TOOL-001",
        "notes": "retain manual values",
    }
    with csv_path.open("w", encoding=encoding, newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "path",
                "test_area",
                "line_count",
                "test_function_count",
                "uses_pytest_mark",
                "uses_parametrization",
                *MANUAL_COLUMNS,
            ]
        )
        writer.writerow(
            [
                path_value,
                "unit",
                1,
                1,
                "false",
                "false",
                *(manual_values[column] for column in MANUAL_COLUMNS),
            ]
        )
    return manual_values


def _read_inventory_row(csv_path: Path, path_value: str) -> dict[str, str]:
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames is not None
        for row in reader:
            if row.get("path") == path_value:
                return row
    raise AssertionError(f"missing row for path '{path_value}'")


@pytest.mark.parametrize("existing_encoding", ("utf-8", "utf-8-sig"))
def test_inventory_generator_preserves_manual_columns_with_or_without_bom(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, existing_encoding: str
) -> None:
    module = _load_tool_module()
    repo_root = tmp_path
    tests_root = repo_root / "tests" / "unit"
    docs_root = repo_root / "docs" / "testing"
    tests_root.mkdir(parents=True, exist_ok=True)

    test_file = tests_root / "test_example_inventory.py"
    test_source = (
        "import pytest\n\n"
        "@pytest.mark.unit\n"
        "def test_one() -> None:\n"
        "    assert True\n\n"
        "def test_two() -> None:\n"
        "    assert True\n"
    )
    test_file.write_text(test_source, encoding="utf-8")

    csv_path = docs_root / "test_inventory.csv"
    expected_path = "tests/unit/test_example_inventory.py"
    manual_values = _write_existing_inventory(
        csv_path,
        encoding=existing_encoding,
        path_value=expected_path,
    )

    monkeypatch.setattr(module, "REPO_ROOT", repo_root)
    monkeypatch.setattr(module, "TESTS_ROOT", repo_root / "tests")
    monkeypatch.setattr(module, "OUTPUT_DIR", docs_root)
    monkeypatch.setattr(module, "CSV_OUTPUT", csv_path)
    monkeypatch.setattr(module, "MARKDOWN_OUTPUT", docs_root / "test_inventory.md")

    module.main()

    row = _read_inventory_row(csv_path, expected_path)
    for column, expected_value in manual_values.items():
        assert row[column] == expected_value

    expected_line_count = len(test_source.splitlines())
    assert row["line_count"] == str(expected_line_count)
    assert row["test_function_count"] == "2"
    assert row["uses_pytest_mark"] == "true"
    assert row["uses_parametrization"] == "false"

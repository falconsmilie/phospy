from __future__ import annotations

from pathlib import Path

import pytest

from phospy.errors import RequestValidationError
from phospy.validation.requests import validate_existing_file_path


def test_validate_existing_file_path_returns_existing_file(tmp_path: Path) -> None:
    path = tmp_path / "input.tsv"
    path.write_text("a\tb\n1\t2\n")

    validated = validate_existing_file_path(path, context="input table path")

    assert validated == path


def test_validate_existing_file_path_rejects_missing_path(tmp_path: Path) -> None:
    path = tmp_path / "missing.tsv"

    with pytest.raises(RequestValidationError, match="Path does not exist"):
        validate_existing_file_path(path, context="input table path")


def test_validate_existing_file_path_rejects_directory(tmp_path: Path) -> None:
    path = tmp_path / "input_dir"
    path.mkdir()

    with pytest.raises(RequestValidationError, match="Path is not a file"):
        validate_existing_file_path(path, context="input table path")

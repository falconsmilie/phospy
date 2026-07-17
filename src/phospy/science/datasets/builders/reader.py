"""Internal source-reader collaborator for dataset builder inputs."""

from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import Protocol, cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.builders.contracts import DatasetInput


class DatasetInputSourceValidatorProtocol(Protocol):
    """Validate one dataset builder input source before reading."""

    def run(
        self,
        value: object | None,
        *,
        field_name: str,
        allow_none: bool = False,
    ) -> object: ...


class DatasetInputSourceValidator:
    """Validate supported source types for dataset build inputs."""

    def run(
        self,
        source: object | None,
        *,
        field_name: str,
        allow_none: bool = False,
    ) -> pd.DataFrame | str | Path | None:
        if source is None:
            if allow_none:
                return None
            raise PhosPyInputError(f"dataset build request {field_name} is required")
        if isinstance(source, pd.DataFrame):
            return source
        if isinstance(source, (str, Path, PathLike)):
            if isinstance(source, str) and not source.strip():
                raise PhosPyInputError(
                    f"dataset build request {field_name} path cannot be empty"
                )
            if isinstance(source, (str, Path)):
                return source
            return Path(cast(PathLike[str], source))
        raise PhosPyInputError(
            f"dataset build request {field_name} must be a pandas DataFrame or a "
            "file path (str/pathlib.Path)"
        )


class DatasetPathTableReader(Protocol):
    """Reader protocol for dataset builder local table paths."""

    def run(self, path: Path, *, field_name: str) -> pd.DataFrame: ...


class DatasetInputReader:
    """Resolve one supported builder input source into a DataFrame."""

    def __init__(
        self,
        *,
        source_validator: DatasetInputSourceValidatorProtocol | None = None,
        path_reader: DatasetPathTableReader | None = None,
    ) -> None:
        self._source_validator = source_validator or DatasetInputSourceValidator()
        self._path_reader = path_reader

    def run(self, source: DatasetInput, *, field_name: str) -> pd.DataFrame:
        validated_source = self._source_validator.run(source, field_name=field_name)
        if isinstance(validated_source, pd.DataFrame):
            return validated_source
        if not isinstance(validated_source, (str, Path, PathLike)):
            raise PhosPyInputError(
                f"dataset build request {field_name} source validator produced "
                f"unsupported value type {type(validated_source).__name__}; expected "
                "a pandas DataFrame or file path"
            )
        return self._read_from_path(validated_source, field_name=field_name)

    def _read_from_path(
        self,
        source: str | Path | PathLike,
        *,
        field_name: str,
    ) -> pd.DataFrame:
        if self._path_reader is None:
            raise PhosPyInputError(
                "dataset build request path inputs require an injected "
                "DatasetPathTableReader"
            )
        path = Path(source.strip()) if isinstance(source, str) else Path(source)
        return self._path_reader.run(path, field_name=field_name)

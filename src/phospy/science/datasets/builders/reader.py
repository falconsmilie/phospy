"""Internal source-reader collaborator for dataset builder inputs."""

from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import Protocol

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.builders.contracts import DatasetInput
from phospy.validation.datasets.inputs import DatasetInputSourceValidator


class DatasetPathTableReader(Protocol):
    """Reader protocol for dataset builder local table paths."""

    def run(self, path: Path, *, field_name: str) -> pd.DataFrame: ...


class DatasetInputReader:
    """Resolve one supported builder input source into a DataFrame."""

    def __init__(
        self,
        *,
        source_validator: DatasetInputSourceValidator | None = None,
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

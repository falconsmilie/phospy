"""Internal source-reader collaborator for dataset builder inputs."""

from __future__ import annotations

from os import PathLike
from pathlib import Path

import pandas as pd

from phospy.datasets.builders.contracts import DatasetInput
from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError
from phospy.validation.datasets.inputs import DatasetInputSourceValidator


class DatasetInputReader:
    """Resolve one supported builder input source into a DataFrame."""

    def __init__(
        self, *, source_validator: DatasetInputSourceValidator | None = None
    ) -> None:
        self._source_validator = source_validator or DatasetInputSourceValidator()

    def run(self, source: DatasetInput, *, field_name: str) -> pd.DataFrame:
        validated_source = self._source_validator.run(source, field_name=field_name)
        if isinstance(validated_source, pd.DataFrame):
            return validated_source.copy(deep=True)
        assert isinstance(validated_source, (str, Path, PathLike))
        return self._read_from_path(validated_source, field_name=field_name)

    @staticmethod
    def _read_from_path(
        source: str | Path | PathLike, *, field_name: str
    ) -> pd.DataFrame:
        from phospy.io.readers.tables import read_table

        path = Path(source.strip()) if isinstance(source, str) else Path(source)
        try:
            return read_table(path)
        except UnsupportedInputFormatError as exc:
            raise UnsupportedInputFormatError(
                f"dataset build request {field_name} input format is unsupported at "
                f"'{path}'. supported formats: csv, tsv, parquet"
            ) from exc
        except PhosPyInputError as exc:
            raise PhosPyInputError(
                f"failed to read dataset build request {field_name} from '{path}': {exc}"
            ) from exc

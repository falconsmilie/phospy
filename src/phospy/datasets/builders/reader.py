"""Internal source-reader collaborator for dataset builder inputs."""

from __future__ import annotations

from os import PathLike
from pathlib import Path

import pandas as pd

from phospy.datasets.builders.contracts import DatasetInput
from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError


class DatasetInputReader:
    """Resolve one supported builder input source into a DataFrame."""

    def run(self, source: DatasetInput, *, field_name: str) -> pd.DataFrame:
        if isinstance(source, pd.DataFrame):
            return source.copy(deep=True)
        if isinstance(source, str | Path | PathLike):
            return self._read_from_path(source, field_name=field_name)
        raise UnsupportedInputFormatError(
            f"dataset build request {field_name} must be a pandas DataFrame or a file "
            "path (str/pathlib.Path)"
        )

    @staticmethod
    def _read_from_path(
        source: str | Path | PathLike, *, field_name: str
    ) -> pd.DataFrame:
        from phospy.io.tables import read_table

        path = Path(source)
        if isinstance(source, str):
            raw = source.strip()
            if not raw:
                raise UnsupportedInputFormatError(
                    f"dataset build request {field_name} path cannot be empty"
                )
            path = Path(raw)
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

"""Internal source-reader collaborator for dataset builder inputs."""

from __future__ import annotations

from os import PathLike
from pathlib import Path

import pandas as pd

from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError
from phospy.science.datasets.builders.contracts import DatasetInput
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
            return validated_source
        if not isinstance(validated_source, (str, Path, PathLike)):
            raise PhosPyInputError(
                f"dataset build request {field_name} source validator produced "
                f"unsupported value type {type(validated_source).__name__}; expected "
                "a pandas DataFrame or file path"
            )
        return self._read_from_path(validated_source, field_name=field_name)

    @staticmethod
    def _read_from_path(
        source: str | Path | PathLike, *, field_name: str
    ) -> pd.DataFrame:
        from phospy.io.readers.tables import (
            read_phospho_matrix,
            read_sample_metadata,
            read_site_metadata,
            read_table,
            read_total_matrix,
            supported_table_input_formats,
        )

        path = Path(source.strip()) if isinstance(source, str) else Path(source)
        reader_by_field = {
            "phospho": read_phospho_matrix,
            "site_metadata": read_site_metadata,
            "sample_metadata": read_sample_metadata,
            "total": read_total_matrix,
        }
        read_from_path = reader_by_field.get(field_name, read_table)
        try:
            return read_from_path(path)
        except UnsupportedInputFormatError as exc:
            raise UnsupportedInputFormatError(
                f"dataset build request {field_name} has unsupported file format at "
                f"'{path}'. supported formats: {supported_table_input_formats()}"
            ) from exc
        except PhosPyInputError as exc:
            raise PhosPyInputError(
                f"failed to read dataset build request {field_name} from '{path}': {exc}"
            ) from exc

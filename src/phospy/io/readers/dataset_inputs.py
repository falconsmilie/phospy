"""Concrete dataset-builder path reader adapters."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError
from phospy.io.readers.tables import (
    read_phospho_matrix,
    read_sample_metadata,
    read_site_metadata,
    read_table,
    read_total_matrix,
    supported_table_input_formats,
)


class DatasetPathTableReader:
    """Read dataset-builder local table paths with field-specific schemas."""

    def run(self, path: Path, *, field_name: str) -> pd.DataFrame:
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


__all__ = ["DatasetPathTableReader"]

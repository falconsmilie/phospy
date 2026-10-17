from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from phospy.advanced import (
    DatasetIntensityTransformConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.api import DatasetPreprocessingConfig
from phospy.datasets.construction.input_resolution import (
    DatasetInputReader,
    _DatasetInputSourceValidator,
)
from phospy.datasets.contracts import DatasetBuildRequest
from phospy.datasets.validation.inputs import DatasetInputSourceValidator
from phospy.datasets.validation.request import DatasetBuildRequestValidator
from phospy.errors import PhosPyInputError, UnsupportedInputFormatError
from phospy.io.importers.contracts import PhosphositeImportRequest
from phospy.io.importers.validation import (
    PhosphositeImportRequestValidator,
    _PhosphositeImportSourceValidator,
)


class _GenericPathLike(PathLike[str]):
    def __fspath__(self) -> str:
        return "generic.csv"


class _InjectedFailure(RuntimeError):
    pass


class _RaisingSourceValidator:
    def run(
        self,
        value: object | None,
        *,
        field_name: str,
        allow_none: bool = False,
    ) -> object:
        raise _InjectedFailure(f"source:{field_name}")


class _RaisingPathReader:
    def run(self, path: Path, *, field_name: str) -> pd.DataFrame:
        raise _InjectedFailure(f"path:{field_name}:{path}")


class _FalseValuedSourceValidator:
    def __init__(self) -> None:
        self.calls: list[tuple[object | None, str, bool]] = []

    def __bool__(self) -> bool:
        return False

    def run(
        self,
        value: object | None,
        *,
        field_name: str,
        allow_none: bool = False,
    ) -> object:
        self.calls.append((value, field_name, allow_none))
        return value


class _FalseValuedPreprocessingValidator:
    def __init__(self) -> None:
        self.calls = 0

    def __bool__(self) -> bool:
        return False

    def run(self, value: object) -> object:
        self.calls += 1
        return value


class _RaisingPreprocessingValidator:
    def run(self, value: object) -> object:
        raise _InjectedFailure("preprocessing")


def _valid_dataset_request(
    *,
    preprocessing_config: DatasetPreprocessingConfig | None = None,
    total: pd.DataFrame | None = None,
) -> DatasetBuildRequest:
    return DatasetBuildRequest(
        phospho=pd.DataFrame({"sample": [1.0]}),
        site_metadata=pd.DataFrame({"site": ["S1"]}),
        preprocessing_config=(
            DatasetPreprocessingConfig()
            if preprocessing_config is None
            else preprocessing_config
        ),
        total=total,
    )


@pytest.fixture(
    params=(
        (DatasetInputSourceValidator, UnsupportedInputFormatError),
        (_DatasetInputSourceValidator, PhosPyInputError),
        (_PhosphositeImportSourceValidator, UnsupportedInputFormatError),
    ),
    ids=("dataset-request", "construction-reader", "importer"),
)
def source_adapter(
    request: pytest.FixtureRequest,
) -> tuple[Any, type[PhosPyInputError]]:
    adapter_type, error_type = request.param
    return adapter_type(), error_type


def test_source_adapters_preserve_dataframe_string_and_path_identity(
    source_adapter: tuple[Any, type[PhosPyInputError]],
) -> None:
    adapter, _error_type = source_adapter
    frame = pd.DataFrame({"sample": [1.0]})
    string_path = "".join(("table", ".csv"))
    path = Path("table.csv")

    assert adapter.run(frame, field_name="phospho") is frame
    assert adapter.run(string_path, field_name="phospho") is string_path
    assert adapter.run(path, field_name="phospho") is path


def test_source_adapters_convert_generic_pathlike_to_path(
    source_adapter: tuple[Any, type[PhosPyInputError]],
) -> None:
    adapter, _error_type = source_adapter

    resolved = adapter.run(_GenericPathLike(), field_name="phospho")

    assert isinstance(resolved, Path)
    assert resolved == Path("generic.csv")


def test_source_adapters_apply_none_and_allow_none_policy(
    source_adapter: tuple[Any, type[PhosPyInputError]],
) -> None:
    adapter, error_type = source_adapter

    assert adapter.run(None, field_name="phospho", allow_none=True) is None
    with pytest.raises(error_type, match="dataset build request phospho is required"):
        adapter.run(None, field_name="phospho")


@pytest.mark.parametrize("value", ("", "  \t", object()))
def test_source_adapters_reject_blank_and_unsupported_inputs_with_owned_error(
    source_adapter: tuple[Any, type[PhosPyInputError]],
    value: object,
) -> None:
    adapter, error_type = source_adapter
    expected = "path cannot be empty" if isinstance(value, str) else "must be a pandas"

    with pytest.raises(error_type, match=expected) as exc_info:
        adapter.run(value, field_name="phospho")

    assert type(exc_info.value) is error_type


def test_dataset_request_and_construction_adapters_keep_distinct_error_classes() -> (
    None
):
    with pytest.raises(UnsupportedInputFormatError) as request_error:
        DatasetInputSourceValidator().run(object(), field_name="phospho")
    with pytest.raises(PhosPyInputError) as construction_error:
        _DatasetInputSourceValidator().run(object(), field_name="phospho")

    assert type(request_error.value) is UnsupportedInputFormatError
    assert type(construction_error.value) is PhosPyInputError


def test_dataset_input_reader_propagates_injected_source_and_path_exceptions() -> None:
    source_reader = DatasetInputReader(
        source_validator=_RaisingSourceValidator(),  # type: ignore[arg-type]
    )
    with pytest.raises(_InjectedFailure, match="source:phospho"):
        source_reader.run(pd.DataFrame({"sample": [1.0]}), field_name="phospho")

    path_reader = DatasetInputReader(path_reader=_RaisingPathReader())
    with pytest.raises(_InjectedFailure, match=r"path:phospho:table\.csv"):
        path_reader.run("table.csv", field_name="phospho")


def test_dataset_input_reader_replaces_false_valued_validator_with_default() -> None:
    validator = _FalseValuedSourceValidator()
    reader = DatasetInputReader(source_validator=validator)  # type: ignore[arg-type]

    with pytest.raises(
        PhosPyInputError,
        match="phospho must be a pandas DataFrame or a file path",
    ):
        reader.run(object(), field_name="phospho")  # type: ignore[arg-type]

    assert validator.calls == []


def test_dataset_request_validator_replaces_false_valued_source_with_default() -> None:
    validator = _FalseValuedSourceValidator()
    request = DatasetBuildRequest(phospho=object(), site_metadata=pd.DataFrame())

    with pytest.raises(
        UnsupportedInputFormatError,
        match="phospho must be a pandas DataFrame or a file path",
    ):
        DatasetBuildRequestValidator(
            source_validator=validator,  # type: ignore[arg-type]
        ).run(request)

    assert validator.calls == []


def test_dataset_request_validator_replaces_false_valued_preprocessing_validator() -> (
    None
):
    validator = _FalseValuedPreprocessingValidator()
    invalid_config = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(policy="identity"),
        total_protein_correction=DatasetTotalProteinCorrectionConfig(
            policy="subtract_log_total"
        ),
    )
    request = _valid_dataset_request(
        preprocessing_config=invalid_config,
        total=pd.DataFrame({"sample": [2.0]}),
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "total_protein_correction.policy='subtract_log_total' requires "
            "log2-scale phospho and total values"
        ),
    ):
        DatasetBuildRequestValidator(
            preprocessing_validator=validator,  # type: ignore[arg-type]
        ).run(request)

    assert validator.calls == 0


def test_dataset_request_validator_propagates_truthy_injected_exceptions() -> None:
    request = _valid_dataset_request()

    with pytest.raises(_InjectedFailure, match="source:phospho"):
        DatasetBuildRequestValidator(
            source_validator=_RaisingSourceValidator(),  # type: ignore[arg-type]
        ).run(request)

    with pytest.raises(_InjectedFailure, match="preprocessing"):
        DatasetBuildRequestValidator(
            preprocessing_validator=_RaisingPreprocessingValidator(),  # type: ignore[arg-type]
        ).run(request)


def test_importer_validation_runs_source_then_mapping_then_column_fields() -> None:
    source_first = PhosphositeImportRequest(
        source=" ",
        sample_intensity_columns=(),
        gene_symbol_column="",
        source_name="",
    )
    with pytest.raises(
        UnsupportedInputFormatError, match="source path cannot be empty"
    ):
        PhosphositeImportRequestValidator().run(source_first)

    mapping_second = PhosphositeImportRequest(
        source=pd.DataFrame(),
        sample_intensity_columns=(),
        gene_symbol_column="",
        source_name="",
    )
    with pytest.raises(PhosPyInputError, match="must not be empty"):
        PhosphositeImportRequestValidator().run(mapping_second)


def test_importer_propagates_injected_source_exception_before_later_validation() -> (
    None
):
    request = PhosphositeImportRequest(
        source=pd.DataFrame(),
        sample_intensity_columns=(),
        gene_symbol_column="",
        source_name="",
    )
    validator = PhosphositeImportRequestValidator(
        source_validator=_RaisingSourceValidator(),  # type: ignore[arg-type]
    )

    with pytest.raises(_InjectedFailure, match="source:source"):
        validator.run(request)


def test_importer_replaces_false_valued_source_collaborator_with_default() -> None:
    source_validator = _FalseValuedSourceValidator()
    request = PhosphositeImportRequest(
        source=" ",
        sample_intensity_columns=("intensity",),
    )
    validator = PhosphositeImportRequestValidator(
        source_validator=source_validator,  # type: ignore[arg-type]
    )

    with pytest.raises(
        UnsupportedInputFormatError, match="source path cannot be empty"
    ):
        validator.run(request)

    assert source_validator.calls == []

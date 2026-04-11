from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import PhospyValidationError, TableSchemaError
from phospy.validation.schema.frames import require_dataframe
from phospy.validation.values.collections import resolve_required_columns
from phospy.validation.values.enums import validate_svm_mode
from phospy.validation.values.identifiers import require_splitable_gene_p_site
from phospy.validation.values.numeric import validate_fraction


def test_validate_fraction_rejects_non_numeric_values() -> None:
    with pytest.raises(PhospyValidationError, match="between 0 and 1"):
        validate_fraction("bad", name="threshold")  # type: ignore[arg-type]


def test_validate_svm_mode_accepts_supported_value() -> None:
    assert validate_svm_mode("default") == "default"


def test_resolve_required_columns_rejects_empty_collection() -> None:
    with pytest.raises(PhospyValidationError, match="requires at least one column"):
        resolve_required_columns([], argument_name="value_cols", context="example")


def test_require_dataframe_rejects_non_dataframe_inputs() -> None:
    with pytest.raises(TableSchemaError, match="must be a pandas DataFrame"):
        require_dataframe([1, 2, 3], context="example frame")  # type: ignore[arg-type]


def test_require_splitable_gene_p_site_rejects_malformed_identifiers() -> None:
    series = pd.Series(["GENE__S1"])

    with pytest.raises(TableSchemaError, match="single underscore"):
        require_splitable_gene_p_site(series, context="example identifiers")

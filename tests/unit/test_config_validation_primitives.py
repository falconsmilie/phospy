from __future__ import annotations

import pytest

from phospy.errors.input import PhosPyInputError
from phospy.policy_models import IntensityTransformPolicy
from phospy.validation.common.config_values import (
    coerce_policy_enum,
    require_supported_literal,
)
from phospy.validation.common.numbers import (
    require_optional_int_at_least,
    require_optional_real_between,
)
from phospy.validation.common.paths import require_local_filesystem_path
from phospy.validation.configs.preprocessing import (
    validate_preprocessing_section_type,
)


def test_require_supported_literal_reports_sorted_supported_values() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="preprocessing_config.normalisation.policy must be one of: median_center, none",
    ):
        require_supported_literal(
            "bad",
            field_name="preprocessing_config.normalisation.policy",
            supported_values=frozenset({"none", "median_center"}),
            error_type=PhosPyInputError,
        )


def test_coerce_policy_enum_trims_and_parses_supported_values() -> None:
    parsed = coerce_policy_enum(
        IntensityTransformPolicy,
        " log2 ",
        field_name="preprocessing_config.intensity_transform.policy",
        error_type=PhosPyInputError,
    )
    assert parsed is IntensityTransformPolicy.LOG2


def test_coerce_policy_enum_rejects_unsupported_values_with_field_name() -> None:
    with pytest.raises(
        PhosPyInputError,
        match=(
            "preprocessing_config.intensity_transform.policy must be one of: "
            "identity, log2; got 'bad'"
        ),
    ):
        coerce_policy_enum(
            IntensityTransformPolicy,
            "bad",
            field_name="preprocessing_config.intensity_transform.policy",
            error_type=PhosPyInputError,
        )


def test_require_local_filesystem_path_rejects_remote_urls() -> None:
    with pytest.raises(PhosPyInputError, match="must be a local filesystem path"):
        require_local_filesystem_path(
            "https://example.org/ref.fasta",
            field_name="preprocessing_config.site_sequence_resolution.fasta_path",
            error_type=PhosPyInputError,
            when_provided=True,
        )


@pytest.mark.parametrize(
    ("value", "minimum", "expected", "pattern"),
    [
        pytest.param(None, 1, None, None, id="none-allowed"),
        pytest.param(True, 1, None, "config.seed must be an int", id="wrong-type-bool"),
        pytest.param(
            0, 1, None, "config.seed must be greater than or equal to 1", id="zero"
        ),
        pytest.param(
            -1,
            1,
            None,
            "config.seed must be greater than or equal to 1",
            id="negative",
        ),
        pytest.param(5, 1, 5, None, id="valid-positive"),
    ],
)
def test_require_optional_int_at_least_matrix(
    value: object,
    minimum: int,
    expected: int | None,
    pattern: str | None,
) -> None:
    # Matrix keeps before/after coverage obvious for optional-positive integer cases.
    if pattern is None:
        assert (
            require_optional_int_at_least(
                value,
                field_name="config.seed",
                minimum=minimum,
                error_type=PhosPyInputError,
            )
            == expected
        )
        return

    with pytest.raises(PhosPyInputError, match=pattern):
        require_optional_int_at_least(
            value,
            field_name="config.seed",
            minimum=minimum,
            error_type=PhosPyInputError,
        )


def test_optional_numeric_primitives_accept_none_and_validate_ranges() -> None:
    assert (
        require_optional_real_between(
            None,
            field_name="config.threshold",
            minimum=0.0,
            maximum=1.0,
            error_type=PhosPyInputError,
        )
        is None
    )


def test_validate_preprocessing_section_type_reuses_field_specific_message() -> None:
    class DummyConfig:
        pass

    with pytest.raises(
        PhosPyInputError,
        match=(
            "dataset build request preprocessing_config.missing_data "
            "must be a DummyConfig"
        ),
    ):
        validate_preprocessing_section_type(
            "not-a-config",
            field_name="dataset build request preprocessing_config.missing_data",
            expected_type=DummyConfig,
        )

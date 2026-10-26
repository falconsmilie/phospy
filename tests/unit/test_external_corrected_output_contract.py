from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from phospy.errors.input import PhosPyInputError
from phospy.preprocessing.batch_correction.contracts.output import (
    _CALLER_ATTESTED_EXTERNAL_CORRECTION_ORIGIN,
    _CALLER_ATTESTED_EXTERNAL_CORRECTION_VERIFICATION,
    CorrectedPreprocessingOutput,
    create_external_corrected_preprocessing_output,
)
from phospy.preprocessing.batch_correction.validation.output import (
    validate_corrected_preprocessing_output,
)
from phospy.preprocessing.control_sites.contracts import (
    ControlSiteSet,
    ControlSiteSourceMetadata,
)
from phospy.provenance.models.batch_correction import BatchCorrectionProvenance


def _matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_1": (1.0, 2.0), "sample_2": (3.0, 4.0)},
        index=pd.Index(("site_1", "site_2"), name="site_key"),
    )


def _output() -> CorrectedPreprocessingOutput:
    before = _matrix()
    return create_external_corrected_preprocessing_output(
        pre_correction_matrix=before,
        corrected_matrix=before + 0.5,
        method="sps_ruv_style",
        source_software="external-tool",
        source_version="1.4",
        control_site_set=_control_set(before),
        n_unwanted_factors=1,
        batch_column="batch",
        condition_columns=("condition",),
        parameters={"centering": "external-tool-default"},
        diagnostics={"converged": True},
    )


def _control_set(matrix: pd.DataFrame) -> ControlSiteSet:
    return ControlSiteSet.from_site_keys(
        tuple(matrix.index.astype(str)),
        source_metadata=ControlSiteSourceMetadata(
            organism="rat",
            identifier_namespace="site_key",
            source_name="unit-test-controls",
            source_version="1",
            license="test-only",
            redistribution="not redistributed",
        ),
    )


def test_external_factory_builds_untrusted_typed_provenance() -> None:
    output = _output()

    validate_corrected_preprocessing_output(output)
    provenance = output.provenance
    assert isinstance(provenance, BatchCorrectionProvenance)
    assert (
        provenance.resolved_parameters["execution_origin"]
        == _CALLER_ATTESTED_EXTERNAL_CORRECTION_ORIGIN
    )
    assert (
        provenance.resolved_parameters["execution_verification"]
        == _CALLER_ATTESTED_EXTERNAL_CORRECTION_VERIFICATION
    )
    assert provenance.selected_site_key_rows == ("site_1", "site_2")
    assert provenance.control_site_source["selection_performed_by_phospy"] is False
    assert provenance.phospy_version == "external_execution_not_observed_by_phospy"
    assert provenance.python_version == "external_execution_not_observed_by_phospy"
    assert output.batch_correction_report.confounding_check_status == "not_checked"


def test_external_validation_rejects_inconsistent_provenance() -> None:
    output = _output()
    provenance = output.provenance
    assert isinstance(provenance, BatchCorrectionProvenance)
    inconsistent = replace(
        output,
        provenance=replace(provenance, requested_method="ruv_iii_style"),
    )

    with pytest.raises(PhosPyInputError, match="provenance method must match"):
        validate_corrected_preprocessing_output(inconsistent)


@pytest.mark.parametrize("method", ("none", "combat", "custom"))
def test_external_factory_rejects_unregistered_methods(method: str) -> None:
    before = _matrix()

    with pytest.raises(PhosPyInputError, match="method"):
        create_external_corrected_preprocessing_output(
            pre_correction_matrix=before,
            corrected_matrix=before + 0.5,
            method=method,
            source_software="external-tool",
            source_version="1.4",
            control_site_set=_control_set(before),
            n_unwanted_factors=1,
            batch_column="batch",
            condition_columns=("condition",),
        )


def test_external_factory_rejects_nonfinite_corrected_values() -> None:
    before = _matrix()
    corrected = before.copy()
    corrected.iat[0, 0] = np.nan

    with pytest.raises(PhosPyInputError, match="finite and missing-value-free"):
        create_external_corrected_preprocessing_output(
            pre_correction_matrix=before,
            corrected_matrix=corrected,
            method="sps_ruv_style",
            source_software="external-tool",
            source_version="1.4",
            control_site_set=_control_set(before),
            n_unwanted_factors=1,
            batch_column="batch",
            condition_columns=("condition",),
        )


def test_external_factory_enforces_control_count_for_unwanted_factors() -> None:
    before = _matrix()
    one_control = ControlSiteSet.from_site_keys(
        (str(before.index[0]),),
        source_metadata=_control_set(before).source_metadata,
    )

    with pytest.raises(PhosPyInputError, match="eligible controls"):
        create_external_corrected_preprocessing_output(
            pre_correction_matrix=before,
            corrected_matrix=before + 0.5,
            method="sps_ruv_style",
            source_software="external-tool",
            source_version="1.4",
            control_site_set=one_control,
            n_unwanted_factors=1,
            batch_column="batch",
            condition_columns=("condition",),
        )


def test_external_ruv_iii_requires_replicate_column() -> None:
    before = _matrix()

    with pytest.raises(PhosPyInputError, match="replicate_column is required"):
        create_external_corrected_preprocessing_output(
            pre_correction_matrix=before,
            corrected_matrix=before + 0.5,
            method="ruv_iii_style",
            source_software="external-tool",
            source_version="1.4",
            control_site_set=_control_set(before),
            n_unwanted_factors=1,
            batch_column="batch",
            condition_columns=("condition",),
        )

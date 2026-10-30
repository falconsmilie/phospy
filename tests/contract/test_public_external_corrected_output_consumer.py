from __future__ import annotations

from typing import get_type_hints

import pandas as pd
import pytest

from phospy.api import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    Organism,
    PhosPyInputError,
)
from phospy.api.advanced import (
    ControlSiteSet,
    ControlSiteSourceMetadata,
    CorrectedPreprocessingOutput,
    DatasetBatchCorrectionConfig,
    create_external_corrected_preprocessing_output,
)

pytestmark = pytest.mark.contract


def _site_key(protein: str, residue: str, position: int) -> str:
    return (
        "phospy:v1|organism=rat|protein_namespace=protein_id|"
        f"protein_identifier={protein}|residue={residue}|position={position}"
    )


def _consumer_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    site_keys = pd.Index(
        (
            _site_key("MAPK14", "Y", 182),
            _site_key("SRC", "Y", 416),
        ),
        name="site_key",
    )
    phospho = pd.DataFrame(
        {
            "sample_1": (10.0, 20.0),
            "sample_2": (12.0, 24.0),
        },
        index=site_keys,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ("MAPK14", "SRC"),
            "protein_id": ("MAPK14", "SRC"),
            "site": ("Y182", "Y416"),
            "site_sequence": (
                ("A" * 15) + "Y" + ("A" * 15),
                ("A" * 15) + "Y" + ("A" * 15),
            ),
            "localisation_confidence": (0.95, 0.98),
        },
        index=site_keys.copy(),
    )
    return phospho, site_metadata


def _external_output(
    phospho: pd.DataFrame,
) -> CorrectedPreprocessingOutput:
    corrected = phospho.subtract(phospho.mean(axis="columns"), axis="index")
    controls = ControlSiteSet.from_site_keys(
        tuple(phospho.index.astype(str)),
        source_metadata=ControlSiteSourceMetadata(
            organism="rat",
            identifier_namespace="site_key",
            source_name="consumer-controls",
            source_version="1",
            license="caller local use",
            redistribution="not redistributed",
        ),
    )
    return create_external_corrected_preprocessing_output(
        pre_correction_matrix=phospho,
        corrected_matrix=corrected,
        method="sps_ruv_style",
        source_software="consumer-correction-tool",
        source_version="2.1",
        control_site_set=controls,
        n_unwanted_factors=1,
        batch_column="batch",
        condition_columns=("condition",),
        parameters={"centering": "external-tool-default"},
        diagnostics={"convergence": "reported_by_external_tool"},
    )


def test_public_consumer_creates_and_submits_external_corrected_output() -> None:
    phospho, site_metadata = _consumer_tables()
    output = _external_output(phospho)

    annotation = get_type_hints(DatasetBuildRequest)["corrected_preprocessing_output"]
    assert annotation == CorrectedPreprocessingOutput | None
    assert isinstance(output, CorrectedPreprocessingOutput)
    assert CorrectedPreprocessingOutput.__module__ == (
        "phospy.preprocessing.batch_correction.contracts.output"
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig.strict(),
            input_intensity_scale="log2",
            corrected_preprocessing_output=output,
        )
    )

    pd.testing.assert_frame_equal(dataset.phospho, output.corrected_matrix)
    assert dataset.preprocessing_report is not None
    assert dataset.preprocessing_report.batch_correction is not None
    assert dataset.preprocessing_report.batch_correction.status == "applied"
    assert (
        output.diagnostics["execution_origin"] == "caller_supplied_external_correction"
    )
    assert (
        output.diagnostics["execution_verification"]
        == "caller_attested_not_verified_by_phospy"
    )


def test_public_factory_rejects_incompatible_axes() -> None:
    phospho, _ = _consumer_tables()
    controls = ControlSiteSet.from_site_keys(
        tuple(phospho.index.astype(str)),
        source_metadata=ControlSiteSourceMetadata(
            organism="rat",
            identifier_namespace="site_key",
            source_name="consumer-controls",
            source_version="1",
            license="caller local use",
            redistribution="not redistributed",
        ),
    )

    with pytest.raises(PhosPyInputError, match="corrected_matrix.index must match"):
        create_external_corrected_preprocessing_output(
            pre_correction_matrix=phospho,
            corrected_matrix=phospho.iloc[::-1],
            method="sps_ruv_style",
            source_software="consumer-correction-tool",
            source_version="2.1",
            control_site_set=controls,
            n_unwanted_factors=1,
            batch_column="batch",
            condition_columns=("condition",),
        )


def test_public_builder_rejects_external_output_for_different_input_state() -> None:
    phospho, site_metadata = _consumer_tables()
    output = _external_output(phospho)
    changed = phospho.copy()
    changed.iat[0, 0] = 11.0

    with pytest.raises(PhosPyInputError, match="fingerprint mismatch"):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=changed,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig.strict(),
                input_intensity_scale="log2",
                corrected_preprocessing_output=output,
            )
        )


def test_public_builder_rejects_double_correction_before_interpretation() -> None:
    phospho, site_metadata = _consumer_tables()
    output = _external_output(phospho)

    with pytest.raises(PhosPyInputError, match="applied exactly once"):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig(
                    batch_correction=DatasetBatchCorrectionConfig(
                        method="linear_residualize_batch",
                    )
                ),
                input_intensity_scale="log2",
                corrected_preprocessing_output=output,
            )
        )

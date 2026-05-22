from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import DatasetBuildRequest
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.builders.contracts import (
    InterpretedDatasetBuildRequest,
    PreprocessedDatasetBuildTables,
)
from phospy.science.datasets.builders.interpreter_collaborators import (
    DatasetBuildSourceResolver,
)
from phospy.science.datasets.builders.normalizer import NormalizedDatasetInputs
from phospy.science.datasets.builders.provenance_assembler import (
    DatasetRunProvenanceAssembler,
)
from phospy.science.datasets.builders.site_sequence_boundary import (
    AnalysisReadySiteSequenceValidator,
)
from phospy.science.datasets.preprocessing.models import PreprocessingPlan


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )


def _site_metadata(include_sequence: bool = True) -> pd.DataFrame:
    data: dict[str, object] = {
        "gene_symbol": ["MAPK14"],
        "site": ["Y182"],
    }
    if include_sequence:
        data["site_sequence"] = ["AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"]
    return pd.DataFrame(data, index=_phospho().index.copy())


def test_source_resolver_reads_and_normalizes_site_level_inputs() -> None:
    reader_calls: list[str] = []
    normalizer_calls: list[tuple[pd.DataFrame, pd.DataFrame]] = []
    phospho = _phospho()
    site_metadata = _site_metadata()
    sample_metadata = pd.DataFrame(
        {"group": ["A"]},
        index=pd.Index(["sample_a"], name="sample"),
    )

    class ReaderSpy:
        def run(self, value: object, *, field_name: str) -> pd.DataFrame:
            reader_calls.append(field_name)
            assert isinstance(value, pd.DataFrame)
            return value

    class NormalizerSpy:
        def run(
            self,
            *,
            phospho: pd.DataFrame,
            site_metadata: pd.DataFrame,
            sample_metadata: pd.DataFrame | None,
            total: pd.DataFrame | None,
        ) -> NormalizedDatasetInputs:
            assert sample_metadata is not None
            assert total is None
            normalizer_calls.append((phospho, site_metadata))
            return NormalizedDatasetInputs(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=None,
                site_identifier_normalisation=None,
            )

    resolved = DatasetBuildSourceResolver(
        reader=ReaderSpy(),
        normalizer=NormalizerSpy(),
    ).run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            input_intensity_scale="linear",
        )
    )

    assert reader_calls == ["sample_metadata", "phospho", "site_metadata"]
    assert len(normalizer_calls) == 1
    assert resolved.site_resolution_mode == "site_level_resolved"
    assert resolved.multi_site_policy is None
    assert resolved.peptide_evidence_resolution is None


def test_site_sequence_boundary_validator_rejects_missing_sequence_column() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="site_sequence is missing, blank, or invalid after builder enrichment",
    ):
        AnalysisReadySiteSequenceValidator().run(
            site_metadata=_site_metadata(include_sequence=False),
            preprocessing_trace=None,
        )


def test_run_provenance_assembler_records_opaque_token_mode() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    plan = PreprocessingPlan.default()
    request = InterpretedDatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        organism=None,
        preprocessing_plan=plan,
        allow_opaque_site_values=True,
    )
    preprocessed = PreprocessedDatasetBuildTables(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        preprocessing_trace=None,
    )

    provenance = DatasetRunProvenanceAssembler().run(
        request=request,
        preprocessed=preprocessed,
        validated_site_metadata=site_metadata,
        resolved_phospho=phospho,
        resolved_total=None,
        preprocessing_trace=None,
        intensity_scale_label="linear",
        intensity_scale_establishment={"source": "test"},
        quantitative_meaning="phosphosite_abundance",
        allow_opaque_site_values=True,
    )

    assert provenance.workflow_parameters["site_token_validation"] == {
        "mode": "opaque_opt_in"
    }

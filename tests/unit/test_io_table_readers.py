from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.api.requests import DatasetBuildRequest
from phospy.errors import PhosPyInputError
from phospy.io.readers.tables import (
    read_contrast_matrix,
    read_design_matrix,
    read_phospho_matrix,
    read_sample_metadata,
    read_site_metadata,
)
from phospy.science.datasets.builders.public import AnalysisReadyDatasetBuilder


def _write_table(path: Path, frame: pd.DataFrame) -> None:
    if path.suffix == ".tsv":
        frame.to_csv(path, sep="\t")
        return
    frame.to_csv(path)


@pytest.mark.parametrize("suffix", (".csv", ".tsv"))
def test_site_metadata_reader_preserves_na_like_identifiers(
    suffix: str, tmp_path: Path
) -> None:
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["NA", "P00001"],
            "site": ["S1", "T308"],
            "site_sequence": ["A" * 31, "B" * 31],
            "protein_id": ["00123", "1E10"],
            "localisation_confidence": [0.95, 0.9],
        },
        index=pd.Index(["NA", "P00001"], name="site_id"),
    )
    path = tmp_path / f"site_metadata{suffix}"
    _write_table(path, site_metadata)

    read_back = read_site_metadata(path)

    assert read_back.loc["NA", "gene_symbol"] == "NA"
    assert read_back.loc["NA", "protein_id"] == "00123"
    assert read_back.loc["P00001", "protein_id"] == "1E10"
    assert list(read_back.index) == ["NA", "P00001"]
    assert read_back.loc["NA", "gene_symbol"] is not pd.NA


@pytest.mark.parametrize("suffix", (".csv", ".tsv"))
def test_sample_metadata_reader_preserves_numeric_looking_identifiers(
    suffix: str, tmp_path: Path
) -> None:
    sample_metadata = pd.DataFrame(
        {
            "comparison_group": ["01", "02", "03"],
            "batch": ["1E10", "7", "8"],
        },
        index=pd.Index(["01", "02", "03"], name="sample_id"),
    )
    path = tmp_path / f"sample_metadata{suffix}"
    _write_table(path, sample_metadata)

    read_back = read_sample_metadata(path)

    assert list(read_back.index) == ["01", "02", "03"]
    assert read_back.loc["01", "comparison_group"] == "01"
    assert read_back.loc["01", "batch"] == "1E10"
    assert read_back.map(lambda value: isinstance(value, str)).all().all()


@pytest.mark.parametrize("suffix", (".csv", ".tsv"))
def test_phospho_matrix_reader_parses_numeric_values(
    suffix: str, tmp_path: Path
) -> None:
    phospho = pd.DataFrame(
        {"01": [1.0, 2.5], "02": [3.0, 4.25]},
        index=pd.Index(["P00001;S1;", "P00002;T2;"], name="site_id"),
    )
    path = tmp_path / f"phospho{suffix}"
    _write_table(path, phospho)

    read_back = read_phospho_matrix(path)

    assert pd.api.types.is_float_dtype(read_back["01"])
    assert read_back.loc["P00001;S1;", "01"] == pytest.approx(1.0)
    assert read_back.loc["P00002;T2;", "02"] == pytest.approx(4.25)


def test_phospho_matrix_reader_rejects_invalid_numeric_cells(tmp_path: Path) -> None:
    phospho = pd.DataFrame(
        {"sample_a": ["3.14"], "sample_b": ["not_a_number"]},
        index=pd.Index(["P00001;S1;"], name="site_id"),
    )
    path = tmp_path / "phospho.csv"
    phospho.to_csv(path)

    with pytest.raises(
        PhosPyInputError,
        match=(
            "failed to parse numeric cell: .*"
            "table_role='phospho_matrix'.*"
            "row_label='P00001;S1;'.*"
            "column_label='sample_b'.*"
            "offending_value='not_a_number'.*"
            "expected_type='finite numeric value or allowed missing marker"
        ),
    ):
        read_phospho_matrix(path)


def test_dataset_build_from_file_paths_uses_schema_aware_readers(
    tmp_path: Path,
) -> None:
    phospho = pd.DataFrame(
        {"01": [10.0]},
        index=pd.Index(["NA;S1;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["NA"],
            "site": ["S1"],
            "site_sequence": ["A" * 31],
            "localisation_confidence": [0.95],
        },
        index=pd.Index(["NA;S1;"], name="site_id"),
    )
    sample_metadata = pd.DataFrame(
        {"comparison_group": ["01"]},
        index=pd.Index(["01"], name="sample_id"),
    )
    phospho_path = tmp_path / "phospho.csv"
    site_metadata_path = tmp_path / "site_metadata.csv"
    sample_metadata_path = tmp_path / "sample_metadata.csv"
    phospho.to_csv(phospho_path)
    site_metadata.to_csv(site_metadata_path)
    sample_metadata.to_csv(sample_metadata_path)

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho_path,
            site_metadata=site_metadata_path,
            sample_metadata=sample_metadata_path,
        )
    )

    assert built.site_metadata.loc["NA;S1;", "gene_symbol"] == "NA"
    assert list(built.sample_metadata.index) == ["01"]
    assert built.sample_metadata.loc["01", "comparison_group"] == "01"


def test_dataset_build_from_file_paths_reports_invalid_numeric_cell_context(
    tmp_path: Path,
) -> None:
    phospho = pd.DataFrame(
        {"sample_a": ["x"]},
        index=pd.Index(["P00001;S1;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["P00001"],
            "site": ["S1"],
            "site_sequence": ["A" * 31],
            "localisation_confidence": [0.95],
        },
        index=pd.Index(["P00001;S1;"], name="site_id"),
    )
    phospho_path = tmp_path / "phospho.csv"
    site_metadata_path = tmp_path / "site_metadata.csv"
    phospho.to_csv(phospho_path)
    site_metadata.to_csv(site_metadata_path)

    with pytest.raises(
        PhosPyInputError,
        match=(
            "failed to read dataset build request phospho from .*"
            "table_role='phospho_matrix'.*"
            "row_label='P00001;S1;'.*"
            "column_label='sample_a'.*"
            "offending_value='x'"
        ),
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho_path,
                site_metadata=site_metadata_path,
            )
        )


def test_design_matrix_reader_rejects_missing_tokens(tmp_path: Path) -> None:
    design = pd.DataFrame(
        {"intercept": ["1"], "treatment": ["NA"]},
        index=pd.Index(["01"], name="sample_id"),
    )
    path = tmp_path / "design.csv"
    design.to_csv(path)

    with pytest.raises(
        PhosPyInputError,
        match=(
            "table_role='design_matrix'.*"
            "row_label='01'.*"
            "column_label='treatment'.*"
            "offending_value='NA'.*"
            "expected_type='finite numeric value'"
        ),
    ):
        read_design_matrix(path)


def test_contrast_matrix_reader_parses_finite_numeric_values(tmp_path: Path) -> None:
    contrasts = pd.DataFrame(
        {"B_vs_A": [1.0, -1.0]},
        index=pd.Index(["intercept", "treatment"], name="coefficient"),
    )
    path = tmp_path / "contrasts.csv"
    contrasts.to_csv(path)

    read_back = read_contrast_matrix(path)
    assert read_back.loc["intercept", "B_vs_A"] == pytest.approx(1.0)
    assert read_back.loc["treatment", "B_vs_A"] == pytest.approx(-1.0)

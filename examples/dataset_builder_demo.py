#!/usr/bin/env python3
"""Build analysis-ready datasets and show site-matrix row-retention behavior."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DatasetBuildRequest,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    Organism,
    UnsupportedInputFormatError,
)


def _example_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 0.7],
            "sample_b": [1.2, 0.8],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
            ],
        },
        index=phospho.index.copy(),
    )
    return phospho, site_metadata


def build_from_dataframes() -> AnalysisReadyPhosphoDataset:
    phospho, site_metadata = _example_tables()
    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        preprocessing_config=DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(policy="forbid"),
        ),
    )
    return AnalysisReadyDatasetBuilder().run(request)


def build_from_file_paths() -> AnalysisReadyPhosphoDataset:
    phospho, site_metadata = _example_tables()
    with TemporaryDirectory(prefix="phospy-builder-demo-") as tmp_dir:
        root = Path(tmp_dir)
        phospho_path = root / "phospho.csv"
        site_metadata_path = root / "site_metadata.csv"
        phospho.to_csv(phospho_path)
        site_metadata.to_csv(site_metadata_path)
        request = DatasetBuildRequest(
            phospho=phospho_path,
            site_metadata=str(site_metadata_path),
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=2,
                ),
            ),
        )
        return AnalysisReadyDatasetBuilder().run(request)


def build_with_site_matrix_from_metadata() -> tuple[
    AnalysisReadyPhosphoDataset,
    list[str],
    list[str],
    str,
]:
    def _request_for(
        phospho_frame: pd.DataFrame, site_metadata_frame: pd.DataFrame
    ) -> DatasetBuildRequest:
        return DatasetBuildRequest(
            phospho=phospho_frame,
            site_metadata=site_metadata_frame,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(policy="forbid"),
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    missing_data_policy="retain_missing",
                ),
            ),
        )

    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 0.7, 0.5],
            "sample_b": [1.2, 0.8, 0.4],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "",
                "RPHFPQFSYSASSTA",
            ],
        },
        index=phospho.index.copy(),
    )
    builder = AnalysisReadyDatasetBuilder()
    try:
        dataset = builder.run(_request_for(phospho, site_metadata))
        return dataset, list(phospho.index), [], ""
    except UnsupportedInputFormatError as error:
        usable_sequence = (
            site_metadata.loc[:, "site_sequence"].astype("string").str.strip() != ""
        )
        excluded_site_ids = site_metadata.index[~usable_sequence].astype(str).tolist()
        filtered_phospho = phospho.loc[usable_sequence]
        filtered_site_metadata = site_metadata.loc[usable_sequence]
        dataset = builder.run(_request_for(filtered_phospho, filtered_site_metadata))
        return dataset, list(phospho.index), excluded_site_ids, str(error)


def main() -> None:
    df_dataset = build_from_dataframes()
    path_dataset = build_from_file_paths()
    (
        site_matrix_dataset,
        site_matrix_input_ids,
        excluded_site_ids,
        rejection_message,
    ) = build_with_site_matrix_from_metadata()
    print("Dataset builder demo")
    print("DataFrame route phospho shape:", df_dataset.phospho.shape)
    print("File-path route phospho shape:", path_dataset.phospho.shape)
    print("Site metadata columns:", list(path_dataset.site_metadata.columns))
    print(
        "Organism:",
        None if path_dataset.organism is None else path_dataset.organism.value,
    )
    print(
        "Transformation state:",
        f"{path_dataset.transformation_state.label} (builder pass-through lane)",
    )
    print()
    print("build_from_metadata row-retention demo")
    print("Input phospho rows:", len(site_matrix_input_ids))
    print("Retained phospho rows:", site_matrix_dataset.phospho.shape[0])
    print("Rows without usable site_sequence:", excluded_site_ids)
    print(
        "Site-matrix row-drop stats:",
        site_matrix_dataset.phospho.attrs.get("site_matrix_row_drop_stats"),
    )
    if rejection_message:
        print("Builder boundary check for unusable site_sequence:", rejection_message)


if __name__ == "__main__":
    main()

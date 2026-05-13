from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.api.builders import AnalysisReadyDatasetBuilder
from phospy.api.configs import (
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetSiteSequenceResolutionConfig,
)
from phospy.api.requests import DatasetBuildRequest
from phospy.errors.input import PhosPyInputError
from phospy.science.references.models import Organism


def _phospho(index: pd.Index | None = None) -> pd.DataFrame:
    site_index = (
        pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id")
        if index is None
        else index
    )
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=site_index.copy(),
    )


def _site_metadata(
    *,
    site_sequences: list[object] | None = None,
    with_sequence_column: bool = True,
) -> pd.DataFrame:
    data: dict[str, object] = {
        "gene_symbol": ["MAPK14", "GSK3B"],
        "site": ["Y182", "S9"],
        "protein_accession": ["P1", "P2"],
        "localisation_confidence": [0.95, 0.9],
    }
    if with_sequence_column:
        if site_sequences is None:
            site_sequences = ["SEQ_A", "SEQ_B"]
        data["site_sequence"] = site_sequences
    return pd.DataFrame(data, index=_phospho().index.copy())


def _write_fasta(tmp_path: Path) -> str:
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_text(
        ">P1 protein_1\nAAAAYAAAA\n>P2 protein_2\nAAAAASAAA\n",
        encoding="utf-8",
    )
    return str(fasta_path)


def _summary_for_built_dataset(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    organism: Organism | None = None,
    preprocessing_config: DatasetPreprocessingConfig | None = None,
):
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=organism,
            preprocessing_config=(
                DatasetPreprocessingConfig()
                if preprocessing_config is None
                else preprocessing_config
            ),
        )
    )
    assert built.preprocessing_report is not None
    summary = built.preprocessing_report.site_sequence_resolution_summary()
    assert summary is not None
    return summary


def test_sequence_resolution_report_counts_all_sequences_provided_by_input() -> None:
    summary = _summary_for_built_dataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(site_sequences=["SEQ_A", "SEQ_B"]),
    )

    assert summary.total_sites == 2
    assert summary.provided_by_input == 2
    assert summary.resolved_from_reference == 0
    assert summary.resolved_from_fasta == 0
    assert summary.unresolved == 0
    assert summary.conflicts == 0
    assert summary.conflict_policy == "not_applied"
    assert summary.final_sequence_complete_sites == 2


def test_sequence_resolution_report_counts_all_sequences_resolved_from_reference() -> (
    None
):
    summary = _summary_for_built_dataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(with_sequence_column=False),
        organism=Organism.RAT,
    )

    assert summary.total_sites == 2
    assert summary.provided_by_input == 0
    assert summary.resolved_from_reference == 2
    assert summary.resolved_from_fasta == 0
    assert summary.unresolved == 0
    assert summary.conflicts == 0
    assert summary.conflict_policy == "not_applied"
    assert summary.final_sequence_complete_sites == 2


def test_sequence_resolution_report_counts_mixed_input_and_reference_resolution() -> (
    None
):
    summary = _summary_for_built_dataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(site_sequences=["SEQ_A", pd.NA]),
        organism=Organism.RAT,
        preprocessing_config=DatasetPreprocessingConfig(
            site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
        ),
    )

    assert summary.total_sites == 2
    assert summary.provided_by_input == 1
    assert summary.resolved_from_reference == 1
    assert summary.resolved_from_fasta == 0
    assert summary.unresolved == 0
    assert summary.conflicts == 0
    assert summary.conflict_policy == "not_applied"
    assert summary.final_sequence_complete_sites == 2


def test_sequence_resolution_report_counts_all_sequences_resolved_from_fasta(
    tmp_path: Path,
) -> None:
    summary = _summary_for_built_dataset(
        phospho=_phospho(pd.Index(["MAPK14;Y5;", "GSK3B;S6;"], name="site_id")),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "GSK3B"],
                "site": ["Y5", "S6"],
                "protein_accession": ["P1", "P2"],
                "site_sequence": [pd.NA, pd.NA],
                "localisation_confidence": [0.95, 0.9],
            },
            index=pd.Index(["MAPK14;Y5;", "GSK3B;S6;"], name="site_id"),
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
                fasta_path=_write_fasta(tmp_path),
                flank_size=2,
            )
        ),
    )

    assert summary.total_sites == 2
    assert summary.provided_by_input == 0
    assert summary.resolved_from_reference == 0
    assert summary.resolved_from_fasta == 2
    assert summary.unresolved == 0
    assert summary.conflicts == 0
    assert summary.conflict_policy == "preserve_existing"
    assert summary.final_sequence_complete_sites == 2


def test_sequence_resolution_report_counts_conflicts_and_records_policy(
    tmp_path: Path,
) -> None:
    summary = _summary_for_built_dataset(
        phospho=_phospho(pd.Index(["MAPK14;Y5;", "GSK3B;S6;"], name="site_id")),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "GSK3B"],
                "site": ["Y5", "S6"],
                "protein_accession": ["P1", "P2"],
                "site_sequence": ["XXXXX", pd.NA],
                "localisation_confidence": [0.95, 0.9],
            },
            index=pd.Index(["MAPK14;Y5;", "GSK3B;S6;"], name="site_id"),
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
                fasta_path=_write_fasta(tmp_path),
                mode="validate_existing_and_fill_missing",
                conflict_policy="preserve_existing",
                flank_size=2,
            )
        ),
    )

    assert summary.total_sites == 2
    assert summary.provided_by_input == 1
    assert summary.resolved_from_reference == 0
    assert summary.resolved_from_fasta == 1
    assert summary.unresolved == 0
    assert summary.conflicts == 1
    assert summary.conflict_policy == "preserve_existing"
    assert summary.final_sequence_complete_sites == 2


def test_unresolved_sequence_blocks_dataset_creation_when_fasta_resolution_fails(
    tmp_path: Path,
) -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["FAKE1;S5;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["FAKE1"],
            "site": ["S5"],
            "protein_accession": ["P404"],
            "site_sequence": [pd.NA],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="site_sequence is missing, blank, or invalid after builder enrichment",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                preprocessing_config=DatasetPreprocessingConfig(
                    site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
                        fasta_path=_write_fasta(tmp_path),
                        mode="fill_missing_only",
                        flank_size=2,
                    )
                ),
            )
        )

from __future__ import annotations

from pathlib import Path

import pandas as pd

from phospy.api.configs import (
    DatasetPreprocessingConfig,
    DatasetSiteSequenceResolutionConfig,
)
from phospy.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.datasets.preprocessing.models import PreprocessingPlan


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [1.5, 2.5]},
        index=pd.Index(["MAPK14;S5;", "GSK3B;T6;"], name="site_id"),
    )


def _site_metadata(*, site_sequences: list[object]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["S5", "T6"],
            "protein_accession": ["P1", "P2"],
            "site_sequence": site_sequences,
        },
        index=_phospho().index.copy(),
    )


def _write_fasta(path: Path) -> str:
    fasta_path = path / "proteins.fasta"
    fasta_path.write_text(
        ">P1 protein_1\nAAAASAAAA\n>P2 protein_2\nCCCCCTCCCC\n",
        encoding="utf-8",
    )
    return str(fasta_path)


def _stage_diagnostics(preprocessed) -> dict[str, object]:
    trace = preprocessed.preprocessing_trace
    assert trace is not None
    matching = [item for item in trace if item.stage == "site_sequence_resolution"]
    assert len(matching) == 1
    return dict(matching[0].diagnostics)


def test_fasta_resolution_disabled_preserves_existing_behavior() -> None:
    preprocessed = DatasetPreprocessor().run(
        phospho=_phospho(),
        site_metadata=_site_metadata(site_sequences=[pd.NA, "CCTCC"]),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(DatasetPreprocessingConfig()),
    )

    assert "site_sequence_resolution" not in {
        item.stage for item in (preprocessed.preprocessing_trace or ())
    }
    assert pd.isna(preprocessed.site_metadata.loc["MAPK14;S5;", "site_sequence"])
    assert preprocessed.site_metadata.loc["GSK3B;T6;", "site_sequence"] == "CCTCC"


def test_fasta_resolution_fills_missing_and_validates_existing(tmp_path: Path) -> None:
    config = DatasetPreprocessingConfig(
        site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
            fasta_path=_write_fasta(tmp_path),
            flank_size=2,
        )
    )
    preprocessed = DatasetPreprocessor().run(
        phospho=_phospho(),
        site_metadata=_site_metadata(site_sequences=[pd.NA, "CCTCC"]),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(config),
    )

    assert preprocessed.site_metadata.loc["MAPK14;S5;", "site_sequence"] == "AASAA"
    assert preprocessed.site_metadata.loc["GSK3B;T6;", "site_sequence"] == "CCTCC"
    diagnostics = _stage_diagnostics(preprocessed)
    assert diagnostics["resolved_site_count"] == 2
    assert diagnostics["unresolved_site_count"] == 0
    assert diagnostics["filled_missing_count"] == 1
    assert diagnostics["preserved_existing_count"] == 1


def test_conflicting_existing_sequence_is_preserved_by_default(tmp_path: Path) -> None:
    config = DatasetPreprocessingConfig(
        site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
            fasta_path=_write_fasta(tmp_path),
            flank_size=2,
        )
    )
    preprocessed = DatasetPreprocessor().run(
        phospho=_phospho().iloc[:1, :].copy(deep=True),
        site_metadata=_site_metadata(site_sequences=["XXXXX", "CCTCC"])
        .iloc[:1, :]
        .copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(config),
    )

    assert preprocessed.site_metadata.loc["MAPK14;S5;", "site_sequence"] == "XXXXX"
    diagnostics = _stage_diagnostics(preprocessed)
    assert diagnostics["existing_sequence_conflict_count"] == 1
    assert diagnostics["replaced_existing_count"] == 0


def test_replace_existing_mode_replaces_conflicting_sequence(tmp_path: Path) -> None:
    config = DatasetPreprocessingConfig(
        site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
            fasta_path=_write_fasta(tmp_path),
            mode="replace_existing",
            flank_size=2,
        )
    )
    preprocessed = DatasetPreprocessor().run(
        phospho=_phospho().iloc[:1, :].copy(deep=True),
        site_metadata=_site_metadata(site_sequences=["XXXXX", "CCTCC"])
        .iloc[:1, :]
        .copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(config),
    )

    assert preprocessed.site_metadata.loc["MAPK14;S5;", "site_sequence"] == "AASAA"
    diagnostics = _stage_diagnostics(preprocessed)
    assert diagnostics["existing_sequence_conflict_count"] == 1
    assert diagnostics["replaced_existing_count"] == 1


def test_unresolved_reason_counts_are_reported(tmp_path: Path) -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0, 3.0]},
        index=pd.Index(["A;S5;", "B;S5;", "C;T5;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["A", "B", "C"],
            "site": ["S5", "S5", "T5"],
            "protein_accession": [pd.NA, "P404", "P1"],
            "site_sequence": [pd.NA, pd.NA, pd.NA],
        },
        index=phospho.index.copy(),
    )
    config = DatasetPreprocessingConfig(
        site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
            fasta_path=_write_fasta(tmp_path),
            mode="fill_missing_only",
            flank_size=2,
        )
    )
    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(config),
    )

    diagnostics = _stage_diagnostics(preprocessed)
    assert diagnostics["resolved_site_count"] == 0
    assert diagnostics["unresolved_site_count"] == 3
    assert diagnostics["unresolved_counts_by_reason"] == {
        "accession_not_found": 1,
        "missing_accession": 1,
        "residue_mismatch": 1,
    }

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.api.configs import (
    DatasetPreprocessingConfig,
    DatasetSiteSequenceResolutionConfig,
)
from phospy.datasets.builders.preprocessing import (
    DatasetPreprocessor,
    build_dataset_processing_state,
)
from phospy.datasets.preprocessing.models import PreprocessingPlan
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.intensity_scale_state import (
    intensity_scale_state_from_payload,
)


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
            "localisation_confidence": [0.95, 0.9],
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


def _run_site_sequence_resolution(
    *,
    tmp_path: Path,
    site_sequences: list[object],
    mode: str,
    conflict_policy: str | None = None,
):
    config = DatasetPreprocessingConfig(
        site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
            fasta_path=_write_fasta(tmp_path),
            mode=mode,
            conflict_policy=conflict_policy,
            flank_size=2,
        )
    )
    return DatasetPreprocessor().run(
        phospho=_phospho(),
        site_metadata=_site_metadata(site_sequences=site_sequences),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(config),
    )


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


@pytest.mark.parametrize(
    "mode",
    [
        "fill_missing_only",
        "validate_existing_only",
        "validate_existing_and_fill_missing",
    ],
)
def test_existing_site_sequence_whitespace_is_preserved_exactly_across_non_replace_modes(
    tmp_path: Path,
    mode: str,
) -> None:
    preprocessed = _run_site_sequence_resolution(
        tmp_path=tmp_path,
        site_sequences=["  AASAA  ", "  CCTCC  "],
        mode=mode,
    )

    assert preprocessed.site_metadata.loc["MAPK14;S5;", "site_sequence"] == "  AASAA  "
    assert preprocessed.site_metadata.loc["GSK3B;T6;", "site_sequence"] == "  CCTCC  "
    diagnostics = _stage_diagnostics(preprocessed)
    assert diagnostics["filled_missing_count"] == 0
    assert diagnostics["replaced_existing_count"] == 0


def test_validate_existing_and_fill_missing_preserves_existing_whitespace_and_fills_only_missing(
    tmp_path: Path,
) -> None:
    preprocessed = _run_site_sequence_resolution(
        tmp_path=tmp_path,
        site_sequences=["  AASAA  ", pd.NA],
        mode="validate_existing_and_fill_missing",
    )

    assert preprocessed.site_metadata.loc["MAPK14;S5;", "site_sequence"] == "  AASAA  "
    assert preprocessed.site_metadata.loc["GSK3B;T6;", "site_sequence"] == "CCTCC"
    diagnostics = _stage_diagnostics(preprocessed)
    assert diagnostics["filled_missing_count"] == 1
    assert diagnostics["preserved_existing_count"] == 1
    assert diagnostics["replaced_existing_count"] == 0


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
    assert diagnostics["conflict_policy"] == "preserve_existing"
    assert diagnostics["existing_sequence_conflict_count"] == 1
    assert diagnostics["replaced_existing_count"] == 0
    assert diagnostics["row_diagnostics"][0]["action"] == "preserve_existing"


def test_preserve_existing_conflict_policy_is_durable_in_processing_state(
    tmp_path: Path,
) -> None:
    config = DatasetPreprocessingConfig(
        site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
            fasta_path=_write_fasta(tmp_path),
            mode="validate_existing_and_fill_missing",
            conflict_policy="preserve_existing",
            flank_size=2,
        )
    )
    plan = PreprocessingPlan.from_config(config)
    preprocessed = DatasetPreprocessor().run(
        phospho=_phospho().iloc[:1, :].copy(deep=True),
        site_metadata=_site_metadata(site_sequences=["XXXXX", "CCTCC"])
        .iloc[:1, :]
        .copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=plan,
    )

    processing_state = build_dataset_processing_state(
        plan=plan,
        intensity_scale_state=intensity_scale_state_from_payload(
            {
                "phospho": {
                    "kind": "linear",
                    "transformed": False,
                    "established_by": "bundle.fixture",
                },
                "total": None,
                "quantity": "phosphosite_abundance",
            }
        ),
        preprocessing_trace=preprocessed.preprocessing_trace,
        final_phospho=preprocessed.phospho,
        final_site_metadata=preprocessed.site_metadata,
        final_sample_metadata=preprocessed.sample_metadata,
    )

    resolution = processing_state.site_sequence_resolution
    assert resolution.conflict_policy == "preserve_existing"
    assert resolution.existing_sequence_conflict_count == 1
    assert len(resolution.row_diagnostics) == 1
    assert resolution.row_diagnostics[0].action == "preserve_existing"
    assert resolution.row_diagnostics[0].fasta_site_sequence == "AASAA"


def test_replace_existing_mode_replaces_conflicting_sequence(tmp_path: Path) -> None:
    preprocessed = _run_site_sequence_resolution(
        tmp_path=tmp_path,
        site_sequences=["XXXXX", "CCTCC"],
        mode="replace_existing",
        conflict_policy="replace_existing",
    )

    assert preprocessed.site_metadata.loc["MAPK14;S5;", "site_sequence"] == "AASAA"
    diagnostics = _stage_diagnostics(preprocessed)
    assert diagnostics["conflict_policy"] == "replace_existing"
    assert diagnostics["existing_sequence_conflict_count"] == 1
    assert diagnostics["replaced_existing_count"] == 1
    assert diagnostics["row_diagnostics"][0]["action"] == "replace_existing"


def test_error_conflict_policy_raises_with_structured_row_diagnostics(
    tmp_path: Path,
) -> None:
    config = DatasetPreprocessingConfig(
        site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
            fasta_path=_write_fasta(tmp_path),
            mode="validate_existing_and_fill_missing",
            conflict_policy="error",
            flank_size=2,
        )
    )

    with pytest.raises(
        PhosPyInputError,
        match="conflict_policy='error'",
    ) as caught:
        DatasetPreprocessor().run(
            phospho=_phospho().iloc[:1, :].copy(deep=True),
            site_metadata=_site_metadata(site_sequences=["XXXXX", "CCTCC"])
            .iloc[:1, :]
            .copy(deep=True),
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan.from_config(config),
        )

    diagnostics = caught.value.diagnostics
    assert isinstance(diagnostics, dict)
    row_diagnostics = diagnostics.get("row_diagnostics")
    assert isinstance(row_diagnostics, list)
    assert len(row_diagnostics) == 1
    assert row_diagnostics[0]["row_index"] == 0
    assert row_diagnostics[0]["site_id"] == "MAPK14;S5;"
    assert row_diagnostics[0]["existing_site_sequence"] == "XXXXX"
    assert row_diagnostics[0]["fasta_site_sequence"] == "AASAA"
    assert row_diagnostics[0]["action"] == "error"
    assert row_diagnostics[0]["conflict_policy"] == "error"
    assert row_diagnostics[0]["resolver_version"] == "phospy.sequences.resolver.v1"
    assert row_diagnostics[0]["fasta_source_path"] == diagnostics.get(
        "fasta_source_path"
    )
    assert row_diagnostics[0]["fasta_sha256"] == diagnostics.get("fasta_sha256")


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
            "localisation_confidence": [0.95, 0.9, 0.92],
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


def test_site_sequence_resolution_processing_state_populates_fasta_provenance(
    tmp_path: Path,
) -> None:
    config = DatasetPreprocessingConfig(
        site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
            fasta_path=_write_fasta(tmp_path),
            mode="replace_existing",
            conflict_policy="replace_existing",
            flank_size=2,
        )
    )
    plan = PreprocessingPlan.from_config(config)
    preprocessed = DatasetPreprocessor().run(
        phospho=_phospho(),
        site_metadata=_site_metadata(site_sequences=[pd.NA, "XXXXX"]),
        sample_metadata=None,
        total=None,
        plan=plan,
    )

    processing_state = build_dataset_processing_state(
        plan=plan,
        intensity_scale_state=intensity_scale_state_from_payload(
            {
                "phospho": {
                    "kind": "linear",
                    "transformed": False,
                    "established_by": "bundle.fixture",
                },
                "total": None,
                "quantity": "phosphosite_abundance",
            }
        ),
        preprocessing_trace=preprocessed.preprocessing_trace,
        final_phospho=preprocessed.phospho,
        final_site_metadata=preprocessed.site_metadata,
        final_sample_metadata=preprocessed.sample_metadata,
    )
    diagnostics = _stage_diagnostics(preprocessed)
    resolution = processing_state.site_sequence_resolution

    assert resolution.configured is True
    assert resolution.mode == "replace_existing"
    assert resolution.flank_size == 2
    assert resolution.fasta_source_path == diagnostics.get("fasta_source_path")
    assert resolution.fasta_source_label == diagnostics.get("fasta_source_label")
    assert resolution.fasta_sha256 == diagnostics.get("fasta_sha256")
    assert resolution.resolver_version == diagnostics.get("resolver_version")
    assert resolution.resolved_site_count == diagnostics.get("resolved_site_count")
    assert resolution.unresolved_site_count == diagnostics.get("unresolved_site_count")
    assert resolution.unresolved_counts_by_reason == diagnostics.get(
        "unresolved_counts_by_reason"
    )
    assert resolution.filled_missing_count == diagnostics.get("filled_missing_count")
    assert resolution.replaced_existing_count == diagnostics.get(
        "replaced_existing_count"
    )
    assert resolution.preserved_existing_count == diagnostics.get(
        "preserved_existing_count"
    )
    assert resolution.existing_sequence_conflict_count == diagnostics.get(
        "existing_sequence_conflict_count"
    )
    assert resolution.conflict_policy == "replace_existing"
    assert len(resolution.row_diagnostics) == 2
    assert resolution.row_diagnostics[1].action == "replace_existing"

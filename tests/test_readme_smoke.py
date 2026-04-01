from __future__ import annotations

from pathlib import Path

import pandas as pd

from phospy import KinaseActivityAnalyzer, PhosphoDataset, PhosRPipeline

EXAMPLE_OUTPUT_FILES = {
    "df_phospho_corrected.csv",
    "df_phospho_filtered.csv",
    "df_total_filtered.csv",
    "df_total_unique.csv",
    "kinase_activity_matrix.csv",
    "kinase_target_counts.csv",
    "kinase_target_table.csv",
    "ksea_counts.csv",
    "ksea_scores.csv",
    "mat_phospho_corrected.csv",
    "phosr_input.csv",
    "site_sequences.csv",
}


def test_readme_example_analyzer_runs_end_to_end(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parents[1]

    dataset = PhosphoDataset.from_files(
        repo_root / "examples" / "data" / "total.tsv",
        repo_root / "examples" / "data" / "phospho.tsv",
        phospho_encoding="utf-16le",
    )
    core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

    analyzer = KinaseActivityAnalyzer()
    result = analyzer.load_and_analyze(
        pred_mat_path=repo_root / "examples" / "data" / "predMat.csv",
        phospho_matrix=core.site_matrix.matrix,
        threshold=0.6,
        min_substrates=1,
        top_n_substrates=1,
    )
    analyzer.write_outputs(result, outdir=tmp_path)

    assert result.target_counts.to_dict() == {"PRKACA": 3, "BTK": 2}
    assert EXAMPLE_OUTPUT_FILES - {
        "df_phospho_corrected.csv",
        "df_phospho_filtered.csv",
        "df_total_filtered.csv",
        "df_total_unique.csv",
        "mat_phospho_corrected.csv",
        "phosr_input.csv",
        "site_sequences.csv",
    } <= {path.name for path in tmp_path.iterdir()}


def test_readme_example_pipeline_runs_end_to_end(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    outdir = tmp_path / "output"

    pipeline = PhosRPipeline.from_files(
        total_path=repo_root / "examples" / "data" / "total.tsv",
        phospho_path=repo_root / "examples" / "data" / "phospho.tsv",
        pred_mat_path=repo_root / "examples" / "data" / "predMat.csv",
        phospho_encoding="utf-16le",
        max_unmatched_fraction=0.1,
    )
    outputs = pipeline.run(outdir=outdir)

    assert outputs.core.site_matrix.matrix.index.tolist() == ["BTK;Y551;"]
    assert outputs.kinase_activity is not None
    assert EXAMPLE_OUTPUT_FILES.issubset({path.name for path in outdir.iterdir()})

    expected_matrix = pd.read_csv(
        repo_root / "examples" / "output" / "mat_phospho_corrected.csv",
        index_col=0,
    )
    pd.testing.assert_frame_equal(
        outputs.core.site_matrix.matrix,
        expected_matrix,
        check_index_type=False,
        check_column_type=False,
    )

    expected_target_counts = pd.read_csv(
        repo_root / "examples" / "output" / "kinase_target_counts.csv"
    )
    actual_target_counts = outputs.kinase_activity.target_counts.rename_axis(
        "kinase"
    ).reset_index(name="n_targets")
    pd.testing.assert_frame_equal(actual_target_counts, expected_target_counts)

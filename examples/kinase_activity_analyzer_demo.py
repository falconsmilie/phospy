#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from phospy.activities import KinaseActivityAnalyzer
from phospy.datasets import PhosphoDataset


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    dataset = PhosphoDataset.from_files(
        repo_root / "examples" / "data" / "total.tsv",
        repo_root / "examples" / "data" / "phospho.tsv",
        phospho_encoding="utf-16le",
    )
    core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

    analyzer = KinaseActivityAnalyzer()
    result = analyzer.run(
        pred_mat=analyzer.load_pred_mat(
            repo_root / "examples" / "data" / "predMat.csv"
        ),
        phospho_matrix=core.site_matrix.matrix,
        threshold=0.6,
        min_substrates=1,
        top_n_substrates=1,
    )

    with TemporaryDirectory(prefix="phospy-kinase-activity-") as tmp_dir:
        outdir = Path(tmp_dir)
        analyzer.write_outputs(result, outdir=outdir)
        print(f"Wrote kinase activity outputs to {outdir}")
        print("Target counts")
        print(result.target_counts)
        print()
        print("KSEA scores")
        print(result.ksea_scores.round(4))


if __name__ == "__main__":
    main()

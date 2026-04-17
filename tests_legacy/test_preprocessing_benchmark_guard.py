from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_preprocessing_copy_churn_benchmark_guard_passes() -> None:
    script = REPO_ROOT / "benchmarks" / "measure_preprocessing_copy_churn.py"
    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            "--repeats",
            "1",
            "--warmups",
            "0",
            "--small-genes",
            "150",
            "--small-sites-per-gene",
            "2",
            "--large-genes",
            "500",
            "--large-sites-per-gene",
            "2",
            "--check",
            "--stdout-only",
        ],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr or completed.stdout

    report = json.loads(completed.stdout)
    assert report["guards"]["all_passed"] is True
    assert "public_safe_path" in report["benchmarks"]
    assert "owned_fast_path" in report["benchmarks"]
    assert "large_matrix_owned_path" in report["benchmarks"]

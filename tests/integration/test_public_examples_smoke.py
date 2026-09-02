from __future__ import annotations

import runpy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_DIR = ROOT / "examples"

pytestmark = pytest.mark.integration


@pytest.mark.parametrize(
    ("script_name", "expected_markers"),
    [
        (
            "dataset_builder_demo.py",
            [
                "Supported dataset-builder workflow",
                "Bundled reference-compatible organism: rat",
                "protein_group_id present for all sites: True",
            ],
        ),
        (
            "kinase_workflow_demo.py",
            [
                "Supported kinase workflow",
                "Reference input: ReferencePreset.AUTO",
                "Resolved reference organism: rat",
            ],
        ),
        (
            "signalome_workflow_demo.py",
            [
                "Supported signalome workflow",
                "protein_group_id present for all sites: True",
                "Resolved reference organism: rat",
            ],
        ),
        (
            "reference_bundle_builder_demo.py",
            [
                "Reference bundle builder demo",
                "Organism: mouse",
                "Deterministic fallback reference version:",
            ],
        ),
        (
            "maxquant_importer_demo.py",
            [
                "MaxQuant importer candidate matrix:",
                "MaxQuant filtering diagnostics:",
            ],
        ),
        (
            "enrichment_provenance_demo.py",
            [
                "Enrichment provenance demo",
                "Statistics unchanged: True",
                "Run-provenance source type: phospy_derived_quantitative",
            ],
        ),
        (
            "differential_protein_aware_demo.py",
            [
                "Experimental protein-aware differential workflow",
                "Method: protein_covariate_adjusted_moderated_linear_model_v1",
                "Claim: experimental",
                "Withheld sites: 1",
                "Fallback policy: no_fallback_to_ordinary_differential_lane",
            ],
        ),
    ],
)
def test_public_example_script_runs(
    script_name: str,
    expected_markers: list[str],
    capsys: pytest.CaptureFixture[str],
) -> None:
    runpy.run_path(str(EXAMPLES_DIR / script_name), run_name="__main__")
    captured = capsys.readouterr()

    for marker in expected_markers:
        assert marker in captured.out

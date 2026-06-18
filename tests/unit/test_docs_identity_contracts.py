from __future__ import annotations

from pathlib import Path

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)

ROOT = Path(__file__).resolve().parents[2]
DOCS_ROOT = ROOT / "docs"


def _public_markdown_paths() -> tuple[Path, ...]:
    docs = [
        path
        for path in DOCS_ROOT.rglob("*.md")
        if not _is_under(path, DOCS_ROOT / "testing")
    ]
    return (ROOT / "README.md", *sorted(docs))


def _public_example_paths() -> tuple[Path, ...]:
    docs = [
        DOCS_ROOT / "quickstart.md",
        DOCS_ROOT / "workflow_contracts.md",
        DOCS_ROOT / "output_bundles.md",
        *sorted((DOCS_ROOT / "api").rglob("*.md")),
    ]
    examples = sorted((ROOT / "examples").rglob("*.py"))
    return (ROOT / "README.md", *docs, *examples)


def _is_under(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _combined_text(paths: tuple[Path, ...]) -> str:
    return "\n".join(path.read_text(encoding="utf-8") for path in paths)


def _normalise_whitespace(text: str) -> str:
    return " ".join(text.split())


def test_docs_do_not_describe_display_id_as_analysis_ready_row_identity() -> None:
    combined = _combined_text(_public_markdown_paths()).lower()

    forbidden_fragments = (
        "`display_id` is the analysis-ready row identity",
        "`display_id` is row identity",
        "`display_id` values are the analysis-ready row identity",
        "display labels are analysis-ready row identity",
        "display-site identity is the analysis-ready row identity",
        "direct analysis-ready construction may use display-indexed rows",
    )

    for fragment in forbidden_fragments:
        assert fragment not in combined


def test_docs_state_duplicate_display_ids_and_duplicate_site_key_policy() -> None:
    text = _normalise_whitespace(_combined_text(_public_markdown_paths()))

    assert "Duplicate `display_id` values remain valid" in text
    assert "when the corresponding `site_key` values differ" in text
    assert "Duplicate rows that resolve to the same `site_key`" in text
    assert "scientific ambiguity" in text
    assert "fail by default" in text or "fails by default" in text
    assert "non-error duplicate-site" in text


def test_docs_state_public_differential_result_identity_contract() -> None:
    text = _normalise_whitespace(
        (DOCS_ROOT / "api" / "differential-analysis.md").read_text(encoding="utf-8")
    )

    assert "`DifferentialAnalysisWorkflow.run(...)` returns" in text
    assert "Each contrast result table is indexed by the input `site_key`" in text
    assert "The minimum public identity columns are" in text
    for column_name in (
        "`site_key`",
        "`display_id`",
        "`organism`",
        "`protein_namespace`",
        "`protein_identifier`",
        "`gene_symbol`",
        "`site`",
    ):
        assert column_name in text
    assert "stat-only computation payload" in text
    assert "not a public scientific result object" in text
    assert "not valid `DifferentialAnalysisResult` tables" in text


def test_docs_state_clean_source_archive_hygiene() -> None:
    text = _normalise_whitespace(
        (DOCS_ROOT / "maintenance.md").read_text(encoding="utf-8")
    )

    assert "Source and Release Archive Hygiene" in text
    assert "Do not include generated artefacts" in text
    for generated_path in (
        "`build/`",
        "`dist/`",
        "`site/`",
        "`__pycache__/`",
        "`app-src.zip`",
    ):
        assert generated_path in text


def test_user_docs_do_not_encourage_private_escape_hatches() -> None:
    combined = _combined_text(_public_example_paths())

    forbidden_fragments = ("_from_owned(", "._borrow", "_borrow_", "._site_metadata")
    for fragment in forbidden_fragments:
        assert fragment not in combined


def test_documented_dataset_builder_identity_output_contract() -> None:
    dataset = _duplicate_display_dataset()

    assert dataset.phospho.index.name == "site_key"
    assert dataset.site_metadata.index.name == "site_key"
    assert dataset.site_metadata.loc[:, "site_key"].tolist() == (
        dataset.phospho.index.tolist()
    )
    assert dataset.site_metadata.loc[:, "site_key"].is_unique
    assert dataset.site_metadata.loc[:, "display_id"].tolist() == [
        "MAPK14;Y182;",
        "MAPK14;Y182;",
    ]


def test_documented_differential_workflow_result_identity_columns() -> None:
    dataset = _duplicate_display_dataset(log2=True)
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="control_1",
                condition="control",
                biological_replicate_id="control_1",
            ),
            SampleDesignRecord(
                sample_id="control_2",
                condition="control",
                biological_replicate_id="control_2",
            ),
            SampleDesignRecord(
                sample_id="treated_1",
                condition="treated",
                biological_replicate_id="treated_1",
            ),
            SampleDesignRecord(
                sample_id="treated_2",
                condition="treated",
                biological_replicate_id="treated_2",
            ),
        )
    )
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=design,
            contrasts=(
                Contrast(
                    name="treated_vs_control",
                    numerator_condition="treated",
                    denominator_condition="control",
                ),
            ),
        )
    )

    table = result.table_for("treated_vs_control")

    assert table.index.name == "site_key"
    assert table.loc[:, "site_key"].tolist() == table.index.tolist()
    assert table.loc[:, "site_key"].is_unique
    assert table.loc[:, "display_id"].tolist() == [
        "MAPK14;Y182;",
        "MAPK14;Y182;",
    ]
    for column_name in (
        "site_key",
        "display_id",
        "organism",
        "protein_namespace",
        "protein_identifier",
        "gene_symbol",
        "site",
    ):
        assert column_name in table.columns


def _duplicate_display_dataset(*, log2: bool = False):
    phospho = pd.DataFrame(
        {
            "control_1": [100.0, 120.0],
            "control_2": [110.0, 125.0],
            "treated_1": [180.0, 121.0],
            "treated_2": [190.0, 126.0],
        },
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_A"],
            "protein_accession": ["P28482-1", "P28482-2"],
            "localisation_confidence": [0.95, 0.95],
        },
        index=phospho.index.copy(),
    )
    preprocessing = DatasetPreprocessingConfig(
        site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata"),
    )
    if log2:
        preprocessing = DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
            site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata"),
        )

    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale=None if log2 else "linear",
            preprocessing_config=preprocessing,
        )
    )

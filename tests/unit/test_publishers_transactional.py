from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from phospy.api import Organism, ReferenceBundle
from phospy.contracts.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
    SignalomeWorkflowResult,
)
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared import transactions
from phospy.io.publishers import workflows as publishers
from phospy.io.publishers.workflows import (
    publish_dataset,
    publish_kinase_workflow,
    publish_signalome_workflow,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns


def test_publish_dataset_new_destination_succeeds_and_existing_default_rejects(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = _dataset(label="original")
    output_root = tmp_path / "published"
    written = publish_dataset(dataset, output_root, output_format="csv")

    assert written["dataset.manifest"] == output_root / "dataset" / "manifest.json"
    assert all(output_root in path.parents for path in written.values())
    before = _tree_bytes(output_root)
    write_called = False

    def fail_if_called(*args: Any, **kwargs: Any) -> None:
        nonlocal write_called
        write_called = True
        raise AssertionError("publisher wrote after destination preflight failed")

    monkeypatch.setattr(publishers, "write_table", fail_if_called)

    with pytest.raises(PhosPyInputError) as exc_info:
        publish_dataset(_dataset(label="replacement"), output_root)

    assert str(output_root) in str(exc_info.value)
    assert "Pass overwrite=True" in str(exc_info.value)
    assert write_called is False
    assert _tree_bytes(output_root) == before
    _assert_no_transaction_leftovers(tmp_path, "published")


def test_publish_dataset_existing_file_requires_deliberate_overwrite(
    tmp_path: Path,
) -> None:
    dataset = _dataset(label="replacement")
    output_root = tmp_path / "published"
    output_root.write_bytes(b"old-file")

    with pytest.raises(PhosPyInputError) as exc_info:
        publish_dataset(dataset, output_root)

    assert str(output_root) in str(exc_info.value)
    assert output_root.read_bytes() == b"old-file"

    written = publish_dataset(dataset, output_root, overwrite=True)

    assert output_root.is_dir()
    assert written["dataset.phospho"] == output_root / "dataset" / "phospho.csv"
    assert (output_root / "dataset" / "manifest.json").is_file()
    _assert_no_transaction_leftovers(tmp_path, "published")


def test_publish_dataset_symlink_overwrite_moves_link_not_target_when_supported(
    tmp_path: Path,
) -> None:
    dataset = _dataset(label="replacement")
    target = tmp_path / "external-target"
    target.mkdir()
    (target / "sentinel.txt").write_text("do not touch", encoding="utf-8")
    output_root = tmp_path / "published"
    try:
        output_root.symlink_to(target, target_is_directory=True)
    except OSError as exc:
        pytest.skip(f"symlink creation is not available in this environment: {exc}")

    publish_dataset(dataset, output_root, overwrite=True)

    assert output_root.is_dir()
    assert not output_root.is_symlink()
    assert (target / "sentinel.txt").read_text(encoding="utf-8") == "do not touch"
    _assert_no_transaction_leftovers(tmp_path, "published")


def test_publish_dataset_successful_overwrite_replaces_complete_tree(
    tmp_path: Path,
) -> None:
    output_root = tmp_path / "published"
    original_dataset = _dataset(label="original")
    replacement_dataset = _dataset(label="replacement")
    publish_dataset(original_dataset, output_root)
    stale_path = output_root / "dataset" / "stale.txt"
    stale_path.write_text("stale", encoding="utf-8")
    expected_root = tmp_path / "expected"
    publish_dataset(replacement_dataset, expected_root)

    written = publish_dataset(replacement_dataset, output_root, overwrite=True)

    assert _tree_bytes(output_root) == _tree_bytes(expected_root)
    assert not stale_path.exists()
    assert all(output_root in path.parents for path in written.values())
    assert not any(
        ".tmp-" in path.name or ".previous-" in path.name for path in written.values()
    )
    _assert_no_transaction_leftovers(tmp_path, "published")


def test_publish_dataset_table_write_failure_preserves_existing_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "published"
    publish_dataset(_dataset(label="original"), output_root)
    before = _tree_bytes(output_root)
    original_write_table = publishers.write_table

    def fail_phospho_write(table: pd.DataFrame, path: Path) -> None:
        if path.name == "phospho.csv":
            raise PhosPyInputError("simulated table-write failure")
        original_write_table(table, path)

    monkeypatch.setattr(publishers, "write_table", fail_phospho_write)

    with pytest.raises(PhosPyInputError, match="simulated table-write failure"):
        publish_dataset(_dataset(label="replacement"), output_root, overwrite=True)

    assert _tree_bytes(output_root) == before
    _assert_no_transaction_leftovers(tmp_path, "published")


def test_publish_dataset_manifest_write_failure_preserves_existing_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "published"
    publish_dataset(_dataset(label="original"), output_root)
    before = _tree_bytes(output_root)

    def fail_manifest(path: Path, payload: dict[str, object]) -> None:
        del path, payload
        raise PhosPyInputError("simulated manifest-write failure")

    monkeypatch.setattr(publishers, "_write_manifest", fail_manifest)

    with pytest.raises(PhosPyInputError, match="simulated manifest-write failure"):
        publish_dataset(_dataset(label="replacement"), output_root, overwrite=True)

    assert _tree_bytes(output_root) == before
    _assert_no_transaction_leftovers(tmp_path, "published")


def test_publish_dataset_first_time_failure_leaves_no_final_partial_root(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "published"

    def fail_manifest(path: Path, payload: dict[str, object]) -> None:
        del path, payload
        raise PhosPyInputError("simulated first-time manifest failure")

    monkeypatch.setattr(publishers, "_write_manifest", fail_manifest)

    with pytest.raises(
        PhosPyInputError,
        match="simulated first-time manifest failure",
    ):
        publish_dataset(_dataset(label="new"), output_root)

    assert not output_root.exists()
    _assert_no_transaction_leftovers(tmp_path, "published")


def test_publish_kinase_failure_after_dataset_staging_does_not_publish_dataset_only_tree(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "published"

    def fail_kinase_scoring(
        result: KinaseWorkflowResult,
        workflow_dir: Path,
        suffix: str,
        written: dict[str, Path],
    ) -> None:
        del result, workflow_dir, suffix, written
        raise PhosPyInputError("simulated kinase-stage failure")

    monkeypatch.setattr(
        publishers, "_write_kinase_scoring_outputs", fail_kinase_scoring
    )

    with pytest.raises(PhosPyInputError, match="simulated kinase-stage failure"):
        publish_kinase_workflow(_kinase_result(label="new"), output_root)

    assert not output_root.exists()
    _assert_no_transaction_leftovers(tmp_path, "published")


def test_publish_kinase_failure_after_dataset_staging_preserves_existing_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "published"
    publish_kinase_workflow(_kinase_result(label="original"), output_root)
    before = _tree_bytes(output_root)

    def fail_kinase_scoring(
        result: KinaseWorkflowResult,
        workflow_dir: Path,
        suffix: str,
        written: dict[str, Path],
    ) -> None:
        del result, workflow_dir, suffix, written
        raise PhosPyInputError("simulated kinase-stage failure")

    monkeypatch.setattr(
        publishers, "_write_kinase_scoring_outputs", fail_kinase_scoring
    )

    with pytest.raises(PhosPyInputError, match="simulated kinase-stage failure"):
        publish_kinase_workflow(
            _kinase_result(label="replacement"),
            output_root,
            overwrite=True,
        )

    assert _tree_bytes(output_root) == before
    _assert_no_transaction_leftovers(tmp_path, "published")


@pytest.mark.parametrize(
    "failure_hook",
    ["kinase", "signalome"],
)
def test_publish_signalome_failure_after_staging_preserves_existing_complete_tree(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure_hook: str,
) -> None:
    output_root = tmp_path / "published"
    publish_signalome_workflow(_signalome_result(label="original"), output_root)
    before = _tree_bytes(output_root)

    if failure_hook == "kinase":

        def fail_kinase_scoring(
            result: KinaseWorkflowResult,
            workflow_dir: Path,
            suffix: str,
            written: dict[str, Path],
        ) -> None:
            del result, workflow_dir, suffix, written
            raise PhosPyInputError("simulated signalome kinase-stage failure")

        monkeypatch.setattr(
            publishers,
            "_write_kinase_scoring_outputs",
            fail_kinase_scoring,
        )
    else:

        def fail_signalome_outputs(
            result: SignalomeWorkflowResult,
            workflow_dir: Path,
            suffix: str,
            written: dict[str, Path],
        ) -> None:
            del result, workflow_dir, suffix, written
            raise PhosPyInputError("simulated signalome-stage failure")

        monkeypatch.setattr(
            publishers,
            "_write_signalome_outputs",
            fail_signalome_outputs,
        )

    with pytest.raises(PhosPyInputError, match="simulated signalome.*failure"):
        publish_signalome_workflow(
            _signalome_result(label="replacement"),
            output_root,
            overwrite=True,
        )

    assert _tree_bytes(output_root) == before
    _assert_no_transaction_leftovers(tmp_path, "published")


def test_publish_signalome_first_time_failure_after_kinase_staging_leaves_no_mixed_tree(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "published"

    def fail_signalome_outputs(
        result: SignalomeWorkflowResult,
        workflow_dir: Path,
        suffix: str,
        written: dict[str, Path],
    ) -> None:
        del result, workflow_dir, suffix, written
        raise PhosPyInputError("simulated signalome-stage failure")

    monkeypatch.setattr(publishers, "_write_signalome_outputs", fail_signalome_outputs)

    with pytest.raises(PhosPyInputError, match="simulated signalome-stage failure"):
        publish_signalome_workflow(_signalome_result(label="new"), output_root)

    assert not output_root.exists()
    _assert_no_transaction_leftovers(tmp_path, "published")


def test_publish_dataset_promotion_failure_rolls_back_existing_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "published"
    publish_dataset(_dataset(label="original"), output_root)
    before = _tree_bytes(output_root)
    _fail_path_replace_when(
        monkeypatch,
        lambda source, target: (
            source.name.startswith(".published.tmp-") and target == output_root
        ),
        PermissionError("simulated promotion failure"),
    )

    with pytest.raises(PhosPyInputError, match="rollback restored"):
        publish_dataset(_dataset(label="replacement"), output_root, overwrite=True)

    assert _tree_bytes(output_root) == before
    _assert_no_transaction_leftovers(tmp_path, "published")


def test_publish_dataset_backup_cleanup_failure_reports_recovery_copy(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "published"
    publish_dataset(_dataset(label="original"), output_root)
    original_tree = _tree_bytes(output_root)
    expected_root = tmp_path / "expected"
    publish_dataset(_dataset(label="replacement"), expected_root)
    expected_tree = _tree_bytes(expected_root)
    _fail_rmtree_when(
        monkeypatch,
        lambda path: path.name.startswith(".published.previous-"),
        PermissionError("simulated backup cleanup failure"),
    )

    with pytest.raises(
        PhosPyInputError,
        match="promoted successfully.*backup cleanup failed",
    ) as exc_info:
        publish_dataset(_dataset(label="replacement"), output_root, overwrite=True)

    backup_root = _single_recovery_backup(tmp_path, output_name="published")
    assert str(backup_root) in str(exc_info.value)
    assert _tree_bytes(output_root) == expected_tree
    assert _tree_bytes(backup_root) == original_tree
    assert not list(tmp_path.glob(".published.tmp-*"))


def _dataset(*, label: str) -> AnalysisReadyPhosphoDataset:
    site_index = _site_index()
    display_ids = _display_ids()
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [3.0, 4.0],
        },
        index=site_index,
        dtype=float,
    )
    if label == "replacement":
        phospho = phospho + 10.0
    elif label == "new":
        phospho = phospho + 20.0
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": _proteins(),
            "site": _sites(),
            "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in _sites()],
            "protein_id": _proteins(),
        },
        index=site_index,
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False,
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _kinase_result(*, label: str) -> KinaseWorkflowResult:
    dataset = _dataset(label=label)
    score_matrix = pd.DataFrame(
        {"K1": [0.25, 0.75]},
        index=dataset.phospho.index.copy(),
        dtype=float,
    )
    if label == "replacement":
        score_matrix = score_matrix + 0.1
    prediction_matrix = score_matrix.clip(lower=0.0, upper=1.0)
    return KinaseWorkflowResult(
        dataset=dataset,
        references=_reference_bundle(),
        scoring_result=KinaseScoringResult(profile_scores=score_matrix),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


def _signalome_result(*, label: str) -> SignalomeWorkflowResult:
    kinase_result = _kinase_result(label=label)
    site_key = str(kinase_result.dataset.phospho.index[0])
    display_id = _display_ids()[0]
    assignment = pd.DataFrame(
        {
            "site_key": [site_key],
            "display_id": [display_id],
            "gene_symbol": [_proteins()[0]],
            "site": [_sites()[0]],
            "protein_group_id": [_proteins()[0]],
            "protein_accession": [""],
            "isoform_id": [""],
            "module_id": [1],
            "top_kinase": ["K1"],
            "top_score": [0.9],
            "top_kinase_candidates": [("K1",)],
            "top_kinase_weights": [(("K1", 1.0),)],
            "top_kinase_tie_count": [1],
            "top_kinase_is_ambiguous": [False],
            "top_kinase_selection_policy": ["max_score_then_lexicographic_tiebreak"],
            "module_top_kinase": ["K1"],
            "module_top_kinase_candidates": [("K1",)],
            "module_top_kinase_tie_count": [1],
            "module_top_kinase_is_ambiguous": [False],
            "module_top_kinase_selection_policy": [
                "max_score_then_lexicographic_tiebreak"
            ],
        },
        index=pd.Index([site_key], name="site_key"),
    )
    modules = pd.DataFrame(
        {"K1": [100.0]},
        index=pd.Index([1], name="module_id"),
        dtype=float,
    )
    edges = pd.DataFrame(
        columns=[
            "source_kinase",
            "target_kinase",
            "correlation",
            "valid_observations",
        ]
    )
    return SignalomeWorkflowResult(
        dataset=kinase_result.dataset,
        kinase_result=kinase_result,
        module_assignments=SignalomeAssignments(table=assignment),
        signalome_modules=SignalomeModules(table=modules),
        kinase_network=KinaseNetwork(edges=edges),
        expanded_signalome=None,
    )


def _reference_bundle() -> ReferenceBundle:
    display_ids = _display_ids()
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K1"],
                "substrate_site": display_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in _sites()]},
            index=pd.Index(display_ids, name="site_id"),
        ),
    )


def _site_index() -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=_proteins(),
        sites=_sites(),
        organism="rat",
        protein_namespace="protein_id",
    )


def _display_ids() -> list[str]:
    return [
        f"{protein};{site};"
        for protein, site in zip(_proteins(), _sites(), strict=True)
    ]


def _proteins() -> list[str]:
    return ["MAPK14", "GSK3B"]


def _sites() -> list[str]:
    return ["Y182", "S9"]


def _tree_bytes(root: Path) -> dict[str, bytes]:
    return {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _assert_no_transaction_leftovers(parent: Path, output_name: str) -> None:
    assert not list(parent.glob(f".{output_name}.tmp-*"))
    assert not list(parent.glob(f".{output_name}.previous-*"))


def _single_recovery_backup(tmp_path: Path, *, output_name: str) -> Path:
    backup_roots = list(tmp_path.glob(f".{output_name}.previous-*"))
    assert len(backup_roots) == 1
    return backup_roots[0]


def _fail_path_replace_when(
    monkeypatch: pytest.MonkeyPatch,
    predicate: Callable[[Path, Path], bool],
    error: OSError,
) -> None:
    original_replace = Path.replace

    def replace(self: Path, target: str | Path) -> Path:
        source = Path(self)
        target_path = Path(target)
        if predicate(source, target_path):
            raise error
        return original_replace(self, target)

    monkeypatch.setattr(Path, "replace", replace)


def _fail_rmtree_when(
    monkeypatch: pytest.MonkeyPatch,
    predicate: Callable[[Path], bool],
    error: OSError,
) -> None:
    original_rmtree = transactions.shutil.rmtree

    def rmtree(path: str | Path, *args: Any, **kwargs: Any) -> None:
        root = Path(path)
        if predicate(root):
            raise error
        original_rmtree(path, *args, **kwargs)

    monkeypatch.setattr(transactions.shutil, "rmtree", rmtree)

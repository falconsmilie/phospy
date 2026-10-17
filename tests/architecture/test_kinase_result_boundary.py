from __future__ import annotations

import ast
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pandas.testing as pdt

import phospy as phospy_root
import phospy.advanced as phospy_advanced
import phospy.api as phospy_api
import phospy.kinase.caveats as kinase_caveats
import phospy.kinase.policies as kinase_policies
import phospy.kinase.scoring.rank_fusion as scoring_rank_fusion
import phospy.kinase.scoring.tables as kinase_tables
from phospy._api_inventory import (
    ADVANCED_PUBLIC_API,
    ORDINARY_PUBLIC_API,
    ROOT_PUBLIC_API,
)
from phospy.io.importers.contracts import (
    PhosphositeImporter,
    PhosphositeImportRequest,
)
from phospy.kinase.prediction.tables import KinasePredictionMatrix
from phospy.kinase.requests import KinaseWorkflowRequest
from phospy.kinase.result_boundary import (
    KINASE_DOWNSTREAM_LIMITATION_CAVEAT_CODES,
    KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE,
    KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE,
    KINASE_SCORING_LIMITATION_CAVEAT_CODE,
    SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY,
    DownstreamScoreSelectionPolicy,
    DownstreamScoreSource,
    KinaseScoreMatrix,
)
from phospy.kinase.scoring.results import KinaseScoringResult
from phospy.provenance.scientific_policy_models import ScientificPolicyId
from phospy.signalome.requests import SignalomeWorkflowRequest
from tests.support.site_keys import protein_site_key_index

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"

_BOUNDARY_MODULE = "phospy.kinase.result_boundary"
_SIGNALOME_SOURCE_PATHS = (
    "src/phospy/signalome/caveats.py",
    "src/phospy/signalome/component_models.py",
    "src/phospy/signalome/contracts.py",
    "src/phospy/signalome/interpreter.py",
    "src/phospy/signalome/matrix_alignment.py",
    "src/phospy/signalome/provenance.py",
    "src/phospy/signalome/score_matrix_selection.py",
    "src/phospy/signalome/validator.py",
)
_GOVERNED_IMPORTS = frozenset(
    {
        (
            "src/phospy/signalome/caveats.py",
            "phospy.science.scoring.policy_models",
        ),
        (
            "src/phospy/signalome/caveats.py",
            "phospy.workflows.kinase.caveats",
        ),
        (
            "src/phospy/signalome/component_models.py",
            "phospy.science.prediction.scoring",
        ),
        (
            "src/phospy/signalome/component_models.py",
            "phospy.science.scoring.policy_models",
        ),
        (
            "src/phospy/signalome/contracts.py",
            "phospy.science.prediction.scoring",
        ),
        (
            "src/phospy/signalome/contracts.py",
            "phospy.science.scoring.policy_models",
        ),
        ("src/phospy/signalome/contracts.py", "phospy.kinase.scoring.tables"),
        (
            "src/phospy/signalome/interpreter.py",
            "phospy.science.prediction.internal_view",
        ),
        (
            "src/phospy/signalome/interpreter.py",
            "phospy.science.prediction.scoring",
        ),
        (
            "src/phospy/signalome/matrix_alignment.py",
            "phospy.science.scoring.policy_models",
        ),
        (
            "src/phospy/signalome/provenance.py",
            "phospy.science.scoring.policy_models",
        ),
        (
            "src/phospy/signalome/score_matrix_selection.py",
            "phospy.science.prediction.internal_view",
        ),
        (
            "src/phospy/signalome/score_matrix_selection.py",
            "phospy.science.prediction.scoring",
        ),
        (
            "src/phospy/signalome/score_matrix_selection.py",
            "phospy.science.scoring.policy_models",
        ),
        (
            "src/phospy/signalome/validator.py",
            "phospy.science.prediction.internal_view",
        ),
        (
            "src/phospy/signalome/validator.py",
            "phospy.science.scoring.policy_models",
        ),
    }
)


def _score_frame(offset: float) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "K1": [0.1 + offset, 0.2 + offset],
            "K2": [0.3 + offset, 0.4 + offset],
        },
        index=protein_site_key_index(
            protein_identifiers=["P1", "P2"],
            sites=["S1", "S2"],
        ),
    )


def test_kinase_result_boundary_owns_each_moved_passive_type() -> None:
    assert DownstreamScoreSource.__module__ == _BOUNDARY_MODULE
    assert DownstreamScoreSelectionPolicy.__module__ == _BOUNDARY_MODULE
    assert KinaseScoreMatrix.__module__ == _BOUNDARY_MODULE
    assert KinasePredictionMatrix.__module__ == "phospy.kinase.prediction.tables"

    assert not hasattr(kinase_policies, "DownstreamScoreSource")
    assert not hasattr(scoring_rank_fusion, "DownstreamScoreSelectionPolicy")
    assert not hasattr(
        scoring_rank_fusion,
        "SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY",
    )
    assert not hasattr(kinase_tables, "KinaseScoreMatrix")
    assert not hasattr(kinase_tables, "KinasePredictionMatrix")
    assert not hasattr(
        kinase_caveats,
        "KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE",
    )
    assert not hasattr(
        kinase_caveats,
        "KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE",
    )
    assert not hasattr(kinase_caveats, "KINASE_SCORING_LIMITATION_CAVEAT_CODE")


def test_request_definitions_have_approved_concrete_modules_and_identity() -> None:
    expected = {
        "PhosphositeImportRequest": PhosphositeImportRequest,
        "PhosphositeImporter": PhosphositeImporter,
        "KinaseWorkflowRequest": KinaseWorkflowRequest,
        "SignalomeWorkflowRequest": SignalomeWorkflowRequest,
    }

    assert KinaseWorkflowRequest.__module__ == "phospy.kinase.requests"
    for name, definition in expected.items():
        expected_module = {
            "KinaseWorkflowRequest": "phospy.kinase.requests",
            "SignalomeWorkflowRequest": "phospy.signalome.requests",
        }.get(name, "phospy.io.importers.contracts")
        assert definition.__module__ == expected_module
        assert getattr(phospy_api, name) is definition


def test_policy_singleton_and_caveat_grouping_are_exact() -> None:
    result = KinaseScoringResult(profile_scores=_score_frame(0.0))
    policy = SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY

    assert result.downstream_score_selection_policy is policy
    assert policy.record.to_payload() == {
        "id": ScientificPolicyId.SIGNALOME_DOWNSTREAM_SCORE_SELECTION.value,
        "name": "signalome_downstream_score_rank_weighted_preferred_v1",
        "version": "1",
        "description": (
            "Select PhosR-inspired rank-weighted fusion scores when available; "
            "otherwise use profile-only scores."
        ),
        "parameters": {
            "preferred_source": "rank_weighted_fusion_scores",
            "fallback_source": "profile_scores",
        },
        "assumptions": [
            "Downstream score-lane selection changes which upstream signal is "
            "interpreted by signalome stages."
        ],
        "output_scale": "Selected downstream kinase-site score matrix.",
        "quantitative_meaning": "relative_downstream_support",
    }
    assert isinstance(KINASE_DOWNSTREAM_LIMITATION_CAVEAT_CODES, frozenset)
    assert KINASE_DOWNSTREAM_LIMITATION_CAVEAT_CODES == {
        KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE,
        KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE,
        KINASE_SCORING_LIMITATION_CAVEAT_CODE,
    }
    assert KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE == (
        "kinase_non_default_reference_source"
    )
    assert KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE == (
        "kinase_reference_score_fallback"
    )
    assert KINASE_SCORING_LIMITATION_CAVEAT_CODE == (
        "kinase_non_phosr_equivalent_scoring"
    )


def test_all_four_authoritative_score_sources_return_defensive_exact_lanes() -> None:
    profile = _score_frame(0.0)
    fusion = _score_frame(0.1)
    motif = _score_frame(0.2)
    combined = _score_frame(0.3)

    cases = (
        (DownstreamScoreSource.PROFILE_SCORES, profile, {}),
        (
            DownstreamScoreSource.RANK_WEIGHTED_FUSION_SCORES,
            fusion,
            {"rank_weighted_fusion_scores": fusion},
        ),
        (
            DownstreamScoreSource.KINASE_LIBRARY_MOTIF_SCORES,
            motif,
            {"kinase_library_motif_scores": motif},
        ),
        (
            DownstreamScoreSource.COMBINED_PROFILE_MOTIF_SCORES,
            combined,
            {"combined_profile_motif_scores": combined},
        ),
    )
    for source, expected, optional_scores in cases:
        result = KinaseScoringResult(
            profile_scores=profile,
            score_source=source,
            **optional_scores,
        )
        first = result.authoritative_scores
        second = result.authoritative_scores
        pdt.assert_frame_equal(first, expected)
        assert first is not second
        first.iloc[0, 0] = -1.0
        pdt.assert_frame_equal(second, expected)
        assert result.score_source == source.value


def test_signalome_sources_have_none_of_the_sixteen_governed_imports() -> None:
    observed: set[tuple[str, str]] = set()
    for relative_path in _SIGNALOME_SOURCE_PATHS:
        tree = ast.parse((PROJECT_ROOT / relative_path).read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                observed.update((relative_path, alias.name) for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module is not None:
                observed.add((relative_path, node.module))

    assert observed.isdisjoint(_GOVERNED_IMPORTS)


def test_boundary_and_signalome_imports_do_not_load_kinase_implementation() -> None:
    modules = (_BOUNDARY_MODULE,) + tuple(
        relative_path.removeprefix("src/").removesuffix(".py").replace("/", ".")
        for relative_path in _SIGNALOME_SOURCE_PATHS
    )
    forbidden = (
        "phospy.kinase.prediction.execution",
        "phospy.kinase.prediction.internal_view",
        "phospy.kinase.caveats",
        "phospy.kinase.executor",
        "phospy.kinase.interpreter",
        "phospy.kinase.provenance",
        "phospy.kinase.result_assembly",
        "phospy.kinase.scoring.execution",
        "phospy.kinase.validator",
    )
    code = "\n".join(
        [
            "import importlib, sys",
            f"modules = {modules!r}",
            f"forbidden = {forbidden!r}",
            "for name in modules: importlib.import_module(name)",
            "loaded = sorted(name for name in sys.modules if name in forbidden)",
            "assert loaded == [], loaded",
        ]
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        (str(SRC_ROOT), environment.get("PYTHONPATH", ""))
    )

    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr


def test_public_inventories_remain_ordered_3_66_124() -> None:
    assert tuple(phospy_root.__all__) == tuple(ROOT_PUBLIC_API)
    assert tuple(phospy_api.__all__) == tuple(ORDINARY_PUBLIC_API)
    assert tuple(phospy_advanced.__all__) == tuple(ADVANCED_PUBLIC_API)
    assert tuple(
        map(
            len,
            (phospy_root.__all__, phospy_api.__all__, phospy_advanced.__all__),
        )
    ) == (3, 66, 124)

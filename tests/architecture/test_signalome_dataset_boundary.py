from __future__ import annotations

from pathlib import Path

import pytest

from tools.architecture.signalome_dataset_boundary import (
    APPROVED_DATASET_IMPORTS,
    audit_signalome_sources,
    audit_signalome_tree,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"


def _audit(
    source: str,
    *,
    consumer: str = "phospy.signalome.results",
    local: dict[str, str] | None = None,
):
    return audit_signalome_sources(
        {consumer: source},
        local_phospy_sources=local,
    )


def _categories(report) -> set[str]:
    return {violation.category for violation in report.violations}


@pytest.mark.parametrize(
    ("consumer", "module", "symbol"),
    sorted(APPROVED_DATASET_IMPORTS),
)
def test_each_approved_dataset_import_is_accepted(
    consumer: str, module: str, symbol: str
) -> None:
    report = _audit(f"from {module} import {symbol}\n", consumer=consumer)

    assert report.violations == ()
    assert report.observed_imports == ((consumer, module, symbol),)


def test_real_signalome_tree_is_fully_visited_with_exact_import_multiset() -> None:
    report = audit_signalome_tree(SRC_ROOT)
    expected_modules = {
        "phospy."
        + str(path.relative_to(SRC_ROOT / "phospy"))
        .replace("\\", ".")
        .replace("/", ".")
        .removesuffix(".py")
        .removesuffix(".__init__")
        for path in (SRC_ROOT / "phospy" / "signalome").rglob("*.py")
    }

    assert set(report.visited_modules) == expected_modules
    assert report.violations == ()
    assert set(report.import_counts) == APPROVED_DATASET_IMPORTS
    assert set(report.import_counts.values()) == {1}


@pytest.mark.parametrize(
    ("consumer", "source", "category"),
    (
        (
            "phospy.signalome.executor",
            "from phospy.datasets.internal_view import DatasetInternalView\n",
            "wrong_consumer",
        ),
        (
            "phospy.signalome.results",
            "from phospy.datasets.models import DatasetPreprocessingReport\n",
            "unapproved_dataset_symbol",
        ),
        (
            "phospy.signalome.results",
            "from phospy.datasets.models import DatasetProcessingState\n",
            "unapproved_dataset_symbol",
        ),
        (
            "phospy.signalome.results",
            "from phospy.datasets.models import AnalysisReadyPhosphoDataset as Dataset\n",
            "aliased_import",
        ),
        (
            "phospy.signalome.results",
            "from phospy.datasets.models import *\n",
            "wildcard_import",
        ),
        (
            "phospy.signalome.executor",
            "import phospy.datasets\n",
            "whole_module_import",
        ),
        (
            "phospy.signalome.executor",
            "import phospy.datasets.models as models\n",
            "whole_module_import",
        ),
        (
            "phospy.signalome.executor",
            "from phospy.datasets import models\n",
            "unapproved_dataset_symbol",
        ),
        (
            "phospy.signalome.results",
            "from ..datasets.models import AnalysisReadyPhosphoDataset\n",
            "relative_import",
        ),
        (
            "phospy.signalome.results",
            "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from phospy.datasets.models import AnalysisReadyPhosphoDataset as Dataset\n",
            "aliased_import",
        ),
    ),
)
def test_unapproved_static_import_forms_are_rejected(
    consumer: str, source: str, category: str
) -> None:
    assert category in _categories(_audit(source, consumer=consumer))


def test_approved_type_checking_import_is_counted_once() -> None:
    source = """
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from phospy.datasets.models import AnalysisReadyPhosphoDataset
"""
    report = _audit(source)

    assert report.violations == ()
    assert (
        report.import_counts[
            (
                "phospy.signalome.results",
                "phospy.datasets.models",
                "AnalysisReadyPhosphoDataset",
            )
        ]
        == 1
    )


def test_duplicate_approved_import_is_rejected_across_type_checking_scope() -> None:
    source = """
from typing import TYPE_CHECKING
from phospy.datasets.models import AnalysisReadyPhosphoDataset
if TYPE_CHECKING:
    from phospy.datasets.models import AnalysisReadyPhosphoDataset
"""
    assert "duplicate_import" in _categories(_audit(source))


@pytest.mark.parametrize(
    "source",
    (
        "import phospy\nvalue = phospy.datasets.models.AnalysisReadyPhosphoDataset\n",
        "import phospy.datasets.models as models\nvalue = models.AnalysisReadyPhosphoDataset\n",
        "import phospy.datasets.models as models\nbox = [models]\nvalue = box[0].DatasetInternalView\n",
        "import phospy.datasets.models as models\nbox = {'m': models}\nvalue = box['m']['AnalysisReadyPhosphoDataset']\n",
        "import phospy.datasets.models as models\nvalue = vars(models)['AnalysisReadyPhosphoDataset']\n",
        "import phospy.datasets.models as models\nvalue = models.__dict__['AnalysisReadyPhosphoDataset']\n",
        "import phospy.datasets.models as models\nvalue = getattr(models, 'AnalysisReadyPhosphoDataset')\n",
        "import phospy.datasets.models as models\nvalue = hasattr(models, 'DatasetInternalView')\n",
        "import operator\nimport phospy.datasets.models as models\nvalue = operator.attrgetter('AnalysisReadyPhosphoDataset')(models)\n",
    ),
)
def test_dataset_module_indirection_is_rejected(source: str) -> None:
    categories = _categories(_audit(source, consumer="phospy.signalome.executor"))
    assert categories & {
        "whole_module_import",
        "module_indirection",
        "reflective_acquisition",
    }


@pytest.mark.parametrize(
    "source",
    (
        "import importlib\nmodule = importlib.import_module('phospy.datasets.models')\n",
        "import importlib as il\nmodule = il.import_module('phospy.datasets.models')\n",
        "from importlib import import_module\nmodule = import_module('phospy.' + 'datasets.models')\n",
        "module = __import__('phospy.datasets.models', fromlist=('AnalysisReadyPhosphoDataset',))\n",
        "loader = __import__\nmodule = loader('phospy.datasets.models')\n",
        "from importlib import import_module\nroot = 'phospy.datasets'\nmodule = import_module(root + suffix)\nvalue = getattr(module, 'DatasetInternalView')\n",
    ),
)
def test_dynamic_dataset_acquisition_is_rejected(source: str) -> None:
    assert "dynamic_import" in _categories(
        _audit(source, consumer="phospy.signalome.executor")
    )


@pytest.mark.parametrize(
    ("route", "local"),
    (
        (
            "from phospy.fixture_one import AnalysisReadyPhosphoDataset\n",
            {
                "phospy.fixture_one": (
                    "from phospy.datasets.models import AnalysisReadyPhosphoDataset\n"
                    "__all__ = ['AnalysisReadyPhosphoDataset']\n"
                )
            },
        ),
        (
            "from phospy.fixture_two import AnalysisReadyPhosphoDataset\n",
            {
                "phospy.fixture_one": (
                    "from phospy.datasets.models import AnalysisReadyPhosphoDataset\n"
                    "__all__ = ['AnalysisReadyPhosphoDataset']\n"
                ),
                "phospy.fixture_two": (
                    "from phospy.fixture_one import AnalysisReadyPhosphoDataset\n"
                    "NAMES = ['AnalysisReadyPhosphoDataset']\n"
                    "__all__ = [] + NAMES\n"
                ),
            },
        ),
    ),
)
def test_local_dataset_reexport_shortcuts_are_rejected(
    route: str, local: dict[str, str]
) -> None:
    assert "reexport_import" in _categories(_audit(route, local=local))


@pytest.mark.parametrize(
    "member",
    (
        "_phospho",
        "_site_metadata",
        "_sample_metadata",
        "_total",
        "_comparisons",
        "_internal_frame_store",
        "_internal_frame_store_for_current_frames",
        "_borrow_phospho",
    ),
)
@pytest.mark.parametrize(
    "template",
    (
        "def read(dataset: AnalysisReadyPhosphoDataset):\n    return dataset.{member}\n",
        "def read(dataset: DatasetInternalView):\n    return getattr(dataset, '{member}')\n",
        "def read(dataset: AnalysisReadyPhosphoDataset):\n    return vars(dataset)['{member}']\n",
        "def read(dataset: AnalysisReadyPhosphoDataset):\n    return dataset.__dict__['{member}']\n",
        "def read(dataset: AnalysisReadyPhosphoDataset):\n    return dataset['{member}']\n",
        "def read(request):\n    dataset = request.kinase_result.dataset\n    return dataset.{member}\n",
    ),
)
def test_private_dataset_storage_is_rejected(member: str, template: str) -> None:
    source = template.format(member=member)
    assert "private_storage" in _categories(
        _audit(source, consumer="phospy.signalome.executor")
    )


@pytest.mark.parametrize(
    "source",
    (
        "from .models import SignalomeAssignments\n",
        "value = result.diagnostics\n",
        "value = dataset.processing_state\n",
        "import numpy\nvalue = numpy.linalg.norm([1.0])\n",
        "from importlib import import_module\nmodule = import_module('numpy.linalg')\n",
        "class Local:\n    pass\nvalue = getattr(Local(), '_phospho')\n",
    ),
)
def test_scope_controls_remain_legal(source: str) -> None:
    assert _audit(source, consumer="phospy.signalome.executor").violations == ()

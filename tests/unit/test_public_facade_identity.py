from __future__ import annotations

import importlib
import operator
import sys
from collections.abc import MutableSequence
from typing import cast

import pytest


def test_root_type_visibility_preserves_lazy_runtime_loading_and_identity() -> None:
    loaded_phospy_modules = {
        name: module
        for name, module in sys.modules.items()
        if name == "phospy" or name.startswith("phospy.")
    }
    try:
        for name in loaded_phospy_modules:
            del sys.modules[name]

        phospy = importlib.import_module("phospy")
        assert "phospy.datasets.builder" not in sys.modules
        assert "phospy.kinase.public" not in sys.modules
        assert "phospy.signalome.public" not in sys.modules
        assert phospy.__all__ == [
            "AnalysisReadyDatasetBuilder",
            "KinaseWorkflow",
            "SignalomeWorkflow",
        ]

        from phospy.datasets.builder import AnalysisReadyDatasetBuilder
        from phospy.kinase.public import KinaseWorkflow
        from phospy.signalome.public import SignalomeWorkflow

        assert phospy.AnalysisReadyDatasetBuilder is AnalysisReadyDatasetBuilder
        assert phospy.KinaseWorkflow is KinaseWorkflow
        assert phospy.SignalomeWorkflow is SignalomeWorkflow
    finally:
        for name in tuple(sys.modules):
            if name == "phospy" or name.startswith("phospy."):
                del sys.modules[name]
        sys.modules.update(loaded_phospy_modules)


def test_advanced_inventory_preserves_authoritative_tuple_ownership() -> None:
    from phospy import _api_inventory, advanced

    assert advanced.__all__ is _api_inventory.ADVANCED_PUBLIC_API
    assert isinstance(advanced.__all__, tuple)
    assert len(advanced.__all__) == 124
    with pytest.raises(TypeError):
        operator.setitem(cast(MutableSequence[object], advanced.__all__), 0, "changed")

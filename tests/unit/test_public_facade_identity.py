from __future__ import annotations

import operator
from collections.abc import MutableSequence
from typing import cast

import pytest


def test_root_exposes_only_version_metadata() -> None:
    import phospy

    assert phospy.__all__ == []
    assert phospy.__version__
    assert "__getattr__" not in vars(phospy)
    for name in ("AnalysisReadyDatasetBuilder", "KinaseWorkflow", "SignalomeWorkflow"):
        assert not hasattr(phospy, name)


def test_advanced_inventory_preserves_authoritative_tuple_ownership() -> None:
    import phospy.api._inventory as _inventory
    import phospy.api.advanced as advanced

    assert advanced.__all__ is _inventory.ADVANCED_PUBLIC_API
    assert isinstance(advanced.__all__, tuple)
    assert len(advanced.__all__) == 126
    with pytest.raises(TypeError):
        operator.setitem(cast(MutableSequence[object], advanced.__all__), 0, "changed")

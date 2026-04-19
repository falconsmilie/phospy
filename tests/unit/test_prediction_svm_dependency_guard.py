from __future__ import annotations

import builtins

import pytest

from phospy.prediction.svm import require_sklearn


def test_require_sklearn_error_reports_broken_standard_install(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original_import = builtins.__import__

    def failing_import(name: str, *args: object, **kwargs: object):
        if name == "sklearn" or name.startswith("sklearn."):
            raise ImportError("simulated sklearn import failure")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", failing_import)

    with pytest.raises(
        ImportError,
        match="part of PhosPy's standard install.*unexpected environment problem",
    ):
        require_sklearn()

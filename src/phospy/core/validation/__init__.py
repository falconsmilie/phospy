"""Internal validation domain for package-owned boundary checks.

Validation modules are support code for dataset construction, importers,
reference handling, and workflow boundaries. They are intentionally not part of
the stable public user API; public validation happens by constructing datasets
through supported builders or by running workflows.
"""

from __future__ import annotations

__all__: list[str] = []

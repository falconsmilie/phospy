"""Internal dataset validation support.

This package is not a public user API. Dataset validators support dataset
construction, importers, preprocessing, and workflow boundaries inside PhosPy.
Public users validate data by constructing datasets with supported builders or
by running workflows, not by importing validators from this package.

Internal code should import concrete validator modules directly, for example
``phospy.validation.datasets.site_metadata``.
"""

from __future__ import annotations

__all__: list[str] = []

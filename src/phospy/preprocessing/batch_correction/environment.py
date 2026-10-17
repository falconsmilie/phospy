"""Environment selection for reproducible batch correction."""

from __future__ import annotations

from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.models.environment import EnvironmentProvenance

BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES = (
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
)


def collect_batch_correction_environment_provenance(
    *,
    use_cache: bool = True,
) -> EnvironmentProvenance:
    """Collect versions relevant to numerical batch-correction provenance."""

    return collect_environment_provenance(
        dependency_names=BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES,
        use_cache=use_cache,
    )


__all__ = [
    "BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES",
    "collect_batch_correction_environment_provenance",
]

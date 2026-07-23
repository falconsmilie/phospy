"""Science-owned configuration policy exports."""
# pyright: reportUnsupportedDunderAll=false

from __future__ import annotations

from phospy.science.configs.dataset import *  # noqa: F403
from phospy.science.configs.dataset import __all__ as _dataset_all
from phospy.science.configs.differential import *  # noqa: F403
from phospy.science.configs.differential import __all__ as _differential_all
from phospy.science.configs.enrichment import *  # noqa: F403
from phospy.science.configs.enrichment import __all__ as _enrichment_all
from phospy.science.configs.kinase import *  # noqa: F403
from phospy.science.configs.kinase import __all__ as _kinase_all
from phospy.science.configs.localisation import *  # noqa: F403
from phospy.science.configs.localisation import __all__ as _localisation_all
from phospy.science.configs.prediction import *  # noqa: F403
from phospy.science.configs.prediction import __all__ as _prediction_all
from phospy.science.configs.preprocessing import *  # noqa: F403
from phospy.science.configs.preprocessing import __all__ as _preprocessing_all
from phospy.science.configs.reference_context import *  # noqa: F403
from phospy.science.configs.reference_context import __all__ as _reference_context_all
from phospy.science.configs.signalome import *  # noqa: F403
from phospy.science.configs.signalome import __all__ as _signalome_all

_exports: tuple[str, ...] = (
    tuple(_dataset_all)
    + tuple(_differential_all)
    + tuple(_enrichment_all)
    + tuple(_kinase_all)
    + tuple(_localisation_all)
    + tuple(_prediction_all)
    + tuple(_preprocessing_all)
    + tuple(_reference_context_all)
    + tuple(_signalome_all)
)

__all__ = _exports

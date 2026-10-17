"""Passive contracts published by the batch-correction capability."""
# pyright: reportUnsupportedDunderAll=false

from phospy.preparation_records.batch_correction import *  # noqa: F403
from phospy.preparation_records.batch_correction import __all__ as _passive_all
from phospy.preprocessing.batch_correction.contracts.configuration import *  # noqa: F403
from phospy.preprocessing.batch_correction.contracts.configuration import (
    __all__ as _configuration_all,
)
from phospy.preprocessing.batch_correction.contracts.missingness import *  # noqa: F403
from phospy.preprocessing.batch_correction.contracts.missingness import (
    __all__ as _missingness_all,
)
from phospy.preprocessing.batch_correction.contracts.output import *  # noqa: F403
from phospy.preprocessing.batch_correction.contracts.output import (
    __all__ as _output_all,
)
from phospy.preprocessing.batch_correction.contracts.request import *  # noqa: F403
from phospy.preprocessing.batch_correction.contracts.request import (
    __all__ as _request_all,
)
from phospy.preprocessing.batch_correction.contracts.results import *  # noqa: F403
from phospy.preprocessing.batch_correction.contracts.results import (
    __all__ as _results_all,
)
from phospy.preprocessing.batch_correction.contracts.stage import *  # noqa: F403
from phospy.preprocessing.batch_correction.contracts.stage import __all__ as _stage_all

__all__ = [
    *_configuration_all,
    *_missingness_all,
    *_request_all,
    *_passive_all,
    *_results_all,
    *_output_all,
    *_stage_all,
]

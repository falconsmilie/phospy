"""Supported public application surface.

This package is the thin public entry point for PhosPy. It should coordinate
domain operations without becoming a second implementation layer. Import the
workflow classes from here; module-local result bundle types remain available
from ``phospy.api.workflows`` for narrow internal and advanced use.
"""

from .workflows import (
    KinaseWorkflow,
    PredMatWorkflow,
    SignalomeWorkflow,
    SimpleKinaseWorkflow,
)

__all__ = [
    "KinaseWorkflow",
    "PredMatWorkflow",
    "SignalomeWorkflow",
    "SimpleKinaseWorkflow",
]

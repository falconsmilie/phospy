"""Supported public application surface.

This package is the thin public entry point for PhosPy. It should coordinate
domain operations without becoming a second implementation layer."""

from .workflows import (
    KinaseWorkflow,
    KinaseWorkflowResult,
    PredMatWorkflow,
    PredMatWorkflowResult,
    SignalomeWorkflow,
    SimpleKinaseWorkflow,
    SimpleKinaseWorkflowResult,
)

__all__ = [
    "KinaseWorkflow",
    "KinaseWorkflowResult",
    "PredMatWorkflow",
    "PredMatWorkflowResult",
    "SignalomeWorkflow",
    "SimpleKinaseWorkflow",
    "SimpleKinaseWorkflowResult",
]

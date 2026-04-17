"""I/O adapters and CLI plumbing for the supported rewrite lane."""

from phospy.io.adapters import DatasetFileInputs, build_dataset_from_files
from phospy.io.cli import main

__all__ = ["DatasetFileInputs", "build_dataset_from_files", "main"]

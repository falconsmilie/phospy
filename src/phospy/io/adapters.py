"""File-to-request adapters for the supported rewrite lane."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from phospy.api.builders import AnalysisReadyDatasetBuilder
from phospy.api.requests import DatasetBuildRequest
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.input import PhosPyInputError
from phospy.io.tables import read_table
from phospy.references.models import Organism, ReferencePreset


@dataclass(frozen=True, slots=True)
class DatasetFileInputs:
    """Supported file-based dataset inputs for the rewrite lane."""

    phospho_path: Path
    site_metadata_path: Path
    sample_metadata_path: Path | None = None
    total_path: Path | None = None
    organism: Organism | str | None = None


def dataset_build_request_from_files(inputs: DatasetFileInputs) -> DatasetBuildRequest:
    """Build a DataFrame-based dataset request from supported file inputs."""

    phospho = read_table(inputs.phospho_path)
    site_metadata = read_table(inputs.site_metadata_path)
    sample_metadata = (
        None
        if inputs.sample_metadata_path is None
        else read_table(inputs.sample_metadata_path)
    )
    total = None if inputs.total_path is None else read_table(inputs.total_path)
    organism = organism_from_value(inputs.organism)
    return DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=total,
        organism=organism,
    )


def build_dataset_from_files(
    inputs: DatasetFileInputs,
    *,
    builder: AnalysisReadyDatasetBuilder | None = None,
) -> AnalysisReadyPhosphoDataset:
    """Build ``AnalysisReadyPhosphoDataset`` from supported file inputs."""

    request = dataset_build_request_from_files(inputs)
    dataset_builder = builder or AnalysisReadyDatasetBuilder()
    return dataset_builder.run(request)


def organism_from_value(value: Organism | str | None) -> Organism | None:
    """Parse an optional organism token into ``Organism``."""

    if value is None:
        return None
    if isinstance(value, Organism):
        return value
    normalized = value.strip().lower()
    for organism in Organism:
        if organism.value == normalized:
            return organism
    supported = ", ".join(member.value for member in Organism)
    raise PhosPyInputError(
        f"unsupported organism '{value}'. supported organisms: {supported}"
    )


def reference_preset_from_value(value: ReferencePreset | str) -> ReferencePreset:
    """Parse a reference preset token into ``ReferencePreset``."""

    if isinstance(value, ReferencePreset):
        return value
    normalized = value.strip().lower()
    for preset in ReferencePreset:
        if preset.value == normalized:
            return preset
    supported = ", ".join(member.value for member in ReferencePreset)
    raise PhosPyInputError(
        f"unsupported reference preset '{value}'. supported presets: {supported}"
    )

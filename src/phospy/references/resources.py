"""Packaged bundled-reference resource loading helpers."""

from __future__ import annotations

from importlib import resources

import pandas as pd

from phospy.errors.references import ReferenceResolutionError, UnsupportedOrganismError
from phospy.references.models import Organism

_BUNDLED_DEFAULTS: dict[Organism, str] = {
    Organism.RAT: "l6_native",
}
_NON_BUNDLED_GUIDANCE = (
    "non-bundled organism lanes require a caller-supplied ReferenceBundle"
)


def supported_bundled_organisms() -> tuple[Organism, ...]:
    """Return organisms with packaged bundled references in this release."""

    return tuple(sorted(_BUNDLED_DEFAULTS, key=lambda organism: organism.value))


def bundled_reference_name_for_organism(organism: Organism) -> str:
    """Resolve one organism to its packaged bundled reference lane."""

    reference_name = _BUNDLED_DEFAULTS.get(organism)
    if reference_name is not None:
        return reference_name
    supported = ", ".join(item.value for item in supported_bundled_organisms())
    raise UnsupportedOrganismError(
        f"no bundled references are available for organism '{organism.value}' in the "
        f"current release; supported bundled organisms: {supported}; "
        f"{_NON_BUNDLED_GUIDANCE}"
    )


def load_bundled_kinase_substrate_map(organism: Organism) -> pd.DataFrame:
    """Load the packaged kinase-substrate map for one supported organism."""

    reference_name = bundled_reference_name_for_organism(organism)
    frame = _read_csv_resource(
        organism=organism,
        reference_name=reference_name,
        filename="substrate_map.csv",
    )
    required_columns = {"kinase", "site_id"}
    missing = required_columns.difference(frame.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ReferenceResolutionError(
            "bundled substrate_map.csv is missing required columns for "
            f"{organism.value}/{reference_name}: {missing_text}"
        )
    substrate_map = frame.loc[:, ["kinase", "site_id"]].rename(
        columns={"site_id": "substrate_site"}
    )
    substrate_map.loc[:, "kinase"] = (
        substrate_map.loc[:, "kinase"].astype(str).str.strip()
    )
    substrate_map.loc[:, "substrate_site"] = (
        substrate_map.loc[:, "substrate_site"].astype(str).str.strip()
    )
    return substrate_map


def load_bundled_site_sequences(organism: Organism) -> pd.DataFrame:
    """Load the packaged site-sequence table for one supported organism."""

    reference_name = bundled_reference_name_for_organism(organism)
    frame = _read_csv_resource(
        organism=organism,
        reference_name=reference_name,
        filename="site_sequences.csv",
    )
    required_columns = {"site_id", "centralized_sequence"}
    missing = required_columns.difference(frame.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ReferenceResolutionError(
            "bundled site_sequences.csv is missing required columns for "
            f"{organism.value}/{reference_name}: {missing_text}"
        )
    cleaned = frame.loc[:, ["site_id", "centralized_sequence"]].copy()
    cleaned.loc[:, "site_id"] = cleaned.loc[:, "site_id"].astype(str).str.strip()
    cleaned.loc[:, "centralized_sequence"] = (
        cleaned.loc[:, "centralized_sequence"].astype(str).str.strip()
    )
    site_sequences = cleaned.set_index("site_id")
    site_sequences.index.name = "site_id"
    return site_sequences.rename(columns={"centralized_sequence": "site_sequence"})


def _read_csv_resource(
    *,
    organism: Organism,
    reference_name: str,
    filename: str,
) -> pd.DataFrame:
    package_root = resources.files("phospy")
    resource = package_root.joinpath(
        "data",
        "reference_bundles",
        organism.value,
        reference_name,
        filename,
    )
    if not resource.is_file():
        raise ReferenceResolutionError(
            "bundled reference resource is missing for "
            f"{organism.value}/{reference_name}: {filename}"
        )
    with resources.as_file(resource) as resolved_path:
        return pd.read_csv(resolved_path)

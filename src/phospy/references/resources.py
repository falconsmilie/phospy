"""Packaged bundled-reference resource loading helpers."""

from __future__ import annotations

from importlib import resources

import pandas as pd

from phospy.errors.references import ReferenceResolutionError, UnsupportedOrganismError
from phospy.references.models import Organism
from phospy.site_ids import canonicalize_site_index, canonicalize_site_series

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
    if (substrate_map.loc[:, "kinase"] == "").any():
        raise ReferenceResolutionError(
            "bundled substrate_map.csv contains blank kinase values for "
            f"{organism.value}/{reference_name}"
        )
    substrate_map.loc[:, "substrate_site"] = canonicalize_site_series(
        substrate_map.loc[:, "substrate_site"],
        field_name="bundled reference substrate_map.csv substrate_site",
        error_type=ReferenceResolutionError,
    )
    return substrate_map.drop_duplicates(
        subset=["kinase", "substrate_site"],
        ignore_index=True,
    )


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
    cleaned.loc[:, "site_id"] = canonicalize_site_series(
        cleaned.loc[:, "site_id"],
        field_name="bundled reference site_sequences.csv site_id",
        error_type=ReferenceResolutionError,
    )
    cleaned.loc[:, "centralized_sequence"] = (
        cleaned.loc[:, "centralized_sequence"].astype(str).str.strip()
    )
    if (cleaned.loc[:, "centralized_sequence"] == "").any():
        raise ReferenceResolutionError(
            "bundled site_sequences.csv contains blank centralized_sequence values "
            f"for {organism.value}/{reference_name}"
        )
    site_sequences = cleaned.set_index("site_id")
    site_sequences.index = canonicalize_site_index(
        site_sequences.index,
        field_name="bundled reference site_sequences.csv site_id",
        error_type=ReferenceResolutionError,
        index_name="site_id",
    )
    return site_sequences.rename(columns={"centralized_sequence": "site_sequence"})


def load_bundled_motif_scores(organism: Organism) -> pd.DataFrame | None:
    """Load optional bundled motif-score table for one supported organism."""

    reference_name = bundled_reference_name_for_organism(organism)
    frame = _read_optional_csv_resource(
        organism=organism,
        reference_name=reference_name,
        filename="motif_scores.csv",
    )
    if frame is None:
        return None
    if frame.empty:
        raise ReferenceResolutionError(
            "bundled motif_scores.csv must be non-empty for "
            f"{organism.value}/{reference_name}"
        )
    if "site_id" in frame.columns:
        score_frame = frame.set_index("site_id")
    else:
        index_name = frame.columns[0]
        score_frame = frame.set_index(index_name)
    score_frame.index = canonicalize_site_index(
        score_frame.index,
        field_name="bundled reference motif_scores.csv site_id",
        error_type=ReferenceResolutionError,
        index_name="site_id",
    )
    score_frame.columns = score_frame.columns.astype(str)
    score_frame.columns.name = "kinase"
    try:
        return score_frame.astype(float)
    except ValueError as exc:
        raise ReferenceResolutionError(
            "bundled motif_scores.csv must contain numeric score values for "
            f"{organism.value}/{reference_name}"
        ) from exc


def load_bundled_motif_sizes(organism: Organism) -> pd.Series | None:
    """Load optional bundled motif-size metadata for one supported organism."""

    reference_name = bundled_reference_name_for_organism(organism)
    frame = _read_optional_csv_resource(
        organism=organism,
        reference_name=reference_name,
        filename="motif_sizes.csv",
    )
    if frame is None:
        return None
    required_columns = {"kinase", "motif_size"}
    missing = required_columns.difference(frame.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ReferenceResolutionError(
            "bundled motif_sizes.csv is missing required columns for "
            f"{organism.value}/{reference_name}: {missing_text}"
        )
    cleaned = frame.loc[:, ["kinase", "motif_size"]].copy()
    cleaned.loc[:, "kinase"] = cleaned.loc[:, "kinase"].astype(str).str.strip()
    if (cleaned.loc[:, "kinase"] == "").any():
        raise ReferenceResolutionError(
            "bundled motif_sizes.csv contains blank kinase values for "
            f"{organism.value}/{reference_name}"
        )
    try:
        sizes = cleaned.set_index("kinase").loc[:, "motif_size"].astype(float)
    except ValueError as exc:
        raise ReferenceResolutionError(
            "bundled motif_sizes.csv must contain numeric motif_size values for "
            f"{organism.value}/{reference_name}"
        ) from exc
    sizes.index.name = "kinase"
    sizes.name = "motif_size"
    return sizes


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


def _read_optional_csv_resource(
    *,
    organism: Organism,
    reference_name: str,
    filename: str,
) -> pd.DataFrame | None:
    package_root = resources.files("phospy")
    resource = package_root.joinpath(
        "data",
        "reference_bundles",
        organism.value,
        reference_name,
        filename,
    )
    if not resource.is_file():
        return None
    with resources.as_file(resource) as resolved_path:
        return pd.read_csv(resolved_path)

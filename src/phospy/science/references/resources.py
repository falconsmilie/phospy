"""Packaged bundled-reference resource loading helpers."""

from __future__ import annotations

import json
from functools import cache
from importlib import resources

import pandas as pd

from phospy.errors.references import ReferenceResolutionError, UnsupportedOrganismError
from phospy.science.references.models import (
    BundledReferenceLane,
    Organism,
    ReferenceManifest,
)
from phospy.science.references.validation import load_reference_manifest
from phospy.science.sites.identifiers import (
    canonicalize_site_index,
    canonicalize_site_series,
)

_BUNDLED_DEFAULTS: dict[Organism, str] = {
    Organism.RAT: "l6_native",
}
_MANIFEST_FILENAME = "manifest.json"
_REFERENCE_BUNDLE_DOCS_URL = "https://phospy.com/docs/api/guide/#references"
_EXPLICIT_REFERENCE_BUNDLE_GUIDANCE = (
    "provide an explicit ReferenceBundle by passing "
    "ReferenceBundle(organism=..., kinase_substrate_map=..., site_sequences=...) "
    "as the references input"
)


def supported_bundled_organisms() -> tuple[Organism, ...]:
    """Return organisms with packaged bundled references in this release."""

    return tuple(sorted(_BUNDLED_DEFAULTS, key=lambda organism: organism.value))


def available_bundled_reference_lanes() -> tuple[BundledReferenceLane, ...]:
    """Return inventory metadata for packaged bundled reference lanes."""

    lanes: list[BundledReferenceLane] = []
    for organism in supported_bundled_organisms():
        manifest = load_bundled_reference_manifest(organism)
        lanes.append(
            BundledReferenceLane(
                organism=organism,
                bundle_id=manifest.bundle_id,
                source_name=manifest.source_name,
                source_version=manifest.reference_version,
                retrieved_at=manifest.retrieved_at,
                redistribution_status=manifest.redistribution_status.value,
                supports=manifest.supports,
                limitations=manifest.limitations,
            )
        )
    return tuple(lanes)


def bundled_reference_name_for_organism(organism: Organism) -> str:
    """Resolve one organism to its packaged bundled reference lane."""

    reference_name = _BUNDLED_DEFAULTS.get(organism)
    if reference_name is not None:
        return reference_name
    supported = _format_supported_bundled_organisms()
    raise UnsupportedOrganismError(
        f"no bundled references are available for organism '{organism.value}' in the "
        f"current release; supported bundled organisms: {supported}; "
        f"{_EXPLICIT_REFERENCE_BUNDLE_GUIDANCE}; "
        f"reference-bundle docs: {_REFERENCE_BUNDLE_DOCS_URL}"
    )


def clear_bundled_reference_manifest_cache() -> None:
    """Clear cached bundled manifest parsing (primarily for tests)."""

    _load_bundled_reference_manifest_cached.cache_clear()


def load_bundled_reference_manifest(organism: Organism) -> ReferenceManifest:
    """Load and validate the machine-readable bundled manifest for an organism."""

    return _load_bundled_reference_manifest_cached(organism)


@cache
def _load_bundled_reference_manifest_cached(organism: Organism) -> ReferenceManifest:
    reference_name = bundled_reference_name_for_organism(organism)
    bundle_resource = _bundled_reference_bundle_resource(
        organism=organism,
        reference_name=reference_name,
    )
    if not bundle_resource.is_dir():
        raise ReferenceResolutionError(
            "bundled reference bundle directory is missing for "
            f"{organism.value}/{reference_name}"
        )
    with resources.as_file(bundle_resource) as bundle_root:
        return load_reference_manifest(
            bundle_root.joinpath(_MANIFEST_FILENAME),
            bundle_root=bundle_root,
            bundled=False,
            require_redistribution_allowed=False,
            require_all_files_listed=True,
        )


def load_bundled_kinase_substrate_map(organism: Organism) -> pd.DataFrame:
    """Load the packaged kinase-substrate map for one supported organism."""

    load_bundled_reference_manifest(organism)
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

    load_bundled_reference_manifest(organism)
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

    load_bundled_reference_manifest(organism)
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

    load_bundled_reference_manifest(organism)
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


def _format_supported_bundled_organisms() -> str:
    supported = tuple(item.value for item in supported_bundled_organisms())
    if not supported:
        return "(none)"
    return ", ".join(supported)


def _read_csv_resource(
    *,
    organism: Organism,
    reference_name: str,
    filename: str,
) -> pd.DataFrame:
    resource = _bundled_reference_resource(
        organism=organism,
        reference_name=reference_name,
        filename=filename,
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
    resource = _bundled_reference_resource(
        organism=organism,
        reference_name=reference_name,
        filename=filename,
    )
    if not resource.is_file():
        return None
    with resources.as_file(resource) as resolved_path:
        return pd.read_csv(resolved_path)


def _read_json_resource(
    *,
    organism: Organism,
    reference_name: str,
    filename: str,
) -> object:
    resource = _bundled_reference_resource(
        organism=organism,
        reference_name=reference_name,
        filename=filename,
    )
    if not resource.is_file():
        raise ReferenceResolutionError(
            "bundled reference resource is missing for "
            f"{organism.value}/{reference_name}: {filename}"
        )
    with resources.as_file(resource) as resolved_path:
        try:
            return json.loads(resolved_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise ReferenceResolutionError(
                "bundled reference manifest is not valid JSON for "
                f"{organism.value}/{reference_name}: {filename}"
            ) from exc


def _bundled_reference_resource(
    *,
    organism: Organism,
    reference_name: str,
    filename: str,
):
    package_root = resources.files("phospy")
    return (
        package_root.joinpath("data")
        .joinpath("reference_bundles")
        .joinpath(organism.value)
        .joinpath(reference_name)
        .joinpath(filename)
    )


def _bundled_reference_bundle_resource(
    *,
    organism: Organism,
    reference_name: str,
):
    package_root = resources.files("phospy")
    return (
        package_root.joinpath("data")
        .joinpath("reference_bundles")
        .joinpath(organism.value)
        .joinpath(reference_name)
    )

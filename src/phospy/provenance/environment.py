"""Environment provenance construction."""

from __future__ import annotations

import hashlib
import importlib
import locale
import os
import platform
from collections.abc import Mapping
from datetime import datetime
from functools import lru_cache
from importlib import metadata
from pathlib import Path
from types import ModuleType
from typing import cast

from phospy.provenance.models import (
    ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2,
    EnvironmentProvenance,
    JsonValue,
)

CORE_ENVIRONMENT_DEPENDENCIES = ("numpy", "pandas", "scipy", "scikit-learn")
BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES = CORE_ENVIRONMENT_DEPENDENCIES
OPTIONAL_ENVIRONMENT_DEPENDENCIES = ("pyarrow", "openpyxl")
DEFAULT_ENVIRONMENT_DEPENDENCIES = (
    *CORE_ENVIRONMENT_DEPENDENCIES,
    *OPTIONAL_ENVIRONMENT_DEPENDENCIES,
)
THREAD_ENVIRONMENT_VARIABLES = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
CONSTRAINTS_FILE_PATTERNS = (
    "constraints.txt",
    "requirements.txt",
    "requirements-dev.txt",
    "requirements-ci.txt",
    "constraints/*.txt",
    "constraints/*.in",
    "poetry.lock",
    "Pipfile.lock",
    "uv.lock",
    "environment.yml",
    "conda-lock.yml",
)
_CONSTRAINTS_HASH_ALGORITHM = "sha256"


def collect_environment_provenance(
    *,
    package_name: str = "phospy",
    dependency_names: tuple[str, ...] = DEFAULT_ENVIRONMENT_DEPENDENCIES,
    use_cache: bool = True,
) -> EnvironmentProvenance:
    """Collect package/python/dependency versions for run provenance."""

    if use_cache:
        return _collect_environment_provenance_cached(package_name, dependency_names)
    return _collect_environment_provenance_uncached(
        package_name=package_name,
        dependency_names=dependency_names,
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


def clear_environment_provenance_cache() -> None:
    """Clear process-local cached environment provenance snapshot."""

    _collect_environment_provenance_cached.cache_clear()


@lru_cache(maxsize=16)
def _collect_environment_provenance_cached(
    package_name: str,
    dependency_names: tuple[str, ...],
) -> EnvironmentProvenance:
    return _collect_environment_provenance_uncached(
        package_name=package_name,
        dependency_names=dependency_names,
    )


def _collect_environment_provenance_uncached(
    *,
    package_name: str,
    dependency_names: tuple[str, ...],
) -> EnvironmentProvenance:
    return EnvironmentProvenance(
        schema_version=ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2,
        package_name=package_name,
        package_version=_package_version(package_name),
        python_version=platform.python_version(),
        dependency_versions={
            dependency: _distribution_version(dependency)
            for dependency in dependency_names
        },
        platform=_platform_provenance(),
        blas_lapack=_blas_lapack_provenance(),
        thread_environment=_thread_environment_provenance(),
        timezone=_timezone_provenance(),
        locale=_locale_provenance(),
        constraints_fingerprint=_constraints_fingerprint(),
    )


def _distribution_version(distribution_name: str) -> str | None:
    try:
        return metadata.version(distribution_name)
    except metadata.PackageNotFoundError:
        return None


def _package_version(package_name: str) -> str:
    return (
        _distribution_version(package_name)
        or _project_metadata_version(package_name)
        or "unavailable"
    )


def _project_metadata_version(package_name: str) -> str | None:
    toml_parser = _toml_parser()
    if toml_parser is None:
        return None
    try:
        pyproject = _project_root_from_module() / "pyproject.toml"
        payload = cast(
            Mapping[str, object],
            toml_parser.loads(pyproject.read_text(encoding="utf-8")),
        )
    except Exception:
        return None
    project = payload.get("project")
    if not isinstance(project, Mapping):
        return None
    if _normalize_optional_value(project.get("name")) != package_name:
        return None
    return _normalize_optional_value(project.get("version"))


def _toml_parser() -> ModuleType | None:
    for module_name in ("tomllib", "tomli"):
        try:
            return importlib.import_module(module_name)
        except ModuleNotFoundError:
            continue
    return None


def _platform_provenance() -> dict[str, str]:
    return {
        "platform": _normalize_platform_value(platform.platform()),
        "system": _normalize_platform_value(platform.system()),
        "release": _normalize_platform_value(platform.release()),
        "version": _normalize_platform_value(platform.version()),
        "machine": _normalize_platform_value(platform.machine()),
        "processor": _normalize_platform_value(platform.processor()),
        "python_implementation": _normalize_platform_value(
            platform.python_implementation()
        ),
    }


def _normalize_platform_value(value: str) -> str:
    normalized = str(value).strip()
    if normalized:
        return normalized
    return "unknown"


def _normalize_optional_value(value: object) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip()
    if normalized:
        return normalized
    return None


def _blas_lapack_provenance() -> dict[str, JsonValue]:
    payload = {
        "source": "unavailable",
        "blas_name": None,
        "blas_version": None,
        "blas_detection_method": None,
        "blas_openblas_configuration": None,
        "lapack_name": None,
        "lapack_version": None,
        "lapack_detection_method": None,
        "lapack_openblas_configuration": None,
    }
    try:
        import numpy as np
    except Exception:
        payload["source"] = "numpy_unavailable"
        return payload

    payload["source"] = "numpy_config"
    numpy_config = getattr(np, "__config__", None)
    if numpy_config is None:
        payload["source"] = "numpy_config_missing"
        return payload
    config_payload = getattr(numpy_config, "CONFIG", None)
    if isinstance(config_payload, Mapping):
        build_dependencies = config_payload.get("Build Dependencies")
        if isinstance(build_dependencies, Mapping):
            _update_backend_payload(
                payload=payload,
                backend_payload=build_dependencies.get("blas"),
                key_prefix="blas",
            )
            _update_backend_payload(
                payload=payload,
                backend_payload=build_dependencies.get("lapack"),
                key_prefix="lapack",
            )
            return payload

    legacy_get_info = getattr(numpy_config, "get_info", None)
    if callable(legacy_get_info):
        payload["source"] = "numpy_get_info"
        _update_legacy_backend_payload(
            payload=payload,
            backend_payload=legacy_get_info("blas_opt_info"),
            key_prefix="blas",
        )
        _update_legacy_backend_payload(
            payload=payload,
            backend_payload=legacy_get_info("lapack_opt_info"),
            key_prefix="lapack",
        )
    else:
        payload["source"] = "numpy_backend_unavailable"
    return payload


def _update_backend_payload(
    *,
    payload: dict[str, JsonValue],
    backend_payload: object,
    key_prefix: str,
) -> None:
    if not isinstance(backend_payload, Mapping):
        return
    payload[f"{key_prefix}_name"] = _normalize_optional_value(
        backend_payload.get("name")
    )
    payload[f"{key_prefix}_version"] = _normalize_optional_value(
        backend_payload.get("version")
    )
    payload[f"{key_prefix}_detection_method"] = _normalize_optional_value(
        backend_payload.get("detection method")
    )
    payload[f"{key_prefix}_openblas_configuration"] = _normalize_optional_value(
        backend_payload.get("openblas configuration")
    )


def _update_legacy_backend_payload(
    *,
    payload: dict[str, JsonValue],
    backend_payload: object,
    key_prefix: str,
) -> None:
    if not isinstance(backend_payload, Mapping):
        return
    libraries = backend_payload.get("libraries")
    normalized_libraries = None
    if isinstance(libraries, list):
        normalized_libraries = ",".join(
            item
            for item in (
                _normalize_optional_value(candidate) for candidate in libraries
            )
            if item is not None
        )
    payload[f"{key_prefix}_name"] = normalized_libraries
    payload[f"{key_prefix}_version"] = _normalize_optional_value(
        backend_payload.get("version")
    )
    payload[f"{key_prefix}_detection_method"] = _normalize_optional_value(
        backend_payload.get("language")
    )
    payload[f"{key_prefix}_openblas_configuration"] = _normalize_optional_value(
        backend_payload.get("define_macros")
    )


def _thread_environment_provenance() -> dict[str, str | None]:
    return {
        variable: _normalize_optional_value(os.environ.get(variable))
        for variable in THREAD_ENVIRONMENT_VARIABLES
    }


def _timezone_provenance() -> str | None:
    timezone_env = _normalize_optional_value(os.environ.get("TZ"))
    local_timezone_name = _normalize_optional_value(
        datetime.now().astimezone().tzname()
    )
    utc_offset = datetime.now().astimezone().utcoffset()
    if timezone_env is not None:
        return timezone_env
    if local_timezone_name is not None:
        return local_timezone_name
    if utc_offset is not None:
        total_minutes = int(utc_offset.total_seconds() // 60)
        sign = "+" if total_minutes >= 0 else "-"
        absolute_minutes = abs(total_minutes)
        return f"UTC{sign}{absolute_minutes // 60:02d}:{absolute_minutes % 60:02d}"
    return None


def _locale_provenance() -> dict[str, str | None]:
    language_code: str | None
    encoding: str | None
    try:
        language_code, encoding = locale.getlocale()
    except Exception:
        language_code = None
        encoding = None
    try:
        locale_all = _normalize_optional_value(locale.setlocale(locale.LC_ALL, None))
    except Exception:
        locale_all = None
    try:
        preferred_encoding = _normalize_optional_value(
            locale.getpreferredencoding(False)
        )
    except Exception:
        preferred_encoding = None
    return {
        "language_code": _normalize_optional_value(language_code),
        "encoding": _normalize_optional_value(encoding),
        "lc_all": locale_all,
        "preferred_encoding": preferred_encoding,
    }


def _constraints_fingerprint() -> dict[str, str | None]:
    payload: dict[str, str | None] = {
        "algorithm": _CONSTRAINTS_HASH_ALGORITHM,
        "value": None,
        "sources": None,
    }
    discovered_files = _discover_constraints_files()
    if not discovered_files:
        return payload
    hasher = hashlib.new(_CONSTRAINTS_HASH_ALGORITHM)
    source_labels: list[str] = []
    for path in discovered_files:
        try:
            content = path.read_bytes()
        except OSError:
            continue
        source_labels.append(str(path))
        hasher.update(str(path).encode("utf-8", errors="replace"))
        hasher.update(b"\n")
        hasher.update(content)
        hasher.update(b"\n")
    if not source_labels:
        return payload
    payload["value"] = hasher.hexdigest()
    payload["sources"] = ",".join(source_labels)
    return payload


def _discover_constraints_files() -> tuple[Path, ...]:
    roots: list[Path] = [Path.cwd()]
    package_root = _project_root_from_module()
    if package_root not in roots:
        roots.append(package_root)

    discovered: set[Path] = set()
    for root in roots:
        for pattern in CONSTRAINTS_FILE_PATTERNS:
            discovered.update(candidate.resolve() for candidate in root.glob(pattern))
    return tuple(sorted(path for path in discovered if path.is_file()))


def _project_root_from_module() -> Path:
    # src/phospy/provenance/environment.py -> repository root
    return Path(__file__).resolve().parents[3]


__all__ = [
    "BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES",
    "CORE_ENVIRONMENT_DEPENDENCIES",
    "DEFAULT_ENVIRONMENT_DEPENDENCIES",
    "OPTIONAL_ENVIRONMENT_DEPENDENCIES",
    "THREAD_ENVIRONMENT_VARIABLES",
    "clear_environment_provenance_cache",
    "collect_batch_correction_environment_provenance",
    "collect_environment_provenance",
]

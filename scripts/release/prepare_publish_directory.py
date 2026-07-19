from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

ATTESTATION_SCHEMA = "phospy.release-attestation/v1"


class PublishPreparationError(RuntimeError):
    """Raised when publication inputs do not match the release attestation."""


def prepare_publish_directory(
    *,
    attestation_path: Path,
    dist_dir: Path,
    output_dir: Path,
) -> Path:
    attestation = _read_json_object(attestation_path)
    _require(
        attestation.get("schema") == ATTESTATION_SCHEMA,
        "release attestation schema mismatch",
    )
    _require(
        attestation.get("status") == "success", "release attestation is not success"
    )
    artifacts = _attested_artifacts(attestation)
    _reject_unattested_distribution_artifacts(dist_dir, artifacts)
    _prepare_empty_directory(output_dir)
    for kind, record in sorted(artifacts.items()):
        source = dist_dir / record["filename"]
        _require(source.is_file(), f"attested {kind} artifact is missing: {source}")
        observed_sha256 = _file_sha256(source)
        _require(
            observed_sha256 == record["sha256"],
            f"attested {kind} digest changed: expected {record['sha256']}, "
            f"got {observed_sha256}",
        )
        shutil.copy2(source, output_dir / source.name)
    return output_dir


def _attested_artifacts(attestation: Mapping[str, object]) -> dict[str, dict[str, str]]:
    artifacts = attestation.get("artifacts")
    _require(isinstance(artifacts, Mapping), "attestation artifacts must be an object")
    result: dict[str, dict[str, str]] = {}
    for kind, suffix in {"wheel": ".whl", "sdist": ".tar.gz"}.items():
        record = artifacts.get(kind)
        _require(isinstance(record, Mapping), f"attested {kind} record is missing")
        filename = _required_text(record.get("filename"), field_name=f"{kind} filename")
        _require(
            filename.endswith(suffix), f"attested {kind} filename has wrong suffix"
        )
        sha256 = _required_text(record.get("sha256"), field_name=f"{kind} sha256")
        _require_sha256_hex(sha256, field_name=f"{kind} sha256")
        result[kind] = {"filename": filename, "sha256": sha256}
    filenames = [record["filename"] for record in result.values()]
    _require(
        len(filenames) == len(set(filenames)), "attested artifact filenames duplicate"
    )
    return result


def _reject_unattested_distribution_artifacts(
    dist_dir: Path,
    artifacts: Mapping[str, Mapping[str, str]],
) -> None:
    _require(dist_dir.is_dir(), f"distribution directory is missing: {dist_dir}")
    attested_filenames = {record["filename"] for record in artifacts.values()}
    present_distribution_filenames = {
        path.name
        for path in dist_dir.iterdir()
        if path.is_file()
        and (path.name.endswith(".whl") or path.name.endswith(".tar.gz"))
    }
    unknown = sorted(present_distribution_filenames - attested_filenames)
    _require(
        unknown == [],
        "unattested distribution artifacts present: " + ", ".join(unknown),
    )
    missing = sorted(attested_filenames - present_distribution_filenames)
    _require(
        missing == [], "attested distribution artifacts missing: " + ", ".join(missing)
    )


def _prepare_empty_directory(output_dir: Path) -> None:
    if output_dir.exists():
        _require(
            output_dir.is_dir(),
            f"publish output exists and is not a directory: {output_dir}",
        )
        for child in output_dir.iterdir():
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()
    else:
        output_dir.mkdir(parents=True)


def _read_json_object(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise PublishPreparationError(
            f"attestation JSON is malformed: {path.as_posix()}: {exc}"
        ) from exc
    _require(
        isinstance(payload, dict), f"attestation must contain a JSON object: {path}"
    )
    return payload


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _required_text(value: object, *, field_name: str) -> str:
    _require(
        isinstance(value, str) and bool(value.strip()),
        f"{field_name} must be non-empty text",
    )
    return value.strip()


def _require_sha256_hex(value: str, *, field_name: str) -> None:
    _require(
        len(value) == 64
        and all(character in "0123456789abcdef" for character in value),
        f"{field_name} must be 64 lowercase hexadecimal characters",
    )


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PublishPreparationError(message)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    path = prepare_publish_directory(
        attestation_path=Path(args.attestation),
        dist_dir=Path(args.dist_dir),
        output_dir=Path(args.output_dir),
    )
    print(path.as_posix())
    return 0


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Copy exactly attested distributions into a clean publish directory."
    )
    parser.add_argument(
        "--attestation",
        default="build/release/release-attestation.json",
        help="Successful release attestation JSON.",
    )
    parser.add_argument(
        "--dist-dir",
        default="dist",
        help="Directory containing distribution artifacts.",
    )
    parser.add_argument(
        "--output-dir",
        default="build/publish",
        help="Clean output directory for the publishing action.",
    )
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(main())

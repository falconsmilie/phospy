from __future__ import annotations

import argparse
import hashlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

ATTESTATION_SCHEMA = "phospy.release-attestation/v1"


class ReleaseAttestationVerificationError(RuntimeError):
    """Raised when an existing release attestation no longer matches evidence."""


def verify_release_attestation(
    *,
    attestation_path: Path,
    evidence_root: Path,
    artifact_dir: Path,
) -> None:
    attestation = _read_json_object(attestation_path)
    _require(
        attestation.get("schema") == ATTESTATION_SCHEMA,
        "release attestation schema mismatch",
    )
    _require(
        attestation.get("status") == "success", "release attestation is not success"
    )

    _verify_file_record(attestation["source_identity"], evidence_root=evidence_root)
    _verify_file_record(
        attestation["release_evidence_policy"],
        evidence_root=evidence_root,
    )
    _verify_file_record(attestation["build_manifest"], evidence_root=evidence_root)

    source_checks = attestation.get("source_checks")
    _require(isinstance(source_checks, list), "source_checks must be an array")
    for record in source_checks:
        _verify_file_record(record, evidence_root=evidence_root)

    matrix = attestation.get("artifact_verification_matrix")
    _require(isinstance(matrix, list), "artifact_verification_matrix must be an array")
    for record in matrix:
        _verify_file_record(record, evidence_root=evidence_root)

    artifacts = attestation.get("artifacts")
    _require(isinstance(artifacts, Mapping), "artifacts must be an object")
    for kind in ("wheel", "sdist"):
        record = artifacts.get(kind)
        _require(isinstance(record, Mapping), f"attested {kind} record is missing")
        filename = _required_text(record.get("filename"), field_name=f"{kind} filename")
        expected_sha256 = _required_text(
            record.get("sha256"),
            field_name=f"{kind} sha256",
        )
        _require_sha256_hex(expected_sha256, field_name=f"{kind} sha256")
        path = (artifact_dir / filename).resolve()
        _require(path.is_file(), f"attested {kind} artifact is missing: {filename}")
        actual_sha256 = _file_sha256(path)
        _require(
            actual_sha256 == expected_sha256,
            f"attested {kind} digest changed: expected {expected_sha256}, "
            f"got {actual_sha256}",
        )


def _verify_file_record(record: object, *, evidence_root: Path) -> None:
    _require(isinstance(record, Mapping), "attested file record must be an object")
    filename = _required_text(record.get("filename"), field_name="attested filename")
    expected_digest = _required_text(record.get("digest"), field_name="attested digest")
    _require_sha256_digest(expected_digest, field_name=f"{filename} digest")
    path = (evidence_root / filename).resolve()
    _require(path.is_file(), f"attested evidence file is missing: {filename}")
    actual_digest = "sha256:" + _file_sha256(path)
    _require(
        actual_digest == expected_digest,
        f"attested evidence digest changed for {filename}: expected "
        f"{expected_digest}, got {actual_digest}",
    )


def _read_json_object(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ReleaseAttestationVerificationError(
            f"JSON file is malformed: {path.as_posix()}: {exc}"
        ) from exc
    _require(isinstance(payload, dict), f"JSON file must contain an object: {path}")
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


def _require_sha256_digest(value: str, *, field_name: str) -> None:
    _require(value.startswith("sha256:"), f"{field_name} must use sha256: prefix")
    _require_sha256_hex(value.removeprefix("sha256:"), field_name=field_name)


def _require_sha256_hex(value: str, *, field_name: str) -> None:
    _require(
        len(value) == 64
        and all(character in "0123456789abcdef" for character in value),
        f"{field_name} must be 64 lowercase hexadecimal characters",
    )


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ReleaseAttestationVerificationError(message)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    verify_release_attestation(
        attestation_path=Path(args.attestation),
        evidence_root=Path(args.evidence_root),
        artifact_dir=Path(args.artifact_dir),
    )
    print("release attestation verified")
    return 0


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify a release attestation against retained evidence files."
    )
    parser.add_argument(
        "--attestation",
        default="build/release/release-attestation.json",
        help="Release attestation JSON path.",
    )
    parser.add_argument(
        "--evidence-root",
        default=".",
        help="Root used to resolve attested evidence filenames.",
    )
    parser.add_argument(
        "--artifact-dir",
        default="dist",
        help="Directory containing the attested wheel and sdist.",
    )
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(main())

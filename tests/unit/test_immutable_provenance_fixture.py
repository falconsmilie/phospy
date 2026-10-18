from __future__ import annotations

import hashlib
from pathlib import Path

FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "architecture"
    / "batch_correction_provenance_regression.json"
)
EXPECTED_BYTE_LENGTH = 4283
EXPECTED_SHA256 = "868f68b01f48f59ec4c79837d21d756b307bab96503bfa077d8e293ee017cae7"
EXPECTED_GIT_BLOB_SHA1 = "54ee5cd27e5e1ccf5cf6bd13cffff94b97164f79"


def test_provenance_fixture_retains_pre_regeneration_git_identity() -> None:
    content = FIXTURE.read_bytes()
    git_object = b"blob " + str(len(content)).encode("ascii") + b"\0" + content

    assert len(content) == EXPECTED_BYTE_LENGTH
    assert hashlib.sha256(content).hexdigest() == EXPECTED_SHA256
    assert hashlib.sha1(git_object, usedforsecurity=False).hexdigest() == (
        EXPECTED_GIT_BLOB_SHA1
    )

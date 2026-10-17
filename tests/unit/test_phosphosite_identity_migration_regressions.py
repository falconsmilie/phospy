from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from phospy.errors.validation import PhosPyValidationError
from phospy.phosphosites.site_keys import (
    ProteinScopedPhosphositeKey,
    decode_site_key,
    encode_site_key,
)
from phospy.references.context_compatibility import (
    REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE,
    ReferenceContextCompatibilityWarning,
    validate_reference_context_compatibility,
)
from phospy.references.provenance.context import ReferenceContext
from phospy.taxonomy import Organism

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
LEFT_CONTEXT_ID = (
    "reference-context-v1:"
    "5dca77948e02b21e5359ed97dd54aac6d4cda0e5ef7a37445a15a2989a9e4c5a"
)
RIGHT_CONTEXT_ID = (
    "reference-context-v1:"
    "4d090c4841597048b3fb378bcb0ca666df49669bbb76c61c2fad83f78f9f145c"
)


def _left_context() -> ReferenceContext:
    return ReferenceContext(
        organism="rat",
        protein_namespace="gene_symbol",
        source_name="unit-reference",
        source_version="v1",
        proteome_version="rnor6",
        reference_table_sha256="a" * 64,
    )


def _right_context() -> ReferenceContext:
    return ReferenceContext(
        organism="mouse",
        protein_namespace="uniprot_accession",
        source_name="other-reference",
        source_version="v2",
        proteome_version="gencode-m1",
        reference_table_sha256="b" * 64,
    )


@pytest.mark.parametrize(
    ("isoform_id", "expected"),
    [
        (
            None,
            b"phospy:v1|organism=rat|protein_namespace=ns%7C%3D%2F|"
            b"protein_identifier=P%7C%3D%2F|residue=S|position=42",
        ),
        (
            "iso|=/2",
            b"phospy:v1|organism=rat|protein_namespace=ns%7C%3D%2F|"
            b"protein_identifier=P%7C%3D%2F|residue=S|position=42|"
            b"isoform_id=iso%7C%3D%2F2",
        ),
    ],
)
def test_site_key_encoding_bytes_field_order_and_round_trip(
    isoform_id: str | None,
    expected: bytes,
) -> None:
    key = ProteinScopedPhosphositeKey(
        organism=Organism.RAT,
        protein_namespace="ns|=/",
        protein_identifier="P|=/",
        residue="S",
        position=42,
        isoform_id=isoform_id,
    )

    encoded = encode_site_key(key)

    assert encoded.encode("ascii") == expected
    decoded = decode_site_key(
        encoded,
        field_name="site_key",
        error_type=ValueError,
    )
    assert decoded == key
    assert encode_site_key(decoded) == encoded


def test_reference_context_equality_and_complete_mismatch_message() -> None:
    left = _left_context()
    assert left.reference_context_id == LEFT_CONTEXT_ID
    assert (
        validate_reference_context_compatibility(
            left,
            _left_context(),
            operation="unit equality",
        )
        is None
    )

    with pytest.raises(PhosPyValidationError) as exc_info:
        validate_reference_context_compatibility(
            left,
            _right_context(),
            operation="  unit full mismatch  ",
        )

    assert str(exc_info.value) == (
        "reference-context compatibility failed for operation='unit full mismatch': "
        "incompatible reference contexts; mismatched_fields=organism, "
        "protein_namespace, source_name, source_version, proteome_version, "
        "reference_table_sha256; "
        f"left={{reference_context_id='{LEFT_CONTEXT_ID}', "
        "organism=<Organism.RAT: 'rat'>, protein_namespace='gene_symbol', "
        "source_name='unit-reference', source_version='v1', "
        "proteome_version='rnor6', "
        f"reference_table_sha256='{'a' * 64}'}}; "
        f"right={{reference_context_id='{RIGHT_CONTEXT_ID}', "
        "organism=<Organism.MOUSE: 'mouse'>, "
        "protein_namespace='uniprot_accession', source_name='other-reference', "
        "source_version='v2', proteome_version='gencode-m1', "
        f"reference_table_sha256='{'b' * 64}'}}"
    )


@pytest.mark.parametrize(
    ("left", "right", "missing_sides", "left_id", "right_id"),
    [
        (None, None, ("left", "right"), None, None),
        (None, _right_context(), ("left",), None, RIGHT_CONTEXT_ID),
        (_left_context(), None, ("right",), LEFT_CONTEXT_ID, None),
    ],
)
def test_reference_context_required_and_unknown_warning_payloads(
    left: ReferenceContext | None,
    right: ReferenceContext | None,
    missing_sides: tuple[str, ...],
    left_id: str | None,
    right_id: str | None,
) -> None:
    missing_text = ", ".join(missing_sides)
    with pytest.raises(PhosPyValidationError) as exc_info:
        validate_reference_context_compatibility(
            left,
            right,
            operation="unit unknown",
        )
    assert str(exc_info.value) == (
        "reference-context compatibility failed for operation='unit unknown': "
        f"unknown reference context for {missing_text}; configure an explicit "
        "workflow policy only when unknown reference context is scientifically "
        "acceptable"
    )

    warning = validate_reference_context_compatibility(
        left,
        right,
        operation="  unit unknown  ",
        allow_unknown=True,
    )

    assert isinstance(warning, ReferenceContextCompatibilityWarning)
    assert warning.code == REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE
    assert warning.missing_contexts == missing_sides
    assert warning.left_reference_context_id == left_id
    assert warning.right_reference_context_id == right_id
    noun = "contexts are" if len(missing_sides) > 1 else "context is"
    assert warning.to_payload() == {
        "code": "reference_context_unknown",
        "severity": "warning",
        "message": (
            "Reference-context compatibility could not be proven because "
            f"{missing_text} {noun} unknown for operation='unit unknown'."
        ),
        "operation": "unit unknown",
        "missing_contexts": list(missing_sides),
        "left_reference_context_id": left_id,
        "right_reference_context_id": right_id,
    }


@pytest.mark.parametrize("operation", ["", "   ", None, 1])
def test_reference_context_operation_validation(operation: object) -> None:
    with pytest.raises(
        PhosPyValidationError,
        match=(
            "^reference-context compatibility operation must be a non-empty string$"
        ),
    ):
        validate_reference_context_compatibility(
            _left_context(),
            _left_context(),
            operation=operation,  # type: ignore[arg-type]
        )


@pytest.mark.parametrize(
    "imports",
    [
        (
            "from phospy.phosphosites.site_keys import encode_site_key\n"
            "from phospy.references.context_compatibility import "
            "validate_reference_context_compatibility"
        ),
        (
            "from phospy.references.context_compatibility import "
            "validate_reference_context_compatibility\n"
            "from phospy.phosphosites.site_keys import encode_site_key"
        ),
    ],
)
def test_identity_and_reference_compatibility_import_orders_are_cycle_free(
    imports: str,
) -> None:
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(SRC_ROOT), environment.get("PYTHONPATH", "")]
    )
    completed = subprocess.run(
        [sys.executable, "-c", imports],
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr

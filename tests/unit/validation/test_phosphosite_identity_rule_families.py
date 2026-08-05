from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import PhosPyInputError
from phospy.errors.validation import PhosPyValidationError
from phospy.science.references.models import ReferenceContext
from phospy.science.sites.identity_rules.contracts import (
    ANALYSIS_READY_DATASET_IDENTITY_CONTRACT,
    REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE,
    SITE_SEQUENCE_COLUMN,
    PhosphositeIdentityContract,
    ReferenceContextCompatibilityWarning,
    SequenceContextRequirement,
)
from phospy.science.sites.identity_rules.dataset_identity import (
    enforce_analysis_ready_site_key_index,
    enforce_display_id_column,
    enforce_phosphosite_identity_contract,
    enforce_site_key_column,
    enforce_site_key_column_matches_index,
)
from phospy.science.sites.identity_rules.parsing import (
    looks_like_display_site_index,
    parse_row_site_token,
    parse_site_token,
    resolve_explicit_position,
    resolve_row_position,
    resolve_row_residue,
)
from phospy.science.sites.identity_rules.reference_context import (
    validate_reference_context_compatibility,
)
from phospy.science.sites.identity_rules.result_identity import (
    enforce_result_table_identity_contract,
)
from tests.support.site_keys import (
    protein_site_key,
    site_key_context_columns,
)


def _reference_context(**overrides: object) -> ReferenceContext:
    values = {
        "organism": "rat",
        "protein_namespace": "protein_id",
        "source_name": "unit-reference",
        "source_version": "v1",
        "proteome_version": "proteome-1",
        "reference_table_sha256": "a" * 64,
    }
    values.update(overrides)
    return ReferenceContext(**values)


def _site_key(
    protein_identifier: str = "P28482",
    site: str = "Y182",
) -> str:
    return protein_site_key(protein_identifier=protein_identifier, site=site)


def _dataset_identity_frame() -> pd.DataFrame:
    site_keys = pd.Index(
        [
            _site_key("P28482", "Y182"),
            _site_key("Q99999", "Y182"),
        ],
        name="site_key",
    )
    return pd.DataFrame(
        {
            "site_key": site_keys.astype(str).tolist(),
            "display_id": ["MAPK14;Y182;", "MAPK14;Y182;"],
            **site_key_context_columns(site_keys),
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["AAAAAAAYAAAAAAA", "AAAAAAAYAAAAAAA"],
        },
        index=site_keys.copy(),
    )


def _result_identity_table() -> pd.DataFrame:
    metadata = _dataset_identity_frame()
    return pd.DataFrame(
        {
            "site_key": metadata.loc[:, "site_key"].astype(str).tolist(),
            "display_id": metadata.loc[:, "display_id"].astype(str).tolist(),
            "organism": metadata.loc[:, "organism"].astype(str).tolist(),
            "protein_namespace": metadata.loc[:, "protein_namespace"]
            .astype(str)
            .tolist(),
            "protein_identifier": metadata.loc[:, "protein_identifier"]
            .astype(str)
            .tolist(),
            "gene_symbol": metadata.loc[:, "gene_symbol"].astype(str).tolist(),
            "site": metadata.loc[:, "site"].astype(str).tolist(),
            "logFC": [1.0, -1.0],
        },
        index=metadata.index.copy(),
    )


def _parsing_frame(**overrides: object) -> pd.DataFrame:
    values = {"site": "Y182", "residue": "Y", "position": 182}
    values.update(overrides)
    return pd.DataFrame(values, index=pd.Index(["row-1"], name="row_id"))


def test_reference_context_known_matching_contexts_pass() -> None:
    assert (
        validate_reference_context_compatibility(
            _reference_context(),
            _reference_context(),
            operation="unit matching",
        )
        is None
    )


def test_reference_context_known_mismatches_report_fields_and_summaries() -> None:
    with pytest.raises(PhosPyValidationError) as exc_info:
        validate_reference_context_compatibility(
            _reference_context(source_version="v1"),
            _reference_context(source_version="v2"),
            operation="unit mismatch",
        )

    message = str(exc_info.value)
    assert "operation='unit mismatch'" in message
    assert "mismatched_fields=source_version" in message
    assert "left={reference_context_id=" in message
    assert "right={reference_context_id=" in message


@pytest.mark.parametrize(
    ("left", "right", "missing_contexts"),
    [
        (_reference_context(), None, ("right",)),
        (None, _reference_context(), ("left",)),
        (None, None, ("left", "right")),
    ],
)
def test_reference_context_missing_sides_are_detected(
    left: ReferenceContext | None,
    right: ReferenceContext | None,
    missing_contexts: tuple[str, ...],
) -> None:
    with pytest.raises(PhosPyValidationError, match="unknown reference context"):
        validate_reference_context_compatibility(
            left,
            right,
            operation="unit missing disallowed",
        )

    warning = validate_reference_context_compatibility(
        left,
        right,
        operation="unit missing allowed",
        allow_unknown=True,
    )

    assert isinstance(warning, ReferenceContextCompatibilityWarning)
    assert warning.missing_contexts == missing_contexts


@pytest.mark.parametrize(
    ("allow_unknown", "should_warn"),
    [(False, False), (True, True)],
)
def test_reference_context_supported_unknown_policies(
    allow_unknown: bool,
    should_warn: bool,
) -> None:
    if not should_warn:
        with pytest.raises(PhosPyValidationError):
            validate_reference_context_compatibility(
                _reference_context(),
                None,
                operation="unit policy",
                allow_unknown=allow_unknown,
            )
        return

    warning = validate_reference_context_compatibility(
        _reference_context(),
        None,
        operation="unit policy",
        allow_unknown=allow_unknown,
    )
    assert isinstance(warning, ReferenceContextCompatibilityWarning)


def test_reference_context_warning_payload_contents() -> None:
    left = _reference_context()

    warning = validate_reference_context_compatibility(
        left,
        None,
        operation=" unit payload ",
        allow_unknown=True,
    )

    assert isinstance(warning, ReferenceContextCompatibilityWarning)
    payload = warning.to_payload()
    assert payload["code"] == REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE
    assert payload["severity"] == "warning"
    assert payload["operation"] == "unit payload"
    assert payload["missing_contexts"] == ["right"]
    assert payload["left_reference_context_id"] == left.reference_context_id
    assert payload["right_reference_context_id"] is None
    assert "right context is unknown" in str(payload["message"])


def test_dataset_identity_accepts_valid_site_key_index() -> None:
    frame = _dataset_identity_frame()

    observed = enforce_analysis_ready_site_key_index(
        frame.index,
        field_name="unit.site_metadata.index",
        error_type=PhosPyValidationError,
    )

    assert observed is frame.index


def test_dataset_identity_rejects_duplicate_site_keys() -> None:
    site_key = _site_key()
    duplicate_index = pd.Index([site_key, site_key], name="site_key")

    with pytest.raises(PhosPyValidationError, match="duplicate_site_key_values"):
        enforce_analysis_ready_site_key_index(
            duplicate_index,
            field_name="unit.site_metadata.index",
            error_type=PhosPyValidationError,
        )


def test_dataset_identity_rejects_missing_site_key_column() -> None:
    frame = _dataset_identity_frame().drop(columns=["site_key"])

    with pytest.raises(
        PhosPyValidationError, match="missing required columns: site_key"
    ):
        enforce_site_key_column(
            site_metadata=frame,
            field_name="unit.site_metadata",
            error_type=PhosPyValidationError,
        )


def test_dataset_identity_rejects_site_key_column_index_mismatch() -> None:
    frame = _dataset_identity_frame()
    frame.loc[frame.index[0], "site_key"] = _site_key("P31749", "T308")

    with pytest.raises(PhosPyValidationError, match="must exactly match"):
        enforce_site_key_column_matches_index(
            site_metadata=frame,
            field_name="unit.site_metadata",
            error_type=PhosPyValidationError,
        )


def test_dataset_identity_rejects_missing_display_id_column() -> None:
    frame = _dataset_identity_frame().drop(columns=["display_id"])

    with pytest.raises(
        PhosPyValidationError,
        match="missing required columns: display_id",
    ):
        enforce_display_id_column(
            site_metadata=frame,
            field_name="unit.site_metadata",
            error_type=PhosPyValidationError,
        )


def test_dataset_identity_rejects_invalid_display_id_values() -> None:
    frame = _dataset_identity_frame()
    frame.loc[frame.index[0], "display_id"] = "MAPK14 Y182"

    with pytest.raises(PhosPyValidationError, match="GENE;SITE"):
        enforce_display_id_column(
            site_metadata=frame,
            field_name="unit.site_metadata",
            error_type=PhosPyValidationError,
        )


def test_dataset_identity_requires_protein_identity_metadata() -> None:
    frame = _dataset_identity_frame().drop(columns=["protein_identifier"])

    with pytest.raises(
        PhosPyValidationError,
        match="missing required columns: protein_identifier",
    ):
        enforce_phosphosite_identity_contract(
            site_metadata=frame,
            field_name="unit.site_metadata",
            contract=ANALYSIS_READY_DATASET_IDENTITY_CONTRACT,
            error_type=PhosPyValidationError,
        )


def test_dataset_identity_accepts_repeated_display_ids_with_unique_site_keys() -> None:
    frame = _dataset_identity_frame()

    enforce_phosphosite_identity_contract(
        site_metadata=frame,
        field_name="unit.site_metadata",
        contract=ANALYSIS_READY_DATASET_IDENTITY_CONTRACT,
        error_type=PhosPyValidationError,
    )


def test_dataset_identity_accepts_opaque_site_values_where_explicitly_permitted() -> (
    None
):
    site_key = _site_key("P28482", "S1")
    frame = pd.DataFrame(
        {
            "site_key": [site_key],
            "display_id": ["MAPK14;FOO;"],
            "site": ["FOO"],
            SITE_SEQUENCE_COLUMN: ["AAAAAAASAAAAAAA"],
        },
        index=pd.Index([site_key], name="site_key"),
    )
    contract = PhosphositeIdentityContract(
        contract_id="unit_opaque_sequence_identity",
        required_columns=("site_key", "display_id", "site", SITE_SEQUENCE_COLUMN),
        require_protein_context=False,
        require_site_key_metadata_coherence=False,
        sequence_context=SequenceContextRequirement.CENTRED,
        sequence_context_contract=None,
    )

    enforce_phosphosite_identity_contract(
        site_metadata=frame,
        field_name="unit.site_metadata",
        contract=contract,
        error_type=PhosPyValidationError,
        allow_opaque_site_values=True,
    )


def test_result_identity_accepts_valid_result_metadata() -> None:
    enforce_result_table_identity_contract(
        table=_result_identity_table(),
        field_name="unit.result",
        error_type=PhosPyInputError,
    )


def test_result_identity_rejects_missing_result_site_keys() -> None:
    table = _result_identity_table().drop(columns=["site_key"])

    with pytest.raises(PhosPyInputError, match="missing required columns: site_key"):
        enforce_result_table_identity_contract(
            table=table,
            field_name="unit.result",
            error_type=PhosPyInputError,
        )


def test_result_identity_rejects_display_id_mismatch() -> None:
    table = _result_identity_table()
    table.loc[table.index[0], "display_id"] = "AKT1;Y182;"

    with pytest.raises(
        PhosPyInputError,
        match=r"display_id does not match gene_symbol \+ site",
    ):
        enforce_result_table_identity_contract(
            table=table,
            field_name="unit.result",
            error_type=PhosPyInputError,
        )


def test_result_identity_rejects_result_index_mismatch() -> None:
    table = _result_identity_table()
    mismatched_index = pd.Index(
        [_site_key("P31749", "T308"), table.index[1]],
        name="site_key",
    )
    table = table.set_axis(mismatched_index, axis="index")

    with pytest.raises(PhosPyInputError, match="must exactly match"):
        enforce_result_table_identity_contract(
            table=table,
            field_name="unit.result",
            error_type=PhosPyInputError,
        )


def test_result_identity_preserves_operation_specific_error_context() -> None:
    table = _result_identity_table()
    table.loc[table.index[0], "site"] = "Y999"

    with pytest.raises(PhosPyInputError) as exc_info:
        enforce_result_table_identity_contract(
            table=table,
            field_name="unit.operation.result",
            error_type=PhosPyInputError,
            context_label="Unit result identity metadata",
        )

    message = str(exc_info.value)
    assert message.startswith(
        "Unit result identity metadata is inconsistent with site_key"
    )
    assert "unit.operation.result" in message
    assert "site_key encodes Y182 but row metadata site is 'Y999'" in message


@pytest.mark.parametrize(
    ("value", "residue", "position"),
    [("S1", "S", 1), ("t308", "T", 308), (" y182 ", "Y", 182)],
)
def test_parsing_accepts_valid_residues_and_normalized_site_tokens(
    value: str,
    residue: str,
    position: int,
) -> None:
    parsed = parse_site_token(
        value,
        field_name="unit.site",
        error_type=PhosPyValidationError,
    )

    assert parsed.residue == residue
    assert parsed.position == position


@pytest.mark.parametrize("value", ["A123", "S0", "S", "123", "S-1"])
def test_parsing_rejects_malformed_residues_and_positions(value: object) -> None:
    with pytest.raises(PhosPyValidationError, match="strict 'S/T/Y<position>'"):
        parse_site_token(
            value,
            field_name="unit.site",
            error_type=PhosPyValidationError,
        )


@pytest.mark.parametrize("value", [0, -1, "182", 1.5, True])
def test_parsing_rejects_invalid_explicit_positions(value: object) -> None:
    frame = _parsing_frame(position=value)

    with pytest.raises(PhosPyValidationError):
        resolve_explicit_position(
            site_metadata=frame,
            row_id="row-1",
            field_name="unit.site_metadata",
            error_type=PhosPyValidationError,
        )


def test_parsing_resolves_row_residue_and_position_from_site_token() -> None:
    frame = _parsing_frame(residue=None).drop(columns=["position"])

    assert (
        resolve_row_residue(
            site_metadata=frame,
            row_id="row-1",
            field_name="unit.site_metadata",
            error_type=PhosPyValidationError,
        )
        == "Y"
    )
    assert (
        resolve_row_position(
            site_metadata=frame,
            row_id="row-1",
            field_name="unit.site_metadata",
            error_type=PhosPyValidationError,
        )
        == 182
    )


def test_parsing_rejects_ambiguous_or_inconsistent_tokens() -> None:
    frame = _parsing_frame(site="Y182/T183")

    with pytest.raises(PhosPyValidationError, match="strict 'S/T/Y<position>'"):
        parse_row_site_token(
            site_metadata=frame,
            row_id="row-1",
            field_name="unit.site_metadata",
            error_type=PhosPyValidationError,
        )

    inconsistent = _parsing_frame(site="Y182", residue="S")
    with pytest.raises(PhosPyValidationError, match="inconsistent residue metadata"):
        resolve_row_residue(
            site_metadata=inconsistent,
            row_id="row-1",
            field_name="unit.site_metadata",
            error_type=PhosPyValidationError,
        )


@pytest.mark.parametrize("missing_value", [None, pd.NA, ""])
def test_parsing_reports_missing_values_for_identity_derivation(
    missing_value: object,
) -> None:
    frame = _parsing_frame(site=missing_value, residue=None).drop(columns=["position"])

    with pytest.raises(PhosPyValidationError, match="requires residue metadata"):
        resolve_row_residue(
            site_metadata=frame,
            row_id="row-1",
            field_name="unit.site_metadata",
            error_type=PhosPyValidationError,
        )
    with pytest.raises(PhosPyValidationError, match="requires position metadata"):
        resolve_row_position(
            site_metadata=frame,
            row_id="row-1",
            field_name="unit.site_metadata",
            error_type=PhosPyValidationError,
        )


def test_parsing_recognizes_valid_display_site_indexes_only() -> None:
    assert looks_like_display_site_index(
        pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id")
    )
    assert not looks_like_display_site_index(
        pd.Index(["rat|protein_id|P28482|Y182"], name="site_key")
    )
    assert not looks_like_display_site_index(pd.Index(["MAPK14 Y182"], name="site_id"))

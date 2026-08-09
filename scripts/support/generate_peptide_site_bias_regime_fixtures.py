#!/usr/bin/env python3
"""Build independent peptide-to-site synthetic known-truth fixture payloads."""

from __future__ import annotations

from collections import defaultdict
from typing import Any

SAMPLE_COLUMNS: tuple[str, ...] = ("sample_1", "sample_2")
CLASSIFICATION = "synthetic_validation"
DECLARED_ESTIMAND = (
    "site-level sample abundance from fractional peptide-row allocation "
    "followed by a finite-value arithmetic mean per resolved site and sample"
)


def build_peptide_site_bias_fixture_payload() -> dict[str, Any]:
    """Return raw inputs and independently calculated expected outputs."""

    evidence_rows = _peptide_evidence_rows()
    mapping_rows = _site_mapping_rows()
    known_truth_rows = _known_truth_rows()
    expected_estimate_rows = _expected_site_estimate_rows(
        evidence_rows=evidence_rows,
        mapping_rows=mapping_rows,
    )
    bias_rows = _bias_rows(
        expected_estimate_rows=expected_estimate_rows,
        known_truth_rows=known_truth_rows,
        mapping_rows=mapping_rows,
        evidence_rows=evidence_rows,
    )
    return {
        "peptide_evidence_rows": evidence_rows,
        "peptide_evidence_columns": _peptide_evidence_columns(),
        "site_mapping_rows": mapping_rows,
        "site_mapping_columns": _site_mapping_columns(),
        "known_truth_rows": known_truth_rows,
        "known_truth_columns": _known_truth_columns(),
        "expected_site_estimate_rows": expected_estimate_rows,
        "expected_site_estimate_columns": _expected_site_estimate_columns(),
        "bias_rows": bias_rows,
        "bias_columns": _bias_columns(),
        "expected_bias_summary": _expected_bias_summary(bias_rows),
        "policy_assumptions": _policy_assumptions(),
    }


def _peptide_evidence_columns() -> tuple[str, ...]:
    return (
        "regime",
        "case_id",
        "peptide_row_id",
        "site_id",
        "unique_feature_id",
        "gene_symbol",
        "protein_accession",
        "site_string",
        "sample_1",
        "sample_2",
        "peptide_sequence",
        "modified_peptide_sequence",
        "multi_site",
        "provenance_source",
        "site_sequence",
        "localisation_confidence",
        "evidence_classification",
    )


def _site_mapping_columns() -> tuple[str, ...]:
    return (
        "regime",
        "case_id",
        "peptide_row_id",
        "site_id",
        "mapping_weight",
        "mapping_uncertainty",
        "multi_site_policy",
        "is_multi_site",
        "evidence_classification",
    )


def _known_truth_columns() -> tuple[str, ...]:
    return (
        "regime",
        "case_id",
        "site_id",
        "sample_id",
        "true_site_abundance",
        "truth_status",
        "truth_source",
        "evidence_classification",
    )


def _expected_site_estimate_columns() -> tuple[str, ...]:
    return (
        "regime",
        "case_id",
        "site_id",
        "sample_id",
        "expected_site_abundance",
        "estimator_output_status",
        "finite_allocated_observation_count",
        "declared_estimand",
        "evidence_classification",
    )


def _bias_columns() -> tuple[str, ...]:
    return (
        "regime",
        "case_id",
        "site_id",
        "sample_id",
        "true_site_abundance",
        "expected_site_abundance",
        "resolved_site_abundance",
        "signed_bias",
        "absolute_bias",
        "mapping_weight",
        "localisation_confidence",
        "bias_source",
        "estimator_output_status",
        "bias_status",
        "evidence_classification",
    )


def _peptide_evidence_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    rows.extend(
        [
            _evidence_row(
                regime="duplicate_concordant_control",
                case_id="duplicate_concordant_control",
                peptide_row_id="dup_concordant_a",
                site_id="DUPCON;S1;",
                gene_symbol="DUPCON",
                protein_accession="P_SYN_DUPCON",
                site_string="S1",
                sample_1=100.0,
                sample_2=100.0,
                peptide_sequence="AAASAAA",
                site_sequence="AAASAAA",
                localisation_confidence=0.99,
            ),
            _evidence_row(
                regime="duplicate_concordant_control",
                case_id="duplicate_concordant_control",
                peptide_row_id="dup_concordant_b",
                site_id="DUPCON;S1;",
                gene_symbol="DUPCON",
                protein_accession="P_SYN_DUPCON",
                site_string="S1",
                sample_1=100.0,
                sample_2=100.0,
                peptide_sequence="AAASAAA",
                site_sequence="AAASAAA",
                localisation_confidence=0.99,
            ),
            _evidence_row(
                regime="duplicate_discordant",
                case_id="duplicate_discordant",
                peptide_row_id="dup_discordant_low",
                site_id="DUPDIS;S1;",
                gene_symbol="DUPDIS",
                protein_accession="P_SYN_DUPDIS",
                site_string="S1",
                sample_1=100.0,
                sample_2=100.0,
                peptide_sequence="CCCSCCC",
                site_sequence="CCCSCCC",
                localisation_confidence=0.98,
            ),
            _evidence_row(
                regime="duplicate_discordant",
                case_id="duplicate_discordant",
                peptide_row_id="dup_discordant_high",
                site_id="DUPDIS;S1;",
                gene_symbol="DUPDIS",
                protein_accession="P_SYN_DUPDIS",
                site_string="S1",
                sample_1=130.0,
                sample_2=130.0,
                peptide_sequence="CCCSCCC",
                site_sequence="CCCSCCC",
                localisation_confidence=0.98,
            ),
        ]
    )
    rows.append(
        _evidence_row(
            regime="ambiguous_equal_split",
            case_id="ambiguous_equal_split",
            peptide_row_id="amb_equal",
            site_id="AMBEQ;S4;",
            gene_symbol="AMBEQ",
            protein_accession="P_SYN_AMBEQ",
            site_string="S4,T6",
            sample_1=100.0,
            sample_2=100.0,
            peptide_sequence="AAASATAAA",
            multi_site=True,
            localisation_confidence=0.70,
        )
    )
    rows.append(
        _evidence_row(
            regime="ambiguous_unequal_fraction",
            case_id="ambiguous_unequal_fraction",
            peptide_row_id="amb_unequal",
            site_id="AMBUNEQ;S4;",
            gene_symbol="AMBUNEQ",
            protein_accession="P_SYN_AMBUNEQ",
            site_string="S4,T6",
            sample_1=80.0,
            sample_2=80.0,
            peptide_sequence="CCCSCTCCC",
            multi_site=True,
            localisation_confidence=0.65,
        )
    )
    rows.extend(
        [
            _evidence_row(
                regime="missing_observation",
                case_id="missing_value_in_one_sample",
                peptide_row_id="miss_partial_observed",
                site_id="MISSPART;S1;",
                gene_symbol="MISSPART",
                protein_accession="P_SYN_MISSPART",
                site_string="S1",
                sample_1=100.0,
                sample_2=80.0,
                peptide_sequence="MMMSMMM",
                site_sequence="MMMSMMM",
                localisation_confidence=0.95,
            ),
            _evidence_row(
                regime="missing_observation",
                case_id="missing_value_in_one_sample",
                peptide_row_id="miss_partial_missing_sample",
                site_id="MISSPART;S1;",
                gene_symbol="MISSPART",
                protein_accession="P_SYN_MISSPART",
                site_string="S1",
                sample_1=80.0,
                sample_2=None,
                peptide_sequence="NNNSNNN",
                site_sequence="MMMSMMM",
                localisation_confidence=0.95,
            ),
            _evidence_row(
                regime="missing_observation",
                case_id="entirely_missing_peptide_observation",
                peptide_row_id="miss_entire",
                site_id="MISSENT;S1;",
                gene_symbol="MISSENT",
                protein_accession="P_SYN_MISSENT",
                site_string="S1",
                sample_1=None,
                sample_2=None,
                peptide_sequence="GGGSGGG",
                site_sequence="GGGSGGG",
                localisation_confidence=0.90,
            ),
        ]
    )
    rows.extend(
        [
            _evidence_row(
                regime="localisation_error",
                case_id="localisation_assigned_to_wrong_site",
                peptide_row_id="loc_true_zero_anchor",
                site_id="LOCERR;S1;",
                gene_symbol="LOCERR",
                protein_accession="P_SYN_LOCERR",
                site_string="S1",
                sample_1=0.0,
                sample_2=0.0,
                peptide_sequence="LLLSLLL",
                site_sequence="LLLSLLL",
                localisation_confidence=0.20,
            ),
            _evidence_row(
                regime="localisation_error",
                case_id="localisation_assigned_to_wrong_site",
                peptide_row_id="loc_false_assignment",
                site_id="LOCERR;T2;",
                gene_symbol="LOCERR",
                protein_accession="P_SYN_LOCERR",
                site_string="T2",
                sample_1=100.0,
                sample_2=100.0,
                peptide_sequence="LLLTLLL",
                site_sequence="LLLTLLL",
                localisation_confidence=0.80,
            ),
        ]
    )
    return rows


def _evidence_row(
    *,
    regime: str,
    case_id: str,
    peptide_row_id: str,
    site_id: str,
    gene_symbol: str,
    protein_accession: str,
    site_string: str,
    sample_1: float | None,
    sample_2: float | None,
    peptide_sequence: str,
    site_sequence: str = "",
    localisation_confidence: float,
    multi_site: bool = False,
) -> dict[str, object]:
    return {
        "regime": regime,
        "case_id": case_id,
        "peptide_row_id": peptide_row_id,
        "site_id": site_id,
        "unique_feature_id": f"feature_{peptide_row_id}",
        "gene_symbol": gene_symbol,
        "protein_accession": protein_accession,
        "site_string": site_string,
        "sample_1": sample_1,
        "sample_2": sample_2,
        "peptide_sequence": peptide_sequence,
        "modified_peptide_sequence": peptide_sequence.replace("S", "[pS]", 1)
        if "S" in peptide_sequence
        else peptide_sequence.replace("T", "[pT]", 1),
        "multi_site": "true" if multi_site else "false",
        "provenance_source": "synthetic-peptide-site-known-truth-v1",
        "site_sequence": site_sequence,
        "localisation_confidence": localisation_confidence,
        "evidence_classification": CLASSIFICATION,
    }


def _site_mapping_rows() -> list[dict[str, object]]:
    mapping: list[dict[str, object]] = []
    for row in _peptide_evidence_rows():
        peptide_row_id = str(row["peptide_row_id"])
        regime = str(row["regime"])
        case_id = str(row["case_id"])
        if peptide_row_id == "amb_equal":
            mapping.extend(
                [
                    _mapping_row(
                        regime=regime,
                        case_id=case_id,
                        peptide_row_id=peptide_row_id,
                        site_id="AMBEQ;S4;",
                        mapping_weight=0.5,
                        is_multi_site=True,
                    ),
                    _mapping_row(
                        regime=regime,
                        case_id=case_id,
                        peptide_row_id=peptide_row_id,
                        site_id="AMBEQ;T6;",
                        mapping_weight=0.5,
                        is_multi_site=True,
                    ),
                ]
            )
            continue
        if peptide_row_id == "amb_unequal":
            mapping.extend(
                [
                    _mapping_row(
                        regime=regime,
                        case_id=case_id,
                        peptide_row_id=peptide_row_id,
                        site_id="AMBUNEQ;S4;",
                        mapping_weight=0.25,
                        is_multi_site=True,
                    ),
                    _mapping_row(
                        regime=regime,
                        case_id=case_id,
                        peptide_row_id=peptide_row_id,
                        site_id="AMBUNEQ;T6;",
                        mapping_weight=0.75,
                        is_multi_site=True,
                    ),
                ]
            )
            continue
        mapping.append(
            _mapping_row(
                regime=regime,
                case_id=case_id,
                peptide_row_id=peptide_row_id,
                site_id=str(row["site_id"]),
                mapping_weight=1.0,
                is_multi_site=False,
            )
        )
    return mapping


def _mapping_row(
    *,
    regime: str,
    case_id: str,
    peptide_row_id: str,
    site_id: str,
    mapping_weight: float,
    is_multi_site: bool,
) -> dict[str, object]:
    return {
        "regime": regime,
        "case_id": case_id,
        "peptide_row_id": peptide_row_id,
        "site_id": site_id,
        "mapping_weight": mapping_weight,
        "mapping_uncertainty": "true" if is_multi_site else "false",
        "multi_site_policy": "split",
        "is_multi_site": "true" if is_multi_site else "false",
        "evidence_classification": CLASSIFICATION,
    }


def _known_truth_rows() -> list[dict[str, object]]:
    truth_specs = (
        (
            "duplicate_concordant_control",
            "duplicate_concordant_control",
            "DUPCON;S1;",
            (100.0, 100.0),
        ),
        (
            "duplicate_discordant",
            "duplicate_discordant",
            "DUPDIS;S1;",
            (100.0, 100.0),
        ),
        (
            "ambiguous_equal_split",
            "ambiguous_equal_split",
            "AMBEQ;S4;",
            (100.0, 100.0),
        ),
        (
            "ambiguous_equal_split",
            "ambiguous_equal_split",
            "AMBEQ;T6;",
            (0.0, 0.0),
        ),
        (
            "ambiguous_unequal_fraction",
            "ambiguous_unequal_fraction",
            "AMBUNEQ;S4;",
            (80.0, 80.0),
        ),
        (
            "ambiguous_unequal_fraction",
            "ambiguous_unequal_fraction",
            "AMBUNEQ;T6;",
            (0.0, 0.0),
        ),
        (
            "missing_observation",
            "missing_value_in_one_sample",
            "MISSPART;S1;",
            (90.0, 100.0),
        ),
        (
            "missing_observation",
            "entirely_missing_peptide_observation",
            "MISSENT;S1;",
            (100.0, 100.0),
        ),
        (
            "localisation_error",
            "localisation_assigned_to_wrong_site",
            "LOCERR;S1;",
            (100.0, 100.0),
        ),
        (
            "localisation_error",
            "localisation_assigned_to_wrong_site",
            "LOCERR;T2;",
            (0.0, 0.0),
        ),
    )
    rows: list[dict[str, object]] = []
    for regime, case_id, site_id, values in truth_specs:
        for sample_id, value in zip(SAMPLE_COLUMNS, values, strict=True):
            rows.append(
                {
                    "regime": regime,
                    "case_id": case_id,
                    "site_id": site_id,
                    "sample_id": sample_id,
                    "true_site_abundance": value,
                    "truth_status": "finite",
                    "truth_source": "declared_synthetic_known_truth",
                    "evidence_classification": CLASSIFICATION,
                }
            )
    return rows


def _expected_site_estimate_rows(
    *,
    evidence_rows: list[dict[str, object]],
    mapping_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    evidence_by_peptide = {str(row["peptide_row_id"]): row for row in evidence_rows}
    site_context: dict[str, tuple[str, str]] = {}
    allocated: dict[tuple[str, str], list[float]] = defaultdict(list)
    finite_counts: dict[tuple[str, str], int] = defaultdict(int)
    all_sites: set[str] = set()
    for mapping in mapping_rows:
        peptide_row_id = str(mapping["peptide_row_id"])
        site_id = str(mapping["site_id"])
        all_sites.add(site_id)
        site_context[site_id] = (str(mapping["regime"]), str(mapping["case_id"]))
        weight = float(mapping["mapping_weight"])
        evidence = evidence_by_peptide[peptide_row_id]
        for sample_id in SAMPLE_COLUMNS:
            raw_value = evidence[sample_id]
            if raw_value is None:
                continue
            allocated_value = float(raw_value) * weight
            allocated[(site_id, sample_id)].append(allocated_value)
            finite_counts[(site_id, sample_id)] += 1

    rows: list[dict[str, object]] = []
    for site_id in sorted(all_sites):
        regime, case_id = site_context[site_id]
        for sample_id in SAMPLE_COLUMNS:
            values = allocated.get((site_id, sample_id), [])
            expected_value = sum(values) / float(len(values)) if values else None
            rows.append(
                {
                    "regime": regime,
                    "case_id": case_id,
                    "site_id": site_id,
                    "sample_id": sample_id,
                    "expected_site_abundance": expected_value,
                    "estimator_output_status": (
                        "finite_estimate"
                        if expected_value is not None
                        else "missing_no_finite_allocated_evidence"
                    ),
                    "finite_allocated_observation_count": finite_counts.get(
                        (site_id, sample_id),
                        0,
                    ),
                    "declared_estimand": DECLARED_ESTIMAND,
                    "evidence_classification": CLASSIFICATION,
                }
            )
    return rows


def _bias_rows(
    *,
    expected_estimate_rows: list[dict[str, object]],
    known_truth_rows: list[dict[str, object]],
    mapping_rows: list[dict[str, object]],
    evidence_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    expected_by_cell = {
        (str(row["site_id"]), str(row["sample_id"])): row
        for row in expected_estimate_rows
    }
    mapping_weight_by_site = _single_mapping_weight_by_site(mapping_rows)
    localisation_by_site = _mean_localisation_by_site(
        evidence_rows=evidence_rows,
        mapping_rows=mapping_rows,
    )
    rows: list[dict[str, object]] = []
    for truth in known_truth_rows:
        site_id = str(truth["site_id"])
        sample_id = str(truth["sample_id"])
        expected = expected_by_cell[(site_id, sample_id)]
        expected_value = expected["expected_site_abundance"]
        true_value = float(truth["true_site_abundance"])
        signed_bias = (
            float(expected_value) - true_value if expected_value is not None else None
        )
        rows.append(
            {
                "regime": truth["regime"],
                "case_id": truth["case_id"],
                "site_id": site_id,
                "sample_id": sample_id,
                "true_site_abundance": true_value,
                "expected_site_abundance": expected_value,
                "resolved_site_abundance": expected_value,
                "signed_bias": signed_bias,
                "absolute_bias": abs(signed_bias) if signed_bias is not None else None,
                "mapping_weight": mapping_weight_by_site[site_id],
                "localisation_confidence": localisation_by_site[site_id],
                "bias_source": _bias_source(
                    str(truth["regime"]), str(truth["case_id"])
                ),
                "estimator_output_status": expected["estimator_output_status"],
                "bias_status": (
                    "finite"
                    if signed_bias is not None
                    else "nonestimable_missing_output"
                ),
                "evidence_classification": CLASSIFICATION,
            }
        )
    return rows


def _single_mapping_weight_by_site(
    mapping_rows: list[dict[str, object]],
) -> dict[str, float]:
    weights: dict[str, float] = {}
    for row in mapping_rows:
        site_id = str(row["site_id"])
        weights.setdefault(site_id, float(row["mapping_weight"]))
    return weights


def _mean_localisation_by_site(
    *,
    evidence_rows: list[dict[str, object]],
    mapping_rows: list[dict[str, object]],
) -> dict[str, float]:
    evidence_by_peptide = {str(row["peptide_row_id"]): row for row in evidence_rows}
    values_by_site: dict[str, list[float]] = defaultdict(list)
    for mapping in mapping_rows:
        site_id = str(mapping["site_id"])
        evidence = evidence_by_peptide[str(mapping["peptide_row_id"])]
        values_by_site[site_id].append(float(evidence["localisation_confidence"]))
    return {
        site_id: sum(values) / float(len(values))
        for site_id, values in values_by_site.items()
    }


def _bias_source(regime: str, case_id: str) -> str:
    if regime == "duplicate_concordant_control":
        return "duplicate peptide mean control"
    if regime == "duplicate_discordant":
        return "retained duplicate peptide rows with discordant abundances"
    if regime.startswith("ambiguous"):
        return "ambiguous multi-site fractional allocation"
    if case_id == "entirely_missing_peptide_observation":
        return "no finite peptide observation for the site/sample"
    if regime == "missing_observation":
        return "finite-value mean after sample-level missing peptide evidence"
    if regime == "localisation_error":
        return "signal assigned to the wrong resolved site"
    return regime


def _expected_bias_summary(bias_rows: list[dict[str, object]]) -> dict[str, object]:
    regimes = sorted({str(row["regime"]) for row in bias_rows})
    regime_payload: dict[str, object] = {}
    for regime in regimes:
        rows = [row for row in bias_rows if row["regime"] == regime]
        finite_abs = [
            float(row["absolute_bias"])
            for row in rows
            if row["absolute_bias"] is not None
        ]
        signed = [
            float(row["signed_bias"]) for row in rows if row["signed_bias"] is not None
        ]
        missing_count = len(rows) - len(finite_abs)
        regime_payload[regime] = {
            "finite_bias_cell_count": len(finite_abs),
            "missing_estimate_count": missing_count,
            "mean_absolute_bias": (
                sum(finite_abs) / float(len(finite_abs)) if finite_abs else None
            ),
            "maximum_absolute_bias": max(finite_abs) if finite_abs else None,
            "signed_bias_values": sorted(set(signed)),
            "case_ids": sorted({str(row["case_id"]) for row in rows}),
            "supported_interpretation": _supported_interpretation(regime),
        }
    return {
        "classification": CLASSIFICATION,
        "evidence_category": CLASSIFICATION,
        "validation_type": "production-path synthetic known-truth validation",
        "estimand": DECLARED_ESTIMAND,
        "bias_units": "linear_abundance",
        "absolute_tolerance": 1.0e-12,
        "regimes": regime_payload,
        "tuning_policy": (
            "fixture quantifies adverse-design sensitivity and must not be used "
            "to tune production allocation parameters"
        ),
    }


def _supported_interpretation(regime: str) -> str:
    if regime == "duplicate_concordant_control":
        return "concordant duplicate rows leave the retained-row mean unbiased"
    if regime == "duplicate_discordant":
        return (
            "duplicate peptide evidence can bias a site mean when duplicate "
            "peptide observations are discordant"
        )
    if regime == "ambiguous_equal_split":
        return (
            "equal splitting attenuates the true site and creates spurious "
            "support for the alternative site"
        )
    if regime == "ambiguous_unequal_fraction":
        return (
            "unequal explicit fractions allocate different signal magnitudes "
            "while retaining the same allocated-signal mean estimand"
        )
    if regime == "missing_observation":
        return (
            "missing peptide evidence yields finite means from the observed "
            "subset and missing estimates when no finite allocated value exists"
        )
    if regime == "localisation_error":
        return (
            "wrong localisation transfers signal away from the synthetic true "
            "site and into the assigned false site"
        )
    return regime


def _policy_assumptions() -> dict[str, object]:
    return {
        "classification": CLASSIFICATION,
        "evidence_category": CLASSIFICATION,
        "external_reference": None,
        "validation_type": "production-path synthetic known-truth validation",
        "tested_estimator": (
            "peptide_to_site_linear_abundance_fractional_allocation_arithmetic_mean_v1"
        ),
        "input_intensity_scale": "linear",
        "input_quantitative_meaning": "peptide_abundance",
        "output_intensity_scale": "linear",
        "output_quantitative_meaning": "phosphosite_abundance",
        "allocation_domain": "linear_abundance",
        "absolute_tolerance": 1.0e-12,
        "sample_columns": list(SAMPLE_COLUMNS),
        "multi_site_policy": "split",
        "mapping_weight_source": "explicit_mapping_weight",
        "mapping_weight_normalization": "sum_to_one_per_peptide_evidence_row",
        "signal_allocation": "peptide_signal_multiplied_by_mapping_fraction",
        "site_summarisation": "finite_value_arithmetic_mean_per_site_sample",
        "duplicate_evidence_policy": (
            "retain_duplicate_peptide_evidence_rows_as_separate_observations"
        ),
        "missing_value_policy": (
            "mean_finite_allocated_values_per_site_sample_preserve_missing_if_none"
        ),
        "localisation_summary_policy": (
            "descriptive_mean_of_finite_reported_localisation_confidence_values"
        ),
        "truth_source": "declared synthetic known truth",
        "expected_output_generation": (
            "independent standard-library arithmetic over declared raw rows and "
            "explicit mapping weights"
        ),
        "limitations": [
            "synthetic validation only",
            "not external parity",
            "not empirical validation",
            "not biological ground truth",
            "does not validate post-hoc peptide differential estimate combination",
        ],
    }


if __name__ == "__main__":
    raise SystemExit(
        "Import build_peptide_site_bias_fixture_payload() from this module via "
        "the release-validation fixture generator."
    )

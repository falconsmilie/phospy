#!/usr/bin/env python3
"""Generate compact release-validation regression fixtures."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_DIR = ROOT / "tests" / "fixtures" / "release_validation_regression"
DEFAULT_TIMESTAMP = "2026-07-24T00:00:00Z"
DEFAULT_SEED = 20260724
CANONICAL_TEXT_ENCODING = "utf-8"
CANONICAL_TEXT_NEWLINE = "\n"
CANONICAL_TEXT_BYTE_POLICY = "utf-8 LF with final newline"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Regenerate compact PhosPy regression fixtures for release validation."
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory where release-validation regression fixtures are written.",
    )
    parser.add_argument(
        "--timestamp",
        default=DEFAULT_TIMESTAMP,
        help="Reproducible generation timestamp recorded in manifests.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help="Recorded deterministic seed for fixture governance.",
    )
    parser.add_argument(
        "--manifest-outdir-label",
        default=None,
        help=argparse.SUPPRESS,
    )
    return parser.parse_args()


def _canonical_text_bytes(text: str) -> bytes:
    normalized = text.replace("\r\n", CANONICAL_TEXT_NEWLINE).replace(
        "\r", CANONICAL_TEXT_NEWLINE
    )
    if not normalized.endswith(CANONICAL_TEXT_NEWLINE):
        normalized += CANONICAL_TEXT_NEWLINE
    return normalized.encode(CANONICAL_TEXT_ENCODING)


def _write_canonical_text(path: Path, text: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = _canonical_text_bytes(text)
    path.write_bytes(payload)
    return hashlib.sha256(payload).hexdigest()


def _write_csv(
    path: Path, rows: list[dict[str, object]], columns: tuple[str, ...]
) -> str:
    handle = io.StringIO(newline="")
    writer = csv.DictWriter(
        handle,
        fieldnames=list(columns),
        lineterminator=CANONICAL_TEXT_NEWLINE,
    )
    writer.writeheader()
    writer.writerows(rows)
    return _write_canonical_text(path, handle.getvalue())


def _write_json(path: Path, payload: dict[str, Any]) -> str:
    return _write_canonical_text(
        path,
        json.dumps(payload, indent=2, sort_keys=False),
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _script_sha256() -> str:
    return _sha256(Path(__file__).resolve())


def _write_manifest(
    *,
    directory: Path,
    fixture_family: str,
    classification: str,
    generator_command: str,
    seed: int,
    timestamp: str,
    files: tuple[Path, ...],
    notes: str,
    source_policy: str | None = None,
    extra_fields: dict[str, Any] | None = None,
) -> None:
    payload: dict[str, Any] = {
        "manifest_schema_version": "fixture-manifest-v1",
        "fixture_family": fixture_family,
        "classification": classification,
        "generator": "scripts/active/generate_release_validation_regression_fixtures.py",
        "generator_sha256": _script_sha256(),
        "command": generator_command,
        "seed": int(seed),
        "generation_timestamp_utc": timestamp,
        "source_policy": source_policy
        or "deterministic PhosPy regression contract fixture; not external parity",
        "byte_policy": CANONICAL_TEXT_BYTE_POLICY,
        "software_versions": {"python": "runtime-independent"},
        "notes": notes,
        "files": [
            {
                "relative_path": path.relative_to(directory).as_posix(),
                "sha256": _sha256(path),
            }
            for path in files
        ],
    }
    if extra_fields:
        payload.update(extra_fields)
    _write_json(directory / "MANIFEST.json", payload)


def _evidence_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = (
        "peptide_row_id",
        "site_id",
        "unique_feature_id",
        "gene_symbol",
        "protein_accession",
        "site_string",
        "sample_a",
        "sample_b",
        "peptide_sequence",
        "modified_peptide_sequence",
        "multi_site",
        "provenance_source",
        "site_sequence",
        "localisation_confidence",
    )
    rows: list[dict[str, object]] = [
        {
            "peptide_row_id": "pep_equal",
            "site_id": "MAPK1;S10;",
            "unique_feature_id": "feat_equal",
            "gene_symbol": "MAPK1",
            "protein_accession": "P28482",
            "site_string": "S10,T12",
            "sample_a": 10.0,
            "sample_b": 12.0,
            "peptide_sequence": "AAASTTAAAA",
            "modified_peptide_sequence": "AAA[+80]STTAAAA",
            "multi_site": "true",
            "provenance_source": "release-validation-regression",
            "site_sequence": "",
            "localisation_confidence": 0.90,
        },
        {
            "peptide_row_id": "pep_unequal",
            "site_id": "MAPK1;S10;",
            "unique_feature_id": "feat_unequal",
            "gene_symbol": "MAPK1",
            "protein_accession": "P28482",
            "site_string": "S10,T12",
            "sample_a": 20.0,
            "sample_b": 24.0,
            "peptide_sequence": "CCCSATCCCC",
            "modified_peptide_sequence": "CCC[+80]SATCCCC",
            "multi_site": "true",
            "provenance_source": "release-validation-regression",
            "site_sequence": "",
            "localisation_confidence": 0.92,
        },
        {
            "peptide_row_id": "pep_single_low",
            "site_id": "MAPK1;S10;",
            "unique_feature_id": "feat_single_low",
            "gene_symbol": "MAPK1",
            "protein_accession": "P28482",
            "site_string": "S10",
            "sample_a": 30.0,
            "sample_b": 32.0,
            "peptide_sequence": "GGGSGGG",
            "modified_peptide_sequence": "GGG[+80]SGGG",
            "multi_site": "false",
            "provenance_source": "release-validation-regression",
            "site_sequence": "AAASAAA",
            "localisation_confidence": 0.95,
        },
        {
            "peptide_row_id": "pep_single_high",
            "site_id": "MAPK1;S10;",
            "unique_feature_id": "feat_single_high",
            "gene_symbol": "MAPK1",
            "protein_accession": "P28482",
            "site_string": "S10",
            "sample_a": 50.0,
            "sample_b": 52.0,
            "peptide_sequence": "HHHSHHH",
            "modified_peptide_sequence": "HHH[+80]SHHH",
            "multi_site": "false",
            "provenance_source": "release-validation-regression",
            "site_sequence": "AAASAAA",
            "localisation_confidence": 0.97,
        },
    ]
    return rows, columns


def _evidence_mapping_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = ("peptide_row_id", "site_id", "mapping_weight", "mapping_uncertainty")
    rows: list[dict[str, object]] = [
        {
            "peptide_row_id": "pep_equal",
            "site_id": "MAPK1;S10;",
            "mapping_weight": 0.5,
            "mapping_uncertainty": "true",
        },
        {
            "peptide_row_id": "pep_equal",
            "site_id": "MAPK1;T12;",
            "mapping_weight": 0.5,
            "mapping_uncertainty": "true",
        },
        {
            "peptide_row_id": "pep_unequal",
            "site_id": "MAPK1;S10;",
            "mapping_weight": 0.25,
            "mapping_uncertainty": "true",
        },
        {
            "peptide_row_id": "pep_unequal",
            "site_id": "MAPK1;T12;",
            "mapping_weight": 0.75,
            "mapping_uncertainty": "true",
        },
        {
            "peptide_row_id": "pep_single_low",
            "site_id": "MAPK1;S10;",
            "mapping_weight": 1.0,
            "mapping_uncertainty": "false",
        },
        {
            "peptide_row_id": "pep_single_high",
            "site_id": "MAPK1;S10;",
            "mapping_weight": 1.0,
            "mapping_uncertainty": "false",
        },
    ]
    return rows, columns


def _expected_evidence_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = ("site_id", "sample_a", "sample_b")
    rows: list[dict[str, object]] = [
        {"site_id": "MAPK1;S10;", "sample_a": 22.5, "sample_b": 24.0},
        {"site_id": "MAPK1;T12;", "sample_a": 10.0, "sample_b": 12.0},
    ]
    return rows, columns


def _sequence_fixture_rows(
    *,
    invalid_second_sequence: bool,
) -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = (
        "peptide_row_id",
        "site_id",
        "unique_feature_id",
        "gene_symbol",
        "protein_accession",
        "site_string",
        "sample_a",
        "sample_b",
        "peptide_sequence",
        "modified_peptide_sequence",
        "multi_site",
        "provenance_source",
        "site_sequence",
        "localisation_confidence",
    )
    second_sequence = "AAATAAA" if invalid_second_sequence else "CCCSCCC"
    rows: list[dict[str, object]] = []
    for index, sequence in enumerate(("AAASAAA", second_sequence), start=1):
        rows.append(
            {
                "peptide_row_id": f"seq_{index}",
                "site_id": "AKT1;S473;",
                "unique_feature_id": f"feat_seq_{index}",
                "gene_symbol": "AKT1",
                "protein_accession": "P31749",
                "site_string": "S473",
                "sample_a": float(10 * index),
                "sample_b": float(12 * index),
                "peptide_sequence": f"PEPTIDE{index}",
                "modified_peptide_sequence": f"PEP[+80]TIDE{index}",
                "multi_site": "false",
                "provenance_source": "release-validation-regression",
                "site_sequence": sequence,
                "localisation_confidence": 0.9 + 0.01 * index,
            }
        )
    return rows, columns


def generate_evidence_resolution_fixtures(
    *,
    outdir: Path,
    seed: int,
    timestamp: str,
    command: str,
) -> None:
    family_dir = outdir / "evidence_resolution"
    rows, columns = _evidence_rows()
    evidence_path = family_dir / "peptide_evidence.csv"
    _write_csv(evidence_path, rows, columns)
    rows, columns = _evidence_mapping_rows()
    mapping_path = family_dir / "site_mapping.csv"
    _write_csv(mapping_path, rows, columns)
    rows, columns = _expected_evidence_rows()
    expected_path = family_dir / "expected_split_phospho.csv"
    _write_csv(expected_path, rows, columns)
    rows, columns = _sequence_fixture_rows(invalid_second_sequence=False)
    conflict_path = family_dir / "sequence_conflict.csv"
    _write_csv(conflict_path, rows, columns)
    rows, columns = _sequence_fixture_rows(invalid_second_sequence=True)
    mixed_path = family_dir / "mixed_sequence_contexts.csv"
    _write_csv(mixed_path, rows, columns)

    _write_manifest(
        directory=family_dir,
        fixture_family="evidence_resolution",
        classification="regression",
        generator_command=command,
        seed=seed,
        timestamp=timestamp,
        files=(
            evidence_path,
            mapping_path,
            expected_path,
            conflict_path,
            mixed_path,
        ),
        notes=(
            "Peptide-to-site resolution regression cases for equal and unequal "
            "multi-site fractions, multiple peptides per site, sequence conflict, "
            "mixed valid/invalid sequence contexts, and row-order permutation."
        ),
    )


def _site_sequence(site_id: str) -> str:
    residue = site_id.split(";")[1][0].upper()
    return ("A" * 15) + residue + ("A" * 15)


def _kinase_display_ids() -> list[str]:
    return [f"KIN{index};S{index};" for index in range(1, 7)]


def _kinase_phospho_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = ("site_id", "sample_1", "sample_2", "sample_3", "sample_4")
    rows: list[dict[str, object]] = [
        {
            "site_id": site_id,
            "sample_1": sample_1,
            "sample_2": sample_2,
            "sample_3": sample_3,
            "sample_4": sample_4,
        }
        for site_id, sample_1, sample_2, sample_3, sample_4 in zip(
            _kinase_display_ids(),
            (8.0, 7.0, 2.0, 1.0, 3.0, 4.0),
            (7.5, 6.5, 2.2, 1.1, 3.2, 4.1),
            (1.0, 2.0, 7.0, 8.0, 3.0, 4.0),
            (1.2, 2.1, 6.8, 7.9, 3.1, 4.2),
            strict=True,
        )
    ]
    return rows, columns


def _kinase_site_metadata_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = (
        "site_id",
        "gene_symbol",
        "protein_id",
        "site",
        "site_sequence",
        "localisation_confidence",
    )
    rows: list[dict[str, object]] = []
    for index, site_id in enumerate(_kinase_display_ids(), start=1):
        rows.append(
            {
                "site_id": site_id,
                "gene_symbol": site_id.split(";", 1)[0],
                "protein_id": site_id.split(";", 1)[0],
                "site": site_id.split(";")[1],
                "site_sequence": _site_sequence(site_id),
                "localisation_confidence": 0.20 if index == 5 else 0.95,
            }
        )
    return rows, columns


def _kinase_substrate_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = ("kinase", "substrate_site")
    rows: list[dict[str, object]] = [
        {"kinase": "K_BELOW", "substrate_site": "KIN1;S1;"},
        {"kinase": "K_AT", "substrate_site": "KIN1;S1;"},
        {"kinase": "K_AT", "substrate_site": "KIN2;S2;"},
        {"kinase": "K_ABOVE", "substrate_site": "KIN1;S1;"},
        {"kinase": "K_ABOVE", "substrate_site": "KIN2;S2;"},
        {"kinase": "K_ABOVE", "substrate_site": "KIN3;S3;"},
        {"kinase": "K_SPARSE", "substrate_site": "KIN1;S1;"},
        {"kinase": "K_SPARSE", "substrate_site": "OFFSITE1;S1;"},
        {"kinase": "K_SPARSE", "substrate_site": "OFFSITE2;S2;"},
        {"kinase": "K_LOCALISATION", "substrate_site": "KIN5;S5;"},
        {"kinase": "K_LOCALISATION", "substrate_site": "KIN6;S6;"},
        {"kinase": "K_SSGSEA_HIGH", "substrate_site": "KIN1;S1;"},
        {"kinase": "K_SSGSEA_HIGH", "substrate_site": "KIN2;S2;"},
        {"kinase": "K_SSGSEA_LOW", "substrate_site": "KIN3;S3;"},
        {"kinase": "K_SSGSEA_LOW", "substrate_site": "KIN4;S4;"},
    ]
    return rows, columns


def _kinase_sequence_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = ("site_id", "site_sequence")
    site_ids = [
        *_kinase_display_ids(),
        "OFFSITE1;S1;",
        "OFFSITE2;S2;",
    ]
    rows = [
        {"site_id": site_id, "site_sequence": _site_sequence(site_id)}
        for site_id in site_ids
    ]
    return rows, columns


def _kinase_ssgsea_effect_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = ("site_id", "condition_positive", "condition_negative")
    rows = [
        {"site_id": "KIN1;S1;", "condition_positive": 4.0, "condition_negative": 1.0},
        {"site_id": "KIN2;S2;", "condition_positive": 3.0, "condition_negative": 2.0},
        {"site_id": "KIN3;S3;", "condition_positive": 2.0, "condition_negative": 3.0},
        {"site_id": "KIN4;S4;", "condition_positive": 1.0, "condition_negative": 4.0},
    ]
    return rows, columns


def _kinase_expected_contracts() -> dict[str, object]:
    return {
        "classification": "regression",
        "minimum_scoring_substrates": 2,
        "substrate_support_classes": {
            "below_minimum": ["K_BELOW", "K_SPARSE"],
            "at_minimum": ["K_AT"],
            "above_minimum": ["K_ABOVE"],
        },
        "ssgsea_expected_support": {
            "K_SSGSEA_HIGH": 2,
            "K_SSGSEA_LOW": 2,
        },
        "localisation_attrition_site": "KIN5;S5;",
        "external_reference": None,
    }


def generate_kinase_sparse_support_fixtures(
    *,
    outdir: Path,
    seed: int,
    timestamp: str,
    command: str,
) -> None:
    family_dir = outdir / "kinase_sparse_support"
    rows, columns = _kinase_phospho_rows()
    phospho_path = family_dir / "phospho.csv"
    _write_csv(phospho_path, rows, columns)
    rows, columns = _kinase_site_metadata_rows()
    metadata_path = family_dir / "site_metadata.csv"
    _write_csv(metadata_path, rows, columns)
    rows, columns = _kinase_substrate_rows()
    substrate_path = family_dir / "substrate_map.csv"
    _write_csv(substrate_path, rows, columns)
    rows, columns = _kinase_sequence_rows()
    sequence_path = family_dir / "site_sequences.csv"
    _write_csv(sequence_path, rows, columns)
    rows, columns = _kinase_ssgsea_effect_rows()
    effect_path = family_dir / "ssgsea_effect_matrix.csv"
    _write_csv(effect_path, rows, columns)
    expected_path = family_dir / "expected_contracts.json"
    _write_json(expected_path, _kinase_expected_contracts())

    _write_manifest(
        directory=family_dir,
        fixture_family="kinase_sparse_support",
        classification="regression",
        generator_command=command,
        seed=seed,
        timestamp=timestamp,
        files=(
            phospho_path,
            metadata_path,
            substrate_path,
            sequence_path,
            effect_path,
            expected_path,
        ),
        notes=(
            "Kinase sparse-support regression inputs covering below/at/above "
            "substrate minima, sparse reference overlap, leave-one-out versus "
            "self-including profile policy, localisation attrition, production "
            "threshold failure, and ssGSEA activity with/without permutation "
            "significance."
        ),
    )


def _signalome_network_rows(
    *,
    observation_count: int,
) -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = ("site_id", "K1", "K2")
    rows = [
        {"site_id": f"S{index}", "K1": float(index), "K2": float(index)}
        for index in range(1, observation_count + 1)
    ]
    return rows, columns


def _signalome_clustering_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = ("site_id", "K1", "K2", "K3_all_missing", "K4_partial_missing")
    rows = [
        {
            "site_id": "S1",
            "K1": 1.0,
            "K2": 2.0,
            "K3_all_missing": "",
            "K4_partial_missing": 10.0,
        },
        {
            "site_id": "S2",
            "K1": 1.2,
            "K2": 2.2,
            "K3_all_missing": "",
            "K4_partial_missing": "",
        },
        {
            "site_id": "S3",
            "K1": 0.9,
            "K2": 1.8,
            "K3_all_missing": "",
            "K4_partial_missing": 12.0,
        },
        {
            "site_id": "S4",
            "K1": -1.0,
            "K2": -2.0,
            "K3_all_missing": "",
            "K4_partial_missing": -10.0,
        },
        {
            "site_id": "S5",
            "K1": -1.2,
            "K2": -2.2,
            "K3_all_missing": "",
            "K4_partial_missing": "",
        },
        {
            "site_id": "S6",
            "K1": -0.9,
            "K2": -1.8,
            "K3_all_missing": "",
            "K4_partial_missing": -12.0,
        },
    ]
    return rows, columns


def _historical_threshold2_snapshot() -> dict[str, object]:
    return {
        "signalome_config": {
            "scientific": {
                "substrate_support_cutoff": 0.5,
                "assignment_policy": "cutoff_binary",
            },
            "clustering": {
                "module_count": None,
                "module_selection_primary_correlation_threshold": 0.5,
                "module_selection_fallback_correlation_threshold": 0.1,
                "module_selection_max_clusters": 10,
                "candidate_scoring_policy": "full",
                "clustering_engine": "exact_python",
            },
            "validation": {
                "score_preconditioning_policy": "error_on_drop",
                "allow_mixed_total_protein_quantitative_meaning": False,
                "reference_context_compatibility_policy": "require_known_match",
            },
            "output": {
                "network_correlation_threshold": 0.5,
                "network_policy": "signed",
                "network_min_paired_finite_observations": 2,
                "network_min_paired_finite_observations_effective": 2,
            },
            "performance": {
                "max_exact_tree_sites": 2000,
                "max_full_candidate_scoring_sites": 2000,
            },
        }
    }


def generate_signalome_safety_fixtures(
    *,
    outdir: Path,
    seed: int,
    timestamp: str,
    command: str,
) -> None:
    family_dir = outdir / "signalome_safety"
    rows, columns = _signalome_network_rows(observation_count=2)
    network_two_path = family_dir / "network_two_observations.csv"
    _write_csv(network_two_path, rows, columns)
    rows, columns = _signalome_network_rows(observation_count=3)
    network_three_path = family_dir / "network_three_observations.csv"
    _write_csv(network_three_path, rows, columns)
    rows, columns = _signalome_network_rows(observation_count=5)
    network_five_path = family_dir / "network_five_observations.csv"
    _write_csv(network_five_path, rows, columns)
    rows, columns = _signalome_clustering_rows()
    clustering_path = family_dir / "clustering_missing_dimensions.csv"
    _write_csv(clustering_path, rows, columns)
    snapshot_path = family_dir / "historical_threshold2_config.json"
    _write_json(snapshot_path, _historical_threshold2_snapshot())

    _write_manifest(
        directory=family_dir,
        fixture_family="signalome_safety",
        classification="regression",
        generator_command=command,
        seed=seed,
        timestamp=timestamp,
        files=(
            network_two_path,
            network_three_path,
            network_five_path,
            clustering_path,
            snapshot_path,
        ),
        notes=(
            "Signalome safety regression inputs covering paired-observation "
            "network thresholds, all-missing dimension dropping, partial-missing "
            "median imputation, clustering invariance after all-missing dimension "
            "addition, and historical threshold-2 bundle config reconstruction "
            "without allowing new threshold-2 execution."
        ),
    )


def _sps_ruv_phospho_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = ("site_id", "s1", "s2", "s3", "s4", "s5", "s6")
    planted_factor = (-2.0, -2.0, 0.0, 0.0, 2.0, 2.0)
    condition_indicator = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    rows: list[dict[str, object]] = []
    definitions = (
        ("control_1", 10.0, 1.5, 0.0),
        ("control_2", 20.0, -1.0, 0.0),
        ("control_3", 30.0, 0.5, 0.0),
        ("signal_protected", 40.0, 2.0, 5.0),
        ("signal_null", 50.0, 3.0, 0.0),
    )
    for site_id, baseline, factor_loading, condition_effect in definitions:
        values = [
            baseline
            + factor_loading * factor_value
            + condition_effect * condition_value
            for factor_value, condition_value in zip(
                planted_factor,
                condition_indicator,
                strict=True,
            )
        ]
        rows.append(
            {
                "site_id": site_id,
                **{f"s{position}": value for position, value in enumerate(values, 1)},
            }
        )
    return rows, columns


def _sps_ruv_sample_metadata_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = ("sample_id", "batch", "condition", "replicate", "planted_factor")
    rows = [
        {
            "sample_id": "s1",
            "batch": "run_1",
            "condition": "A",
            "replicate": "r1",
            "planted_factor": -2.0,
        },
        {
            "sample_id": "s2",
            "batch": "run_1",
            "condition": "B",
            "replicate": "r2",
            "planted_factor": -2.0,
        },
        {
            "sample_id": "s3",
            "batch": "run_2",
            "condition": "A",
            "replicate": "r3",
            "planted_factor": 0.0,
        },
        {
            "sample_id": "s4",
            "batch": "run_2",
            "condition": "B",
            "replicate": "r4",
            "planted_factor": 0.0,
        },
        {
            "sample_id": "s5",
            "batch": "run_3",
            "condition": "A",
            "replicate": "r5",
            "planted_factor": 2.0,
        },
        {
            "sample_id": "s6",
            "batch": "run_3",
            "condition": "B",
            "replicate": "r6",
            "planted_factor": 2.0,
        },
    ]
    return rows, columns


def _sps_ruv_known_truth() -> dict[str, object]:
    return {
        "classification": "synthetic_validation",
        "external_reference": None,
        "planted_unwanted_factor_by_sample": {
            "s1": -2.0,
            "s2": -2.0,
            "s3": 0.0,
            "s4": 0.0,
            "s5": 2.0,
            "s6": 2.0,
        },
        "control_site_ids": ["control_1", "control_2", "control_3"],
        "protected_signal_site_id": "signal_protected",
        "protected_condition_effect": {
            "contrast": "B_minus_A",
            "expected_difference": 5.0,
            "acceptance_absolute_tolerance": 1.0e-10,
        },
        "unwanted_factor_recovery": {
            "minimum_abs_correlation": 0.999,
            "control_row_max_span_after_correction": 1.0e-10,
        },
        "tuning_policy": (
            "fixture is release validation evidence only; implementation "
            "parameters must not be tuned against this fixture"
        ),
    }


def generate_sps_ruv_planted_unwanted_factor_fixtures(
    *,
    outdir: Path,
    seed: int,
    timestamp: str,
    command: str,
) -> None:
    family_dir = outdir / "sps_ruv_planted_unwanted_factor"
    rows, columns = _sps_ruv_phospho_rows()
    phospho_path = family_dir / "phospho.csv"
    _write_csv(phospho_path, rows, columns)
    rows, columns = _sps_ruv_sample_metadata_rows()
    metadata_path = family_dir / "sample_metadata.csv"
    _write_csv(metadata_path, rows, columns)
    truth_path = family_dir / "known_truth.json"
    _write_json(truth_path, _sps_ruv_known_truth())

    _write_manifest(
        directory=family_dir,
        fixture_family="sps_ruv_planted_unwanted_factor",
        classification="synthetic_validation",
        generator_command=command,
        seed=seed,
        timestamp=timestamp,
        files=(phospho_path, metadata_path, truth_path),
        notes=(
            "Independent synthetic known-truth fixture for native PhosPy "
            "SPS/RUV-style correction. It plants a one-dimensional unwanted "
            "factor in explicit controls and a separate protected condition "
            "effect in a signal row. This is not PhosR parity."
        ),
        extra_fields={
            "known_truth_source": "closed-form synthetic construction",
            "evidence_category": "synthetic_validation",
        },
        source_policy=(
            "deterministic synthetic known-truth validation fixture; not "
            "external parity and not empirical validation"
        ),
    )


def _peptide_bias_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = (
        "regime",
        "site_id",
        "sample_id",
        "true_site_abundance",
        "resolved_site_abundance",
        "signed_bias",
        "absolute_bias",
        "mapping_weight",
        "localisation_confidence",
        "bias_source",
    )
    raw_rows = [
        ("duplicate_discordant", "DUP;S1;", "sample_1", 100.0, 115.0, 1.0, 0.98),
        ("duplicate_discordant", "DUP;S1;", "sample_2", 100.0, 115.0, 1.0, 0.98),
        ("ambiguous_equal_split", "AMB;S1;", "sample_1", 100.0, 50.0, 0.5, 0.70),
        ("ambiguous_equal_split", "AMB;T2;", "sample_1", 0.0, 50.0, 0.5, 0.70),
        ("missing_observation", "MISS;S1;", "sample_1", 100.0, 80.0, 1.0, 0.95),
        ("missing_observation", "MISS;S1;", "sample_2", 100.0, 80.0, 1.0, 0.95),
        (
            "localisation_error",
            "LOC_TRUE;S1;",
            "sample_1",
            100.0,
            0.0,
            0.0,
            0.20,
        ),
        (
            "localisation_error",
            "LOC_FALSE;T2;",
            "sample_1",
            0.0,
            100.0,
            1.0,
            0.80,
        ),
    ]
    rows: list[dict[str, object]] = []
    for (
        regime,
        site_id,
        sample_id,
        truth,
        resolved,
        mapping_weight,
        localisation_confidence,
    ) in raw_rows:
        signed_bias = float(resolved) - float(truth)
        rows.append(
            {
                "regime": regime,
                "site_id": site_id,
                "sample_id": sample_id,
                "true_site_abundance": truth,
                "resolved_site_abundance": resolved,
                "signed_bias": signed_bias,
                "absolute_bias": abs(signed_bias),
                "mapping_weight": mapping_weight,
                "localisation_confidence": localisation_confidence,
                "bias_source": (
                    "duplicate peptide mean"
                    if regime == "duplicate_discordant"
                    else regime.replace("_", " ")
                ),
            }
        )
    return rows, columns


def _peptide_bias_expected_summary() -> dict[str, object]:
    return {
        "classification": "synthetic_validation",
        "estimand": (
            "site-level sample abundance after peptide-to-site signal allocation"
        ),
        "bias_units": "linear_abundance",
        "regimes": {
            "duplicate_discordant": {
                "mean_absolute_bias": 15.0,
                "supported_interpretation": (
                    "duplicate peptide evidence can bias a site mean when "
                    "duplicate peptide observations are discordant"
                ),
            },
            "ambiguous_equal_split": {
                "mean_absolute_bias": 50.0,
                "supported_interpretation": (
                    "equal splitting attenuates the true site and creates "
                    "spurious support for the alternative site"
                ),
            },
            "missing_observation": {
                "mean_absolute_bias": 20.0,
                "supported_interpretation": (
                    "missing peptide evidence shifts the resolved site "
                    "abundance toward the observed subset"
                ),
            },
            "localisation_error": {
                "mean_absolute_bias": 100.0,
                "supported_interpretation": (
                    "wrong localisation transfers signal from the true site "
                    "to the assigned false site"
                ),
            },
        },
        "tuning_policy": (
            "fixture quantifies adverse-design sensitivity and must not be used "
            "to tune production allocation parameters"
        ),
    }


def generate_peptide_site_bias_regime_fixtures(
    *,
    outdir: Path,
    seed: int,
    timestamp: str,
    command: str,
) -> None:
    family_dir = outdir / "peptide_site_bias_regimes"
    rows, columns = _peptide_bias_rows()
    bias_path = family_dir / "bias_regimes.csv"
    _write_csv(bias_path, rows, columns)
    summary_path = family_dir / "expected_bias_summary.json"
    _write_json(summary_path, _peptide_bias_expected_summary())

    _write_manifest(
        directory=family_dir,
        fixture_family="peptide_site_bias_regimes",
        classification="synthetic_validation",
        generator_command=command,
        seed=seed,
        timestamp=timestamp,
        files=(bias_path, summary_path),
        notes=(
            "Closed-form peptide-to-site adverse-regime validation fixture "
            "quantifying bias for duplicate, ambiguous, missing, and "
            "localisation-error cases. This is not external parity."
        ),
        extra_fields={
            "known_truth_source": "closed-form synthetic construction",
            "evidence_category": "synthetic_validation",
        },
        source_policy=(
            "deterministic synthetic known-truth validation fixture; not "
            "external parity and not empirical validation"
        ),
    )


def _kinase_activity_effect_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = ("site_id", "stim_effect", "null_effect")
    rows = [
        {"site_id": "UP1;S1;", "stim_effect": 3.0, "null_effect": 0.0},
        {"site_id": "UP2;S2;", "stim_effect": 2.8, "null_effect": 0.1},
        {"site_id": "UP3;S3;", "stim_effect": 2.6, "null_effect": -0.1},
        {"site_id": "DOWN1;S4;", "stim_effect": -3.0, "null_effect": 0.0},
        {"site_id": "DOWN2;S5;", "stim_effect": -2.8, "null_effect": 0.1},
        {"site_id": "DOWN3;S6;", "stim_effect": -2.6, "null_effect": -0.1},
        {"site_id": "BG1;S7;", "stim_effect": 0.0, "null_effect": 0.0},
        {"site_id": "BG2;S8;", "stim_effect": 0.2, "null_effect": 0.2},
    ]
    return rows, columns


def _kinase_activity_membership_rows() -> tuple[
    list[dict[str, object]], tuple[str, ...]
]:
    columns = ("site_id", "K_UP", "K_DOWN", "K_COVERAGE", "K_SPARSE")
    rows = [
        {
            "site_id": "UP1;S1;",
            "K_UP": 1.0,
            "K_DOWN": 0.0,
            "K_COVERAGE": 1.0,
            "K_SPARSE": 1.0,
        },
        {
            "site_id": "UP2;S2;",
            "K_UP": 1.0,
            "K_DOWN": 0.0,
            "K_COVERAGE": 1.0,
            "K_SPARSE": 0.0,
        },
        {
            "site_id": "UP3;S3;",
            "K_UP": 1.0,
            "K_DOWN": 0.0,
            "K_COVERAGE": 0.0,
            "K_SPARSE": 0.0,
        },
        {
            "site_id": "DOWN1;S4;",
            "K_UP": 0.0,
            "K_DOWN": 1.0,
            "K_COVERAGE": 0.0,
            "K_SPARSE": 0.0,
        },
        {
            "site_id": "DOWN2;S5;",
            "K_UP": 0.0,
            "K_DOWN": 1.0,
            "K_COVERAGE": 0.0,
            "K_SPARSE": 0.0,
        },
        {
            "site_id": "DOWN3;S6;",
            "K_UP": 0.0,
            "K_DOWN": 1.0,
            "K_COVERAGE": 0.0,
            "K_SPARSE": 0.0,
        },
        {
            "site_id": "BG1;S7;",
            "K_UP": 0.0,
            "K_DOWN": 0.0,
            "K_COVERAGE": 0.0,
            "K_SPARSE": 0.0,
        },
        {
            "site_id": "BG2;S8;",
            "K_UP": 0.0,
            "K_DOWN": 0.0,
            "K_COVERAGE": 0.0,
            "K_SPARSE": 0.0,
        },
    ]
    return rows, columns


def _kinase_activity_known_truth() -> dict[str, object]:
    return {
        "classification": "synthetic_validation",
        "external_reference": None,
        "effect_profile": "stim_effect",
        "membership_threshold": 0.5,
        "expected_direction": {
            "K_UP": "positive",
            "K_DOWN": "negative",
        },
        "substrate_coverage_sensitivity": {
            "kinase": "K_COVERAGE",
            "computed_min_substrates": 2,
            "insufficient_min_substrates": 3,
        },
        "sparse_membership": {
            "kinase": "K_SPARSE",
            "substrate_count": 1,
        },
        "supported_interpretation": (
            "known-membership synthetic activity direction and coverage "
            "sensitivity only; not causal kinase activation evidence"
        ),
    }


def generate_kinase_activity_known_membership_fixtures(
    *,
    outdir: Path,
    seed: int,
    timestamp: str,
    command: str,
) -> None:
    family_dir = outdir / "kinase_activity_known_membership"
    rows, columns = _kinase_activity_effect_rows()
    effects_path = family_dir / "phospho_effects.csv"
    _write_csv(effects_path, rows, columns)
    rows, columns = _kinase_activity_membership_rows()
    membership_path = family_dir / "membership_scores.csv"
    _write_csv(membership_path, rows, columns)
    truth_path = family_dir / "known_truth.json"
    _write_json(truth_path, _kinase_activity_known_truth())

    _write_manifest(
        directory=family_dir,
        fixture_family="kinase_activity_known_membership",
        classification="synthetic_validation",
        generator_command=command,
        seed=seed,
        timestamp=timestamp,
        files=(effects_path, membership_path, truth_path),
        notes=(
            "Known-membership kinase activity synthetic validation fixture "
            "covering positive/negative substrate-set direction and sensitivity "
            "to the minimum-substrate coverage rule."
        ),
        extra_fields={
            "known_truth_source": "closed-form synthetic construction",
            "evidence_category": "synthetic_validation",
        },
        source_policy=(
            "deterministic synthetic known-truth validation fixture; not "
            "external parity and not empirical validation"
        ),
    )


def _signalome_planted_rows(
    *,
    perturbed: bool,
) -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = ("site_id", "K1", "K2", "K3", "K4")
    baseline = (
        ("M1_S1", 1.00, 1.00, 0.10, 0.00),
        ("M1_S2", 1.05, 0.95, 0.15, 0.05),
        ("M1_S3", 0.95, 1.05, 0.05, -0.05),
        ("M2_S1", -1.00, -1.00, -0.10, 0.00),
        ("M2_S2", -1.05, -0.95, -0.15, -0.05),
        ("M2_S3", -0.95, -1.05, -0.05, 0.05),
    )
    perturbations = (
        (0.01, 0.00, 0.005, 0.00),
        (-0.01, 0.01, 0.000, 0.00),
        (0.00, -0.01, -0.005, 0.00),
        (0.01, 0.00, 0.005, 0.00),
        (-0.01, 0.01, 0.000, 0.00),
        (0.00, -0.01, -0.005, 0.00),
    )
    rows: list[dict[str, object]] = []
    for values, deltas in zip(baseline, perturbations, strict=True):
        site_id, *scores = values
        if perturbed:
            scores = [
                float(score) + float(delta)
                for score, delta in zip(scores, deltas, strict=True)
            ]
        rows.append(
            {
                "site_id": site_id,
                "K1": scores[0],
                "K2": scores[1],
                "K3": scores[2],
                "K4": scores[3],
            }
        )
    return rows, columns


def _signalome_known_module_rows() -> tuple[list[dict[str, object]], tuple[str, ...]]:
    columns = ("site_id", "planted_module")
    rows = [
        {"site_id": "M1_S1", "planted_module": "module_1"},
        {"site_id": "M1_S2", "planted_module": "module_1"},
        {"site_id": "M1_S3", "planted_module": "module_1"},
        {"site_id": "M2_S1", "planted_module": "module_2"},
        {"site_id": "M2_S2", "planted_module": "module_2"},
        {"site_id": "M2_S3", "planted_module": "module_2"},
    ]
    return rows, columns


def _signalome_known_truth() -> dict[str, object]:
    return {
        "classification": "synthetic_validation",
        "external_reference": None,
        "expected_module_count": 2,
        "module_recovery_metric": "pairwise_coassignment_matches_planted_modules",
        "required_pairwise_accuracy": 1.0,
        "perturbation_stability_metric": "same_pairwise_coassignment_after_perturbation",
        "required_perturbation_stability": 1.0,
        "supported_interpretation": (
            "planted module recovery and deterministic perturbation stability "
            "for a synthetic score matrix only; not causal signalome evidence"
        ),
    }


def generate_signalome_planted_module_fixtures(
    *,
    outdir: Path,
    seed: int,
    timestamp: str,
    command: str,
) -> None:
    family_dir = outdir / "signalome_planted_modules"
    rows, columns = _signalome_planted_rows(perturbed=False)
    baseline_path = family_dir / "score_matrix.csv"
    _write_csv(baseline_path, rows, columns)
    rows, columns = _signalome_planted_rows(perturbed=True)
    perturbed_path = family_dir / "score_matrix_perturbed.csv"
    _write_csv(perturbed_path, rows, columns)
    rows, columns = _signalome_known_module_rows()
    modules_path = family_dir / "planted_modules.csv"
    _write_csv(modules_path, rows, columns)
    truth_path = family_dir / "known_truth.json"
    _write_json(truth_path, _signalome_known_truth())

    _write_manifest(
        directory=family_dir,
        fixture_family="signalome_planted_modules",
        classification="synthetic_validation",
        generator_command=command,
        seed=seed,
        timestamp=timestamp,
        files=(baseline_path, perturbed_path, modules_path, truth_path),
        notes=(
            "Signalome planted-module synthetic validation fixture covering "
            "module recovery, automatic module-count selection, and small "
            "deterministic perturbation stability."
        ),
        extra_fields={
            "known_truth_source": "closed-form synthetic construction",
            "evidence_category": "synthetic_validation",
        },
        source_policy=(
            "deterministic synthetic known-truth validation fixture; not "
            "external parity and not empirical validation"
        ),
    )


def _importer_fixture_index() -> dict[str, object]:
    referenced_paths = (
        "tests/fixtures/maxquant/phospho_sty_sites_standard.txt",
        "tests/fixtures/maxquant/phospho_sty_sites_realistic_variants.txt",
        "tests/fixtures/maxquant/phospho_sty_sites_raw_and_lfq.txt",
        "tests/fixtures/maxquant/phospho_sty_sites_lfq_only.txt",
        "tests/fixtures/maxquant/phospho_sty_sites_multisite.txt",
        "tests/fixtures/fragpipe/ptmprophet_sites.tsv",
        "tests/fixtures/fragpipe/ptmprophet_explicit_site_edge_cases.tsv",
        "tests/fixtures/fragpipe/ptmprophet_peptide_position_edge_cases.tsv",
        "tests/fixtures/fragpipe/ptmprophet_missing_required_start.tsv",
    )
    return {
        "classification": "regression",
        "evidence_category": "regression",
        "external_reference": None,
        "referenced_fixture_files": [
            {
                "relative_path": relative_path,
                "sha256": _sha256(ROOT / relative_path),
            }
            for relative_path in referenced_paths
        ],
        "edge_case_coverage": {
            "maxquant": [
                "standard phospho STY columns",
                "contaminant and reverse filtering or flagging",
                "raw/LFQ intensity ambiguity",
                "LFQ-only intensity detection",
                "protein/gene group collapse diagnostics",
                "multi-site candidates",
                "localisation probability normalization",
            ],
            "fragpipe_ptmprophet": [
                "single-site peptide parsing",
                "multi-site peptide evidence",
                "ambiguous localisation diagnostics",
                "protein-position and peptide-position localisation strings",
                "contaminant and decoy filtering or flagging",
                "protein-group collapse diagnostics",
                "missing required protein start rejection",
            ],
        },
        "supported_interpretation": (
            "targeted importer edge-case regression coverage only; not broad "
            "vendor parity, not Spectronaut/DIA-NN support, and not upstream "
            "statistical-result import"
        ),
    }


def generate_importer_edge_case_manifest(
    *,
    outdir: Path,
    seed: int,
    timestamp: str,
    command: str,
) -> None:
    family_dir = outdir / "importer_edge_cases"
    index_path = family_dir / "fixture_index.json"
    _write_json(index_path, _importer_fixture_index())

    _write_manifest(
        directory=family_dir,
        fixture_family="importer_edge_cases",
        classification="regression",
        generator_command=command,
        seed=seed,
        timestamp=timestamp,
        files=(index_path,),
        notes=(
            "Manifest for targeted importer edge-case regression fixtures. "
            "Referenced files live under tests/fixtures/maxquant and "
            "tests/fixtures/fragpipe and are not external parity evidence."
        ),
        extra_fields={
            "evidence_category": "regression",
        },
    )


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    manifest_outdir_label = args.manifest_outdir_label or Path(args.outdir).as_posix()
    command = (
        "python scripts/active/generate_release_validation_regression_fixtures.py "
        f"--outdir {manifest_outdir_label} "
        f"--timestamp {args.timestamp} --seed {args.seed}"
    )
    generate_evidence_resolution_fixtures(
        outdir=outdir,
        seed=args.seed,
        timestamp=args.timestamp,
        command=command,
    )
    generate_kinase_sparse_support_fixtures(
        outdir=outdir,
        seed=args.seed,
        timestamp=args.timestamp,
        command=command,
    )
    generate_signalome_safety_fixtures(
        outdir=outdir,
        seed=args.seed,
        timestamp=args.timestamp,
        command=command,
    )
    generate_sps_ruv_planted_unwanted_factor_fixtures(
        outdir=outdir,
        seed=args.seed,
        timestamp=args.timestamp,
        command=command,
    )
    generate_peptide_site_bias_regime_fixtures(
        outdir=outdir,
        seed=args.seed,
        timestamp=args.timestamp,
        command=command,
    )
    generate_kinase_activity_known_membership_fixtures(
        outdir=outdir,
        seed=args.seed,
        timestamp=args.timestamp,
        command=command,
    )
    generate_signalome_planted_module_fixtures(
        outdir=outdir,
        seed=args.seed,
        timestamp=args.timestamp,
        command=command,
    )
    generate_importer_edge_case_manifest(
        outdir=outdir,
        seed=args.seed,
        timestamp=args.timestamp,
        command=command,
    )
    print(f"Release-validation regression fixtures written to: {outdir}")


if __name__ == "__main__":
    main()

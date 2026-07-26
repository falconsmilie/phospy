#!/usr/bin/env python3
"""Generate compact release-validation regression fixtures."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_DIR = ROOT / "tests" / "fixtures" / "release_validation_regression"
DEFAULT_TIMESTAMP = "2026-07-24T00:00:00Z"
DEFAULT_SEED = 20260724


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
    return parser.parse_args()


def _write_csv(
    path: Path, rows: list[dict[str, object]], columns: tuple[str, ...]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


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
        "source_policy": (
            "deterministic PhosPy regression contract fixture; not external parity"
        ),
        "software_versions": {
            "python": platform.python_version(),
        },
        "notes": notes,
        "files": [
            {
                "relative_path": path.relative_to(directory).as_posix(),
                "sha256": _sha256(path),
            }
            for path in files
        ],
    }
    (directory / "MANIFEST.json").write_text(
        json.dumps(payload, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )


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
    expected_path.write_text(
        json.dumps(_kinase_expected_contracts(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )

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
    snapshot_path.write_text(
        json.dumps(_historical_threshold2_snapshot(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )

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


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    command = (
        "python scripts/active/generate_release_validation_regression_fixtures.py "
        f"--outdir {Path(args.outdir).as_posix()} "
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
    print(f"Release-validation regression fixtures written to: {outdir}")


if __name__ == "__main__":
    main()

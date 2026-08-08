# Script Layout

Repository scripts are split by maintenance status:

- `active/`: supported maintainer generators used by default workflows
- `support/`: helper modules used by active generators

Default maintainer fixture bootstrap paths are wired through `Makefile`
targets that execute scripts in `active/`.

## Active Maintainer Scripts

The following scripts are protected release infrastructure for parity,
provenance, and public-workflow reference generation.

| Script | Purpose | Expected output location(s) |
| --- | --- | --- |
| `scripts/active/generate_l6_prediction_parity_fixtures.py` | Regenerate Python L6 scoring/prediction parity fixtures. | `tests/fixtures/rewrite_parity/r_reference_l6_prediction/` (`native_profile_scores.csv`, `native_rank_weighted_fusion_scores.csv`, `native_score_fusion_weights.csv`, `predMat.csv`, `native_prediction_top30.csv`, `native_candidate_substrates.csv`) |
| `scripts/active/generate_provenance_goldens.py` | Refresh provenance golden JSON fixtures used by integration contract tests. | `tests/fixtures/public_workflow_reference/kinase_public_predmat_provenance_golden.json`, `tests/fixtures/public_workflow_reference/signalome_l6_provenance_golden.json` |
| `scripts/active/generate_public_predmat_rewrite_reference.py` | Regenerate rewrite-owned public predMat lane references. | `tests/fixtures/public_workflow_reference/predmat_rewrite_stable.csv`, `tests/fixtures/public_workflow_reference/predmat_rewrite_r_parity.csv` |
| `scripts/active/generate_release_validation_regression_fixtures.py` | Regenerate compact PhosPy-owned release-validation regression and synthetic validation fixtures. | `tests/fixtures/release_validation_regression/` including `evidence_resolution/peptide_evidence.csv`, `evidence_resolution/MANIFEST.json`, `kinase_sparse_support/substrate_map.csv`, `signalome_safety/clustering_missing_dimensions.csv`, `sps_ruv_planted_unwanted_factor/`, `peptide_site_bias_regimes/`, `kinase_activity_known_membership/`, `signalome_planted_modules/`, and `importer_edge_cases/` |
| `scripts/active/generate_signalome_public_workflow_reference.py` | Regenerate rewrite-owned signalome public-workflow reference fixtures and contract metadata. | `tests/fixtures/public_workflow_reference/` (`signalome_rewrite_l6_module_assignments.csv`, `signalome_rewrite_l6_modules.csv`, `signalome_rewrite_l6_network_nodes.csv`, `signalome_rewrite_l6_network_edges.csv`, `signalome_rewrite_l6_expanded_signalome.csv`, `signalome_rewrite_l6_contract.json`) |
| `scripts/active/generate_large_differential_limma_trend_fixture.R` | Regenerate the large-feature R/limma trend differential parity fixture. | `tests/fixtures/rewrite_parity/differential_limma_trend_large/` (`matrix.csv`, `limma_B_vs_A.csv`, `MANIFEST.json`) |
| `scripts/active/generate_r_l6_fixtures.R` | Regenerate R/PhosR-side L6 parity fixtures and prediction trace artefacts. | `tests/fixtures/rewrite_parity/r_reference_l6/` and `tests/fixtures/rewrite_parity/r_reference_l6/prediction_trace/` |
| `scripts/run_pyright.py` | Resolve a suitable interpreter and run repository pyright checks. Pair with `tools/testing/pyright_strict_coverage.py --check` for strict-scope policy enforcement. | No fixture output; forwards diagnostics to stdout/stderr. |
| `scripts/verify_installed_distributions.py` | Install and execute exactly one built wheel and one built sdist outside the checkout, including installed import-origin, bundled-resource hash, and representative public workflow checks. | No fixture output; writes only temporary virtual environments and prints JSON status. |

Generators that write manifest-governed text fixtures must emit UTF-8 bytes
with LF line endings and a final newline. They must hash the exact bytes they
write, not platform-translated text-mode output. The text fixture extensions
under `tests/fixtures` are pinned to LF in `.gitattributes` so Git checkout and
archive bytes match the manifest digests without marking binary assets as text.

## Support Modules

- `scripts/support/public_workflow_reference.py`: shared helper for building
  supported L6 dataset inputs used by active public-workflow generators.

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
| `scripts/active/generate_signalome_public_workflow_reference.py` | Regenerate rewrite-owned signalome public-workflow reference fixtures and contract metadata. | `tests/fixtures/public_workflow_reference/` (`signalome_rewrite_l6_module_assignments.csv`, `signalome_rewrite_l6_modules.csv`, `signalome_rewrite_l6_network_nodes.csv`, `signalome_rewrite_l6_network_edges.csv`, `signalome_rewrite_l6_expanded_signalome.csv`, `signalome_rewrite_l6_contract.json`) |
| `scripts/active/generate_r_l6_fixtures.R` | Regenerate R/PhosR-side L6 parity fixtures and prediction trace artefacts. | `tests/fixtures/rewrite_parity/r_reference_l6/` and `tests/fixtures/rewrite_parity/r_reference_l6/prediction_trace/` |
| `scripts/run_pyright.py` | Resolve a suitable interpreter and run repository pyright checks. | No fixture output; forwards diagnostics to stdout/stderr. |

## Support Modules

- `scripts/support/public_workflow_reference.py`: shared helper for building
  supported L6 dataset inputs used by active public-workflow generators.

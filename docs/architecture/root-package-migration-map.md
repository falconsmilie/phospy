# Root Package Migration Map for ADR 0004

This document is the historical migration record for the ADR 0004 domain
refactor. For the current contributor-facing package layout, use
[`docs/architecture/package-layout.md`](./package-layout.md).

## Goals of Ticket 1

- create the target top-level package skeleton
- make package ownership explicit with package-level docstrings
- map the current flat root modules into target domain packages
- identify the files that must be split rather than moved intact

## Top-Level Package Responsibilities

| Package | Responsibility |
| --- | --- |
| `phospy.api` | Supported public entry points and thin orchestration only |
| `phospy.datasets` | Dataset models, dataset loading, builders, and dataset result containers |
| `phospy.preprocessing` | Raw-to-analysis-ready preprocessing, including filtering, correction, normalisation, imputation, and site-matrix preparation |
| `phospy.prediction` | Prediction engines, sampling, scoring components, predMat execution paths, and prediction result models |
| `phospy.activities` | Kinase activity analysis logic and activity result models |
| `phospy.signalomes` | Signalome construction and downstream map and network outputs |
| `phospy.references` | Bundled reference assets, species and reference resolution, substrate maps, motif resources, and site-sequence resources |
| `phospy.io` | Shared structured data access and path-oriented I/O helpers |
| `phospy.validation` | Validation organised by validation type |
| `phospy.errors` | Shared non-validation application error hierarchy as it is extracted from existing modules |
| `phospy.internal` | Truly internal support code only |

## Current Root Module Migration Map

| Current module or path | Target module or package | Action | Notes |
| --- | --- | --- | --- |
| `src/phospy/__init__.py` | `src/phospy/__init__.py` | Keep and slim later | Remains the package root and public export boundary while the domain move proceeds |
| `src/phospy/_dataset_validation.py` | `src/phospy/datasets/validation.py` | Move later | Dataset-loading validation helpers should live with dataset construction, while cross-cutting validation rules stay in `validation/` |
| `src/phospy/_preprocessing_primitives.py` | `src/phospy/preprocessing/primitives.py` | Completed in ticket 4 | Internal preprocessing building blocks now live inside the preprocessing domain; the legacy flat module was removed |
| `src/phospy/_protein_correction.py` | `src/phospy/preprocessing/protein_correction.py` | Completed in ticket 4 | Protein correction now lives inside the preprocessing domain; the legacy flat module was removed |
| `src/phospy/activities` | `src/phospy/activities/` | Expanded in ticket 6 | The activity domain now owns explicit analysis, results, and scoring modules under `activities/` |
| `src/phospy/analysis.py` | `src/phospy/activities/analysis.py` and `src/phospy/activities/results.py` | Completed in ticket 6 | Activity analysis now lives in `activities/analysis.py` and activity result models live in `activities/results.py`; the legacy flat module was removed earlier and the old bundled analyzer module is now gone |
| `src/phospy/cli.py` | `src/phospy/api/cli.py` | Move later | CLI is part of the supported entry surface |
| `src/phospy/constants.py` | `src/phospy/internal/constants.py` | Completed in ticket 9 | Shared constants now live under `internal/`; the legacy flat module was removed |
| `src/phospy/core_processing.py` | `src/phospy/preprocessing/core.py` | Completed in ticket 4 | Core preprocessing configuration, orchestration, and results now live under `preprocessing/`; the legacy flat module was removed |
| `src/phospy/dataset.py` | `src/phospy/datasets/models.py` | Completed in ticket 3 | Dataset models and analysis-ready dataset containers now live under `datasets/`; the legacy flat module was removed |
| `src/phospy/dataset_loader.py` | `src/phospy/datasets/loaders.py` | Completed in ticket 3 | Dataset loading now lives under `datasets/`; the legacy flat module was removed |
| `src/phospy/dataset_preprocessing.py` | `src/phospy/preprocessing/dataset.py` | Completed in ticket 4 | The dataset-bound preprocessing façade now lives under `preprocessing/`; the legacy flat module was removed |
| `src/phospy/dataset_schema.py` | `src/phospy/datasets/schema.py` | Completed in ticket 3 | Dataset schema now lives under `datasets/`; the legacy flat module was removed |
| `src/phospy/dataset_site_matrix.py` | `src/phospy/datasets/builders.py` | Completed in ticket 3 | The dataset-bound site-matrix façade now lives under `datasets/`; the legacy flat module was removed |
| `src/phospy/io` | `src/phospy/io/` | Completed in ticket 9 | Shared table readers, mapping loaders, output publishing, and output writers now live under `io/`; the legacy flat modules were removed |
| `src/phospy/matrices.py` | `src/phospy/preprocessing/site_matrix.py` | Move later | Site-matrix construction is preprocessing behaviour |
| `src/phospy/motifs.py` | `src/phospy/references/` and `src/phospy/prediction/motif_scoring.py` | Split later | This file mixes bundled reference resolution with motif scoring and must be split rather than moved whole |
| `src/phospy/pipeline.py` | `src/phospy/pipeline.py` and `src/phospy/io/` | Partially split in ticket 9 | Pipeline orchestration remains in place for now, but output publishing and writer concerns moved under `io/` |
| `src/phospy/preprocessing` | `src/phospy/preprocessing/` | Expanded in ticket 4 | The preprocessing package now owns core orchestration, services, site-matrix building, primitives, protein correction, step helpers, dataset-bound preprocessing, and explicit analysis-ready modes |
| `src/phospy/preprocessing_services.py` | `src/phospy/preprocessing/services.py` | Completed in ticket 4 | Preprocessing service objects now live inside the preprocessing domain; the legacy flat module was removed |
| `src/phospy/profiles.py` | `src/phospy/prediction/profiles.py` | Move later | Kinase substrate profile generation feeds prediction scoring |
| `src/phospy/publishing.py` | `src/phospy/io/publishing.py` | Completed in ticket 9 | Output publication now lives in `io/publishing.py`; the legacy flat module was removed |
| `src/phospy/scoring.py` | `src/phospy/prediction/scoring.py` | Moved | Kinase scoring now lives in the prediction domain |
| `src/phospy/signalome_assignments.py` | `src/phospy/signalomes/assignments.py` | Moved | Signalome-specific process module |
| `src/phospy/signalome_clustering.py` | `src/phospy/signalomes/clustering.py` | Moved | Signalome-specific process module |
| `src/phospy/signalome_construction.py` | `src/phospy/signalomes/analysis.py` | Moved | Signalome construction and execution now live inside the signalome domain |
| `src/phospy/signalome_maps.py` | `src/phospy/signalomes/maps.py` | Moved | Signalome map outputs belong inside the signalome domain |
| `src/phospy/signalome_models.py` | `src/phospy/signalomes/results.py` | Moved | Signalome result models belong inside the signalome domain |
| `src/phospy/signalome_networks.py` | `src/phospy/signalomes/networks.py` | Moved | Signalome network outputs belong inside the signalome domain |
| `src/phospy/signalome_site_ids.py` | `src/phospy/signalomes/site_ids.py` | Moved | Signalome-specific site identifier parsing and resolution |
| `src/phospy/signalomes` | `src/phospy/signalomes/` | Created now, split later | The package now exists; later tickets can collapse the remaining signalome modules into it |
| `src/phospy/site_matrix_builder.py` | `src/phospy/preprocessing/site_matrix.py` | Completed in ticket 4 | Site-matrix building now lives inside the preprocessing domain; the legacy flat module was removed |
| `src/phospy/types.py` | `src/phospy/internal/types.py` | Completed in ticket 9 | Shared internal type aliases now live under `internal/`; the legacy flat module was removed |
| `src/phospy/workflow.py` | `src/phospy/api/` workflow modules | Completed in ticket 2 | Supported public workflow entry points now live in `api/`; remaining domain helpers were pushed toward `prediction/`, `preprocessing/`, `references/`, `activities/`, and `signalomes/` |
| `src/phospy/writers.py` | `src/phospy/io/writers.py` | Completed in ticket 9 | Output writers now live under `io/writers.py`; the legacy flat module was removed |
| `src/phospy/prediction/models.py` | `src/phospy/prediction/results.py` | Completed in ticket 5 | Prediction result models now live in an explicit `results.py` module |
| `src/phospy/prediction/service.py` | `src/phospy/prediction/engines.py` | Completed in ticket 5 | Prediction execution services now live in `engines.py` with the rest of the execution layer |
| `src/phospy/prediction/workflows.py` | `src/phospy/prediction/engines.py` | Completed in ticket 5 | Workflow-bound prediction execution now lives with the rest of the prediction engine code |
| `src/phospy/prediction/` | `src/phospy/prediction/` | Keep and continue refining | This domain package already exists and should remain the home for prediction logic |
| `src/phospy/validation/` | `src/phospy/validation/` | Keep and continue refining | The validation refactor already established the right top-level direction; shared validation error classes now live under `phospy.errors` |
| `src/phospy/data/` | `src/phospy/references/assets/` | Move later | Bundled biological reference assets should sit with the references domain |

## Worst Mixed-Responsibility Files That Must Be Split

These files should not be relocated intact.

- `src/phospy/pipeline.py`
  - mixes public pipeline entry behaviour, request loading, and output publication coordination
  - should be split between `api/` and `io/`, with delegation into datasets, preprocessing, prediction, and activities
- `src/phospy/dataset_preprocessing.py`
  - still expresses preprocessing through a dataset-bound façade even though preprocessing strategy belongs in `preprocessing/`
  - later work should move this façade under `preprocessing/` while keeping `datasets/` focused on dataset ownership and construction
- `src/phospy/motifs.py`
  - mixes biological reference handling with motif-scoring behaviour used in prediction
  - should be split between `references/` and `prediction/`
- `src/phospy/dataset_site_matrix.py`
  - mixes a dataset-bound façade with site-matrix building that belongs to preprocessing
  - should be split so site-matrix process logic lives under `preprocessing/`

## Notes for Follow-Up Tickets

- Ticket 2 moved the supported public workflows into `src/phospy/api/` workflow modules without adding compatibility shims.
- The package skeleton is now in place so later tickets can move domain code into its long-term home incrementally.
- Existing `prediction/` and `validation/` packages were retained and documented rather than recreated.

- Ticket 4 centralised the preprocessing implementation under `src/phospy/preprocessing/`, removed the flat preprocessing modules, and introduced explicit full-mode and phospho-only analysis-ready preprocessing paths used by the public workflow surface.

## Final Status

The domain refactor is complete enough that this document should be treated as
historical context rather than the main contributor guide. The current package
layout and root-package policy are described in
[`docs/architecture/package-layout.md`](./package-layout.md).

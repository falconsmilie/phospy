# Root package migration map for ADR-0004

This document turns ADR-0004 into a concrete first-pass migration map for the
current `src/phospy/` tree.

## Goals of ticket 1

- create the target top-level package skeleton
- make package ownership explicit with package-level docstrings
- map the current flat root modules into target domain packages
- identify the files that must be split rather than moved intact

## Top-level package responsibilities

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

## Current root module migration map

| Current module or path | Target module or package | Action | Notes |
| --- | --- | --- | --- |
| `src/phospy/__init__.py` | `src/phospy/__init__.py` | Keep and slim later | Remains the package root and public export boundary while the domain move proceeds |
| `src/phospy/_dataset_validation.py` | `src/phospy/datasets/validation.py` | Move later | Dataset-loading validation helpers should live with dataset construction, while cross-cutting validation rules stay in `validation/` |
| `src/phospy/_preprocessing_primitives.py` | `src/phospy/preprocessing/primitives.py` | Move later | Internal preprocessing building blocks belong inside the preprocessing domain |
| `src/phospy/_protein_correction.py` | `src/phospy/preprocessing/protein_correction.py` | Move later | Protein correction is a preprocessing concern |
| `src/phospy/activities` | `src/phospy/activities/` | Created now, split later | The package now exists; later tickets should separate activity calculation modules from any convenience exports |
| `src/phospy/analysis.py` | `src/phospy/activities/analyzer.py` | Move later | `KinaseActivityAnalyzer` is downstream activity orchestration rather than a root-level concern |
| `src/phospy/cli.py` | `src/phospy/api/cli.py` | Move later | CLI is part of the supported entry surface |
| `src/phospy/constants.py` | `src/phospy/internal/constants.py` | Move later | Shared constants are internal support code unless a clearer domain home emerges during later splits |
| `src/phospy/core_processing.py` | `src/phospy/preprocessing/core.py` | Move later | Core preprocessing configuration and results belong in preprocessing |
| `src/phospy/dataset.py` | `src/phospy/datasets/models.py` | Split later | Mixed dataset models and analysis-ready dataset construction need separation |
| `src/phospy/dataset_loader.py` | `src/phospy/datasets/loaders.py` | Move later | Dataset file and frame loading belongs in datasets |
| `src/phospy/dataset_preprocessing.py` | `src/phospy/preprocessing/datasets.py` | Move later | This is preprocessing that is currently expressed through a dataset-shaped façade |
| `src/phospy/dataset_schema.py` | `src/phospy/datasets/schema.py` | Move later | Dataset schema belongs with dataset construction |
| `src/phospy/dataset_site_matrix.py` | `src/phospy/preprocessing/site_matrix_bindings.py` | Split later | The bound façade is dataset-shaped, but the site-matrix process belongs in preprocessing |
| `src/phospy/io` | `src/phospy/io/` | Created now, split later | The package now exists; later tickets can break reading helpers and mapping loaders into clearer submodules |
| `src/phospy/matrices.py` | `src/phospy/preprocessing/site_matrix.py` | Move later | Site-matrix construction is preprocessing behaviour |
| `src/phospy/motifs.py` | `src/phospy/references/` and `src/phospy/prediction/motif_scoring.py` | Split later | This file mixes bundled reference resolution with motif scoring and must be split rather than moved whole |
| `src/phospy/pipeline.py` | `src/phospy/api/pipeline.py` and `src/phospy/io/publishing.py` | Split later | Contains public pipeline orchestration plus output publishing coordination |
| `src/phospy/preprocessing` | `src/phospy/preprocessing/` | Created now, split later | The package now exists; later tickets should separate façade code from stable preprocessing process modules |
| `src/phospy/preprocessing_services.py` | `src/phospy/preprocessing/services.py` | Move later | Preprocessing service objects belong inside the preprocessing domain |
| `src/phospy/profiles.py` | `src/phospy/prediction/profiles.py` | Move later | Kinase substrate profile generation feeds prediction scoring |
| `src/phospy/publishing.py` | `src/phospy/io/publishing.py` | Move later | Output publication is a structured I/O concern |
| `src/phospy/scoring.py` | `src/phospy/prediction/scoring.py` | Move later | Kinase scoring is part of prediction execution |
| `src/phospy/signalome_assignments.py` | `src/phospy/signalomes/assignments.py` | Move later | Signalome-specific process module |
| `src/phospy/signalome_clustering.py` | `src/phospy/signalomes/clustering.py` | Move later | Signalome-specific process module |
| `src/phospy/signalome_construction.py` | `src/phospy/signalomes/construction.py` | Move later | Signalome-specific process module |
| `src/phospy/signalome_maps.py` | `src/phospy/signalomes/maps.py` | Move later | Signalome map outputs belong inside the signalome domain |
| `src/phospy/signalome_models.py` | `src/phospy/signalomes/models.py` | Move later | Signalome result models belong inside the signalome domain |
| `src/phospy/signalome_networks.py` | `src/phospy/signalomes/networks.py` | Move later | Signalome network outputs belong inside the signalome domain |
| `src/phospy/signalome_site_ids.py` | `src/phospy/signalomes/site_ids.py` | Move later | Signalome-specific site identifier parsing and resolution |
| `src/phospy/signalomes` | `src/phospy/signalomes/` | Created now, split later | The package now exists; later tickets can collapse the remaining signalome modules into it |
| `src/phospy/site_matrix_builder.py` | `src/phospy/preprocessing/site_matrix_builder.py` | Move later | Site-matrix building belongs in preprocessing |
| `src/phospy/types.py` | `src/phospy/internal/types.py` | Move later | Shared internal type aliases are infrastructure support rather than a user-facing domain |
| `src/phospy/workflow.py` | `src/phospy/api/workflows.py` | Split later | This is one of the worst mixed-responsibility files because it blends public entry points with domain assembly |
| `src/phospy/writers.py` | `src/phospy/io/writers.py` | Move later | Output writers belong in structured I/O |
| `src/phospy/prediction/` | `src/phospy/prediction/` | Keep and continue refining | This domain package already exists and should remain the home for prediction logic |
| `src/phospy/validation/` | `src/phospy/validation/` | Keep and continue refining | The validation refactor already established the right top-level direction |
| `src/phospy/data/` | `src/phospy/references/assets/` | Move later | Bundled biological reference assets should sit with the references domain |

## Worst mixed-responsibility files that must be split

These files should not be relocated intact.

- `src/phospy/workflow.py`
  - mixes supported public workflows with substantial domain orchestration
  - should end as thin entry points in `api/` delegating into domain packages
- `src/phospy/pipeline.py`
  - mixes public pipeline entry behaviour, request loading, and output publication coordination
  - should be split between `api/` and `io/`, with delegation into datasets, preprocessing, prediction, and activities
- `src/phospy/dataset.py`
  - mixes dataset-shaped models with analysis-ready dataset construction behaviour
  - should be split between dataset models and builder-oriented logic in `datasets/`
- `src/phospy/motifs.py`
  - mixes biological reference handling with motif-scoring behaviour used in prediction
  - should be split between `references/` and `prediction/`
- `src/phospy/dataset_site_matrix.py`
  - mixes a dataset-bound façade with site-matrix building that belongs to preprocessing
  - should be split so site-matrix process logic lives under `preprocessing/`

## Notes for follow-up tickets

- No compatibility shims were added as part of ticket 1.
- The package skeleton is now in place so later tickets can move domain code into its long-term home incrementally.
- Existing `prediction/` and `validation/` packages were retained and documented rather than recreated.

# ADR: Competitive Phosphoproteomics Workflow Coverage Roadmap

## Document Control

- **ADR ID:** ADR-0025
- **Title:** Competitive Phosphoproteomics Workflow Coverage Roadmap
- **Status:** Accepted
- **Date:** 2026-06-11
- **Decision Type:** Architecture and Scientific Roadmap

Complements ADR-0013, ADR-0015, ADR-0019, ADR-0022, and ADR-0024.

## Status

Accepted.

## Context

PhosPy is compared against established phosphoproteomics tools and reference
surfaces such as PhosR, MSstatsPTM, and Kinase Library. Those comparisons are
useful for planning, but they can create implementation drift when the roadmap
lives only in reviews, tickets, or external commentary.

The project needs one internal source of truth that records intended direction
without turning future capabilities into present-tense feature claims.

## Decision

PhosPy will maintain competitive workflow coverage as a phased roadmap governed
by this ADR and surfaced through `docs/scientific-coverage.md`.

This ADR is not an implementation claim. A roadmap item becomes supported only
when code, public contracts, documentation, and tests exist, and when the
scientific coverage matrix is updated to the correct claim category.

## Current State

The supported scientific workflow interface is the Python API.

Current executable analysis/workflow lanes are:

- `AnalysisReadyDatasetBuilder`
- `DifferentialAnalysisWorkflow`
- `KinaseWorkflow`
- optional kinase activity score tables within the kinase workflow
- optional `SignalomeWorkflow` after kinase prediction
- `EnrichmentWorkflow` for offline over-representation analysis (ORA)

Current public input-preparation support includes generic column-mapped
phosphosite import, MaxQuant phosphosite import, and
FragPipe/Philosopher/PTMProphet phosphosite import. These importers emit
`PhosphositeImportResult` candidates and dataset-builder requests; they do not
construct `AnalysisReadyPhosphoDataset` objects, infer sample groups, infer
contrasts, infer batches, infer paired blocks, or bypass builder validation.
They are targeted importer adapters, not broad support for all vendor,
search-engine, or upstream statistical outputs. No dedicated Spectronaut or
DIA-NN semantic importer is current support; those remain future/demand-driven
candidates unless executable adapters, public contracts, docs, and tests are
added.

Current differential support has a parity-protected core lane: two-condition
unpaired simple contrasts with explicit `ExperimentalDesign` and `Contrast`
objects, empirical-Bayes `standard` or `robust` modes, optional trend
moderation, and Benjamini-Hochberg adjustment.

Additional validated PhosPy differential support includes explicit
fixed-effect covariates declared on `ExperimentalDesign`: batch, categorical
covariates, and continuous covariates. It also includes explicit complete
fixed-block paired/block designs through
`paired_design_policy="fixed_block"` when every block has complete condition
coverage for each requested contrast. These ordinary fixed-effect linear-model
terms are validated for required metadata, rank, and contrast estimability
before execution.

Fixed-effect batch terms in differential analysis are model covariates; they
do not correct the input data. They are not ComBat, not RUV, not limma
`removeBatchEffect`, not limma `duplicateCorrelation`, not mixed-effects
modelling, not random effects, and not a general batch-correction method.
Fixed-block terms are ordinary fixed effects over explicit
`SampleDesignRecord.block_id` values; they are not random subject modelling,
mixed-effects modelling, or
`duplicateCorrelation`-style correlated-replicate modelling. Correlated
repeated-measure and mixed-effect differential modelling remain unsupported.

Dataset preprocessing supports one opt-in batch-residualisation method:
`linear_residualize_batch`. It uses fixed-effect residualisation of batch terms
while preserving condition effects by including condition terms in the
residualisation design. It requires explicit `sample_metadata` and
`DatasetBatchCorrectionConfig` settings, rejects confounded or otherwise
inadequate designs before correction, and records typed diagnostics. It is not
ComBat, not RUV, not limma `removeBatchEffect` parity, not
`duplicateCorrelation`, not mixed-effects modelling, and not a solution to all
batch-effect problems.

Dataset preprocessing also supports native SPS/RUV-style correction through
`SpsRuvBatchCorrectionConfig`. This implemented lane requires caller-supplied
control `site_key` annotations, aligned batch and protected condition metadata,
optional replicate metadata for provenance, explicit `CorrectionMissingnessPolicy`,
`n_unwanted_factors`, diagnostics, and provenance. It estimates unwanted factors
from eligible control-site residuals after protecting condition terms and
applies correction before downstream workflows consume the analysis-ready
matrix. It is not PhosR-equivalent SPS/RUV-III parity, not executable RUV-III
support, not ComBat, not limma `removeBatchEffect`, and not hidden online
control fetching or control selection.

Dataset preprocessing also supports total-protein subtraction through
`subtract_log_total` and protein-aware model-input preparation through
`DatasetProteinAwarePreparationConfig(policy="prepare_model_inputs")`.
Protein-aware preparation builds matched phosphosite/protein pairs,
sample-aligned protein covariates, eligibility rows, mapping diagnostics,
sample-alignment diagnostics, transformation-state diagnostics, and explicit
limitations. It does not modify phosphosite values, subtract total protein,
normalise intensities, run differential modelling, or claim MSstatsPTM-style
joint PTM/protein inference. It does not claim MSstatsPTM-style equivalence.
Current `DifferentialAnalysisWorkflow` execution does not consume the prepared
protein covariate matrix.

Current kinase support provides profile/motif scoring, rank-weighted fusion,
deterministic/adaptive prediction, and three explicit activity methods:
`simplified_weighted_substrate_activity_v1`,
`ksea_zscore_activity_v1`, and
`ssgsea_substrate_enrichment_activity_v1`. These outputs are support summaries
or substrate-set enrichment statistics, not calibrated causal inference.

Current Kinase Library-style support includes a pure motif scorer, local
`KinaseLibraryResource` / `KinaseLibraryResourceLoader` contracts, and opt-in
`KinaseWorkflow` scoring modes `kinase_library_motif` and
`combined_profile_motif`. These paths require caller-supplied compatible local
Kinase Library-style resources with explicit provenance. PhosPy does not bundle
official Kinase Library data and does not claim validated Kinase Library parity.
Workflow motif scores are normalized support scores for within-run ranking;
they are not probabilities and do not imply kinase activity unless an explicit
activity method is run.

The ssGSEA-style kinase activity method is a PhosPy rank-walk
substrate-set activity implementation over kinase-substrate membership and
phosphosite effect/statistic values. It is not PTM-SEA parity and is separate
from the native ORA enrichment workflow. When permutations are disabled,
ssGSEA-style statistics record an explicit significance-unavailable status and
leave p-value/q-value fields missing rather than using numeric placeholders.
Seeded permutations provide permutation-significance estimates only; they do
not convert activity-like scores into causal kinase activation evidence.

Current bundled runtime reference data is rat-only. Human and mouse analyses
can be run when the caller supplies an explicit `ReferenceBundle`.

Current signalome support provides score-derived module assignments, module
summaries, clustering diagnostics, and kinase score-profile association network
tables. Signalome network correlations are descriptive associations only. They
are not inferential tests, causal evidence, directionality evidence, activity
claims, or experimentally validated signalling relationships.

Current native enrichment support is offline ORA through `EnrichmentWorkflow`
against caller-supplied `GeneSetCollection`, `PtmSetCollection`, or homogeneous
`EnrichmentSetCollection` inputs with explicit identifier semantics and an
explicit background universe. PhosPy does not bundle GO, KEGG, Reactome,
PTM-SEA, or PTMsigDB resources for this feature, and core workflow execution
does not call Enrichr, gseapy, clusterProfiler, or other online services. ORA
does not imply GSEA, ssGSEA, or PTM-SEA support.

Current core PhosPy does not provide a first-class visualization workflow/API.

Current PhosPy does not support command-line scientific workflow execution.
ADR-0022 remains authoritative: Python API is the supported scientific
workflow interface unless a future ADR changes that boundary.

PhosPy does not claim full PhosR, MSstatsPTM, or Kinase Library equivalence.

## Desired Direction

The desired direction is a competitive but scoped Python phosphoproteomics
workflow stack:

- richer reference handling with explicit provenance, compatibility checks,
  and user-supplied external reference bundles
- additional kinase inference and activity methods with method-specific
  validation and clear output meaning
- additional semantic importers beyond the current generic, MaxQuant, and
  FragPipe/PTMProphet adapters, without bypassing dataset validation or site
  identity contracts
- richer differential designs beyond ordinary fixed covariates and complete
  fixed-block terms, especially correlated repeated-measure,
  `duplicateCorrelation`-style, random-effect, mixed-effect, or additional
  batch-effect methods, only after explicit design/result contracts and parity
  or validation evidence exist
- enrichment workflows and resource integrations beyond current offline ORA,
  kept separate from kinase scoring unless the method is explicitly a kinase
  activity or substrate-set activity method
- visualization adapters that consume validated result objects without becoming
  the source of scientific truth
- possible CLI workflow support only as a thin, validated wrapper over the
  Python API after ADR-0022's reintroduction criteria are satisfied

## Non-Goals

This roadmap does not make PhosPy a clone or full replacement for PhosR,
MSstatsPTM, Kinase Library, MaxQuant, FragPipe, Spectronaut, DIA-NN, or any
other upstream processing or analysis package.

The roadmap does not permit:

- describing planned capabilities as already supported
- broad global parity claims
- hidden sample-name inference for scientific design
- silent reference remapping or protein identity guessing
- bundling reference data without redistribution permission
- adding workflow logic to visualization or CLI layers that diverges from the
  Python API
- treating statistical association, enrichment, or scoring output as causal
  biological proof

## Implementation Phases

Phases are directional, not calendar commitments.

### Phase 0: Guardrails

- Keep `docs/scientific-coverage.md` as the user-facing current support matrix.
- Keep roadmap items separate from current executable support.
- Require tickets for roadmap work to name the target claim category and
  affected workflow contracts.

### Phase 1: Reference and Import Foundations

- Harden `ReferenceBundle` provenance, compatibility diagnostics, and
  external bundle validation.
- Maintain current generic, MaxQuant, and FragPipe/PTMProphet importer
  contracts as dataset-builder input adapters, not hidden analysis engines.
- Add additional semantic importer contracts only when they emit typed tables or
  requests that still pass existing builder/workflow validation.
- Treat Spectronaut and DIA-NN importers as future/demand-driven candidates,
  not current support, until they meet the same contract, documentation, and
  test requirements as existing adapters.
- Keep raw/vendor/search-engine format interpretation out of core workflows
  unless a dedicated importer contract owns it.

### Phase 2: Differential Design Depth

- Preserve current fixed-effect covariate and complete fixed-block support as
  ordinary fixed-term linear-model execution.
- Extend experimental-design contracts before adding correlated repeated-measure
  execution, `duplicateCorrelation`-style modelling, random-effect or
  mixed-effect modelling, or additional batch-effect methods.
- Preserve explicit contrast definitions and provenance.
- Add validation and parity/evidence tests before public support claims.

### Phase 3: Kinase, Activity, and Enrichment Depth

- Add additional kinase inference/activity methods one method at a time with
  stable policy IDs, output-scale definitions, and tests.
- Maintain current offline ORA as a separate enrichment workflow, and keep any
  future broad pathway/gene-set enrichment separate from kinase scoring unless
  the method is explicitly a kinase activity method.
- Keep prediction scores, activity scores, and enrichment statistics labeled by
  their actual statistical meaning.

### Phase 4: Visualization and CLI

- Add visualization as result-consuming adapters, not as hidden analysis
  engines.
- Reintroduce CLI workflow execution only after a separate decision confirms
  Python API validation parity, complete configuration coverage, provenance,
  documentation, and tests.

## Reference-Data Redistribution Rule

PhosPy may bundle reference data only when redistribution is explicitly allowed
by the source license or written permission, and when source provenance is
documented.

If a useful reference source has absent, ambiguous, restrictive, or
non-redistributable terms, PhosPy must not bundle that data. Acceptable
alternatives are:

- document the required schema
- provide validators or adapters for user-supplied local files
- provide scripts that transform user-provided inputs locally
- include only synthetic, minimal, or otherwise redistributable fixtures

Derived reference tables inherit the redistribution limits of their sources
unless the source terms explicitly allow redistribution of derived data.

The rat-only bundled runtime reference status remains the current support
boundary until a future implementation changes the bundled data and passes
license, provenance, documentation, and test review.

## Scientific Claim Categories

Public scientific claims must use the categories maintained in
`docs/scientific-coverage.md`:

- `parity-gated`
- `validated PhosPy implementation`
- `experimental`
- `open gap`
- `deliberate scope difference`
- `not planned`

A roadmap item is not a claim category. A future direction remains an
`open gap` or `deliberate scope difference` until executable support exists. An
`experimental` claim requires executable behavior and explicit caveats. A
`parity-gated` claim requires maintained parity evidence.

## Architecture Responsibility Boundaries

Documentation boundaries:

- `docs/scientific-coverage.md` owns current user-facing support status.
- `docs/workflow_contracts.md` owns executable workflow contracts and known
  limitations.
- ADRs own rationale, direction, and governance.
- README and API guide pages may link to the roadmap but must not expand
  current support claims.

Code ownership boundaries:

- Dataset construction and semantic importers must feed the
  `AnalysisReadyDatasetBuilder` contract rather than bypassing validation.
- Differential extensions belong under differential design, workflow, result,
  and provenance modules; sample-name inference remains out of scope.
- Reference resolution belongs under reference models/resources/resolution and
  must not be embedded ad hoc in workflow code.
- Kinase scoring, prediction, activity, and enrichment methods must own stable
  scientific policy records in their domain modules.
- Visualization layers may read result models but must not mutate scientific
  outputs or implement alternate workflow semantics.
- CLI layers, if reintroduced, must delegate to Python API request/workflow
  objects and preserve the same validation, provenance, and failure behavior.

## Consequences

Future work can be planned against competitive workflow areas without
overclaiming current support.

Review, documentation, and release checks must reject roadmap language that
describes planned capabilities as available behavior.

Tickets that add scientific scope should update the coverage matrix, workflow
contracts, and tests in the same change that adds executable support.

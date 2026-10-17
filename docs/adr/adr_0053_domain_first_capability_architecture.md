# ADR-0053: Domain-First Capability Architecture and Dependency Governance

## Status

- **ADR ID:** ADR-0053
- **Title:** Domain-First Capability Architecture and Dependency Governance
- **Status:** Accepted
- **Date:** 2026-09-26
- **Accepted:** 2026-09-26
- **Amended:** 2026-10-17 by the accepted Wave 2D, Wave 3, Wave 4B, Wave 4C, Wave 5A D05A-1, Wave 5B D05B-1 and D05B-2 decisions, Wave 7B3 ownership decision, Decision D5, and Decision D08-1
- **Decision Type:** Architecture and Package Governance
- **Evidence Base:** Architecture inventory and boundary verification against source snapshot `138b7764bcb75bc3506dc58cdd2d8b044cf70d4c` (PhosPy 1.7.6)
- **Supersedes:** ADR-0010 horizontal package-ownership and top-level dependency-graph decisions
- **Retains:** ADR-0010 module-splitting, single-owner, and facade-without-behaviour rules; compatibility governance is superseded separately by ADR-0054
- **Refines:** ADR-0002, ADR-0003, ADR-0007, ADR-0019, ADR-0020, ADR-0024, ADR-0033, ADR-0034, ADR-0035, ADR-0051, and ADR-0052
- **Public API companion:** ADR-0054

## Context

PhosPy currently has a deliberate and acyclic top-level package graph. Public
workflow shells are thin, scientific domains are separated from orchestration,
validation is explicit, and workflow execution consistently follows the
validator → interpreter → executor lifecycle.

The current physical structure is nevertheless organised primarily by
horizontal responsibility:

- `phospy.contracts` aggregates requests, configuration facades, and results;
- `phospy.validation` contains shared and workflow-specific validation;
- `phospy.workflows` owns orchestration;
- `phospy.science` contains most scientific models and execution;
- `phospy.io` owns readers, bundles, and publishers; and
- `phospy.provenance` supplies shared provenance infrastructure.

This layout protects broad layer direction, but a single scientific capability
is distributed across several top-level packages. It is therefore difficult to
understand Differential, Enrichment, Kinase, Signalome, dataset construction,
or preprocessing primarily from one bounded package.

The top-level dependency graph also treats all of `phospy.science` as one
node. A current-state inventory found a large strongly connected component
among logical science domains even though the physical top-level graph remains
acyclic. Representative causes include:

- dataset preprocessing importing quantification-depth definitions from
  Differential;
- phosphosite identity importing `Organism` through References while
  References imports phosphosite identity;
- general References importing Kinase prediction constants;
- Signalome importing Kinase caveat constants and prediction internals; and
- central table/configuration packages owning domain-specific declarations.

The architecture needs to express the scientific capabilities directly,
retain legitimate shared concepts, and make cross-capability dependencies
explicit and enforceable.

This is an ownership and dependency decision. It is not a scientific rewrite.

## Decision drivers

1. A capability should be understandable mostly from its own package.
2. Configuration, validation, planning, execution, results, diagnostics, and
   domain provenance should have one clear owner.
3. Shared scientific concepts must remain shared where their meaning is
   independent of one workflow.
4. Workflow packages should follow common lifecycle conventions without being
   forced into identical directory structures.
5. The validator → interpreter → executor lifecycle and `run()` naming remain
   standard.
6. `AnalysisReadyPhosphoDataset` remains the strict downstream dataset
   boundary.
7. `site_key` remains the unique protein-scoped analysis identity, while
   `display_id` remains repeatable presentation metadata.
8. Scientific behaviour, validation outcomes, provenance, serialization, and
   dataframe-ownership guarantees must remain unchanged during structural
   migration.
9. Public API policy and physical migration sequencing must remain separate
   decisions.
10. Dependency direction must be checked at logical-domain and contract-path
    levels, not only at broad top-level packages.

## Decision

### Domain-first ownership

PhosPy will use domain-first logical ownership. The major capability domains
are:

1. **Datasets** — dataset construction and the strict analysis-ready dataset
   boundary.
2. **Preprocessing** — data preparation, state transitions, reports, nested
   batch correction, and control-site discovery.
3. **Differential** — experimental design, contrasts, eligibility, fitting,
   moderation, diagnostics, results, and workflow provenance.
4. **Enrichment** — identifier/universe semantics, set filtering, ORA,
   diagnostics, results, and provenance.
5. **Kinase** — reference projection, scoring, prediction, kinase-specific
   activity, eligibility, results, and provenance.
6. **Signalome** — Kinase-result interpretation, score conditioning,
   clustering, module selection, network/context construction, diagnostics,
   results, and provenance.

The shared scientific domains are:

- `activities` for generic activity quantitative semantics, passive activity
  profile and aggregation contracts, generic method quantitative-input
  contracts, contract resolution and rendering, and generic validation of
  those passive contracts;
- `phosphosites` for site identity and sequence-context rules;
- `taxonomy` for the narrow `Organism` vocabulary and normalisation;
- `quantitative` for scale, quantitative meaning, transitions, and generic
  transformation evidence;
- `evidence` for peptide/localisation evidence and reusable site-level
  evidence semantics;
- `references` for general reference identity, manifests, resource
  resolution, and protein-sequence resources;
- `statistics` for genuinely reusable statistical kernels; and
- `provenance` for generic immutable audit, JSON, environment, fingerprint,
  table and scientific-policy primitives; and
- `run_lineage` for cross-capability `RunProvenance` composition and its
  direct codec.

Technical foundations such as errors, dataframe ownership/schema mechanics,
generic scalar parsing, and generic validation mechanics remain shared leaves.
Public facades and IO remain horizontal adapters rather than scientific owners.

### Capability lifecycle ownership

A capability owns the lifecycle artefacts that have capability-specific
meaning:

- requests and configuration;
- boundary and contextual validation;
- interpretation and immutable execution plans;
- executors and scientific kernels;
- result contracts and result construction;
- diagnostics and caveat selection;
- domain-specific provenance assembly; and
- passive serialization schemas or codecs.

This convention does not require identical directories. A compact capability
such as Enrichment may use focused modules, while Kinase and preprocessing may
use nested subpackages. Responsibility and dependency direction define
consistency, not directory count.

Generic validation primitives may remain in a shared validation foundation.
Workflow-specific validation policy and composition belong to the owning
capability. Validation modules must not own orchestration or execution.

This refines ADR-0007's former physical-location rule. “If it is validation,
it belongs in validation” now means that every validator belongs in a clearly
named `validation` module or subpackage. Capability-specific validation lives
under its capability, for example `differential.validation` or
`preprocessing.validation`; only reusable validation mechanics and primitive
invariant checks remain under root `phospy.validation`. This is an ownership
change, not permission to place validation logic in orchestration, contracts,
or numerical kernels.

### Dataset and preprocessing direction

Datasets and preprocessing are separate first-class capabilities with a
one-way dependency:

`datasets → preprocessing`

Dataset construction invokes preprocessing and consumes its passive outputs.
Preprocessing must not import the analysis-ready dataset implementation.
Preprocessing produces processed tables, state, evidence, reports, and
provenance inputs; dataset construction applies final construction invariants
and creates `AnalysisReadyPhosphoDataset`.

Within Datasets, the analysis-ready dataset model is a downstream-consumed
leaf contract and must not depend on preprocessing. Only dataset construction,
planning, and builder composition may depend on preprocessing. Architecture
tests must distinguish these internal roles rather than treating the broad
`datasets` package edge as permission for every dataset module to import
preprocessing.

Downstream analytical workflows depend on the dataset model and published
passive contracts, not builder internals or preprocessing-stage
implementations.

Dataset-build acceptance and workflow-specific scientific eligibility remain
separate decisions.

### Nested preprocessing capabilities

Batch correction remains a nested preprocessing capability even though it has
its own validator → interpreter → executor lifecycle.

SPS discovery is an independently runnable subcapability of
preprocessing/control-site selection. Discovery and correction are separate
scientific operations:

- discovery consumes explicit governed reference datasets;
- discovery returns `SpsDiscoveryResult` and a reusable `ControlSiteSet`;
- correction consumes `ControlSiteSet` without consuming the discovery
  implementation or source reference matrices; and
- `sps_ruv_style` and `ruv_iii_style` remain preprocessing/batch-correction
  methods.

SPS discovery is not a peer top-level analytical capability because the
current implementation has no consumer outside preprocessing and batch
correction.

### Differential ownership

Experimental design, `SampleDesignRecord`, contrasts, design-matrix
construction, and replicate policy are Differential-owned. The current source
has no independent production consumer outside Differential requests,
validation, planning, execution, provenance, and facade exposure.

The public exposure of these types does not make Experimental Design a shared
implementation domain. Facades may re-export the owner-defined object
identities.

Quantification-depth ownership is divided by meaning:

- Evidence or a dataset-evidence subdomain owns the
  `quantification_depth` column contract, supported count kinds, integer/count
  validation, reader/importer handling, and conservative site-level
  aggregation.
- Differential owns the log2 depth-trend transformation, empirical-Bayes
  configuration and interpretation, trend fitting, diagnostics, and
  provenance.

Readers, Evidence, preprocessing, and datasets must not import Differential to
obtain the upstream count-evidence contract.

### Kinase ownership

Kinase is an umbrella capability containing the following subdomains:

- reference projection and Kinase Library interpretation;
- scoring;
- prediction; and
- kinase-specific activity behaviour.

General References may load and identify generic reference resources. It must
not own or import Kinase motif/profile constants, scoring semantics, or
projection behaviour.

Decision D5 assigns the stored kinase–substrate reference-row contract to
References. This includes structural and content validation, identifier
normalisation, duplicate and conflicting-pair classification, validation
reports, and dataframe ownership. Kinase continues to own reference
projection, Kinase Library interpretation, motif/profile meaning, scoring,
prediction, kinase-specific activity, and scientific execution. A stored table
contract does not permit References to import those Kinase implementations.

The accepted Wave 7B3 ownership decision places the passive shared activity
boundary under `phospy.activities`, while kinase-specific configuration,
methods, eligibility, membership, results and activity-stage coordination
remain under `phospy.kinase.activity`. Kinase may consume the passive Activities
boundary, and IO may consume its passive published contracts. These are inbound
uses and do not broaden the Activities outbound dependency allow-list.

### Signalome dependency on Kinase

Signalome consumes the complete `KinaseWorkflowResult` directly. No
`KinaseSignalomeInput`, score-view DTO, or equivalent intermediary result is
introduced.

Signalome may use the dataset, scoring result, prediction result, and other
nested published result contracts carried by `KinaseWorkflowResult`. It must
not import:

- Kinase validators, interpreters, or executors;
- Kinase scoring or prediction implementation modules;
- Kinase internal-view classes;
- Kinase caveat constants; or
- private Kinase assemblers or provenance helpers.

The Kinase result and its nested result contracts must expose the information
that Signalome legitimately requires. The internal copy/borrow mechanism used
to preserve dataframe-ownership and performance guarantees is not decided by
this ADR.

Signalome also has one permanent narrow dependency on the Dataset boundary carried by `KinaseWorkflowResult`. `phospy.signalome.results` may import `phospy.datasets.models.AnalysisReadyPhosphoDataset` and `phospy.datasets.processing_state.DatasetPreprocessingReport`; `phospy.signalome.result_validation` and `phospy.signalome.contracts` may each import `phospy.datasets.models.AnalysisReadyPhosphoDataset` and `phospy.datasets.internal_view.DatasetInternalView`; and `phospy.signalome.validator` may import `phospy.datasets.internal_view.DatasetInternalView`. No other Signalome consumer or Dataset symbol is permitted. `DatasetInternalView` and its frame store remain Dataset-owned. This permission does not extend to Dataset `_internal`, construction, builder or validation modules, private Dataset storage, copied view implementations, public borrowing APIs, replacement DTOs or graph exclusions.

### Phosphosite identity and taxonomy

The shared site-identity package is named `phosphosites`, not the generic
`identity`.

`phosphosites` owns:

- `site_key` construction and parsing;
- display-site identifiers and parsing;
- protein/isoform scope;
- identity columns and identity validation rules; and
- representation and validation of an already supplied site-centred sequence
  window.

General References owns protein-sequence repositories, FASTA parsing, and
sequence retrieval. Evidence or preprocessing owns resolving a phosphosite
against a retrieved protein sequence, checking residue/window agreement, and
recording the resulting evidence and provenance. `phosphosites` must not
import References to perform resource lookup or sequence resolution.

The shared organism package is named `taxonomy`. It is a narrow leaf that owns
only typed organism values and normalisation. It must not become a general
biological-constants package.

References may depend on `phosphosites` and `taxonomy`. `phosphosites` must not
depend on the References implementation merely to obtain `Organism`.

### Tables, configuration, contracts, and results

Generic dataframe ownership and schema mechanics remain in `frames`.
Domain-specific table schemas belong to their owning capability or shared
domain.

Domain-specific configuration, requests, results, validation, and diagnostics
belong to their owning capability. Aggregate facade modules may re-export exact
owner-defined identities but must not become second implementation owners.

Internal code imports from the owning domain rather than through `phospy.api`,
`phospy.advanced`, or an aggregate compatibility facade.

### Provenance and serialization

Shared provenance owns only generic immutable audit, JSON, environment,
fingerprint, table and scientific-policy primitives. Each capability owns the
meaning and assembly of its workflow-specific provenance.

The accepted Wave 4C decision assigns `RunProvenance`, its run-level
scientific-policy validator and its codec to the higher-level
`phospy.run_lineage` composition owner. Run lineage depends one-way on generic
Provenance primitives, References records/codecs, the Evidence historical
peptide codec and Preprocessing stage provenance/codec. References, Evidence
and Preprocessing must not depend back on run lineage. Dataset, Enrichment,
Kinase, Signalome and IO may consume its published composition boundary.

Within protein-aware preparation, Preprocessing owns the passive
`ProteinAwarePreparationStageOutput` produced by its stage. Datasets owns the
durable `ProteinAwarePreparationReport` and `ProteinAwarePreparationResult`
in `phospy.datasets.protein_aware` and composes them one-way from the
stage output. The passive stage output contains no run lineage or Dataset
schema state, and Preprocessing does not import Datasets or `run_lineage`.

The accepted Wave 4B decision makes one narrow audit-contract distinction.
Provenance owns the immutable batch-correction wire records and their direct
historical codec. Preprocessing owns correction policy, scientific validation,
environment selection, active provenance assembly and execution. This does not
generalise to other capability-specific provenance records and introduces no
Provenance-to-Preprocessing dependency.

Capabilities may own passive, versioned serialization schemas or codecs for
their results. Concrete filesystem paths, transactions, integrity checks,
bundle storage, readers, and publishers remain IO adapter responsibilities.
Capabilities must not import concrete IO adapters.

For Enrichment, generic GMT, CSV, and TSV acquisition remains IO-owned. Those
readers may import only Enrichment's passive collection contracts through a
direct owner leaf; they may not import Enrichment validation, planning,
execution, ORA, result assembly, workflow, or provenance modules. This
clarifies the existing IO-to-published-contract rule and does not add a new
logical owner edge.

The Wave 2D preflight decision narrows that reader boundary for one
owner-specific acquisition service. References owns local reference
source-table acquisition when it is the default dependency of the
owner-defined `ReferenceBundleBuilder`. General CSV, TSV, and Parquet readers,
bundle storage, publishers, filesystem integrity checks, and unrelated IO
adapters remain IO-owned. Shared table-format recognition is a technical
foundation; References must not import IO to obtain it.

### Target logical package roots

The target logical ownership roots are:

```text
phospy/
├── datasets/
├── preprocessing/
├── preparation_records/
├── differential/
├── enrichment/
├── kinase/
├── signalome/
├── activities/
├── phosphosites/
├── taxonomy/
├── quantitative/
├── evidence/
├── references/
├── statistics/
├── provenance/
├── run_lineage/
├── frames/
├── validation/       # generic mechanics only
├── policies/         # generic parsing and scalar guards only
├── errors/
├── io/
├── api/
└── advanced/
```

This list records ownership, not public import status or migration order.
Existing `science`, `workflows`, `contracts`, and `tables` packages may exist
temporarily as internal bridge structure during migration, but ADR-0054 does
not support them as external routes. They are not target owners of new domain
behaviour.

### Dependency rules

The steady-state logical graph must be acyclic.

The permitted capability directions are:

- Dataset model and durable-record leaves may depend on
  `preparation_records`, run lineage and allowed shared domains, but not on
  Preprocessing or Dataset builder/composition implementation;
- only Dataset construction, planning and builder composition may consume
  published Preprocessing boundaries;
- preprocessing may depend on `preparation_records`, shared scientific
  domains and technical foundations, but not on Datasets or run lineage;
- Differential may depend on dataset contracts/models, shared domains and the
  exact passive `preparation_records.policies`,
  `preparation_records.protein_aware` and
  `preparation_records.batch_correction` reader leaves accepted by D05A-1;
- Enrichment may depend on its input contracts and shared domains, but does
  not acquire an artificial dataset dependency;
- Kinase may depend on dataset contracts/models and shared domains;
- Kinase has no general dependency on Preprocessing; its two residual direct
  imports are governed exactly to Wave 7C;
- Signalome may depend on `KinaseWorkflowResult`, its nested published result
  contracts, shared domains, and the exact D08-1 Dataset symbols and consumers;
  no other Signalome-to-Datasets import is permitted; and
- public and IO adapters may depend inward on published capability contracts.

The following directions are forbidden:

- shared scientific domains → capability implementations;
- preprocessing or datasets → Differential, Enrichment, Kinase, or Signalome;
- general References → Kinase implementation;
- capability domains → `phospy.api` or `phospy.advanced`;
- capability/shared domains → concrete IO adapters; and
- unlisted capability → capability imports.

The shared-domain and technical-foundation graph is also explicit. An empty
entry means standard-library and third-party dependencies only.

| Source owner | May import from |
| --- | --- |
| `errors` | — |
| `taxonomy` | `errors` |
| `policies` | `errors` |
| `frames` | `errors` |
| root `validation` | `errors`, `frames`, `policies` |
| `provenance` | `errors`, `frames`, `taxonomy` |
| `activities` | `errors`, `frames`, `quantitative` |
| `phosphosites` | `errors`, `frames`, `policies`, `taxonomy`, root `validation` |
| `quantitative` | `errors`, `frames`, `policies`, `provenance`, root `validation` |
| `statistics` | `errors`, `frames`, `policies`, root `validation` |
| `references` | `errors`, `frames`, `phosphosites`, `policies`, `provenance`, `taxonomy`, root `validation` |
| `evidence` | `errors`, `frames`, `phosphosites`, `policies`, `provenance`, `quantitative`, `taxonomy`, root `validation` |
| `preparation_records` | `errors`, `frames`, `policies`, `provenance`, `quantitative`, root `validation` |
| `run_lineage` | `evidence`, `preprocessing`, `provenance`, `references` |

No shared owner may import Datasets, Preprocessing, Differential, Enrichment,
Kinase, Signalome, a public facade, or a concrete IO adapter.

`activities` is a passive shared scientific boundary. It must not import
Kinase, Signalome, Datasets, Preprocessing, Differential, Enrichment,
References, IO, either public facade, or another capability implementation.
Its permitted outbound owner dependencies are exactly `errors`, `frames` and
`quantitative`. Kinase and IO consumption of its published passive contracts
remains dependency flow into Activities and creates no additional outbound
permission.

`preparation_records` additionally may not import run lineage and may contain
no stage, planner, executor, builder, callback or registry. Evidence and
References deliberately do not import one another: preprocessing composes
reference retrieval with evidence resolution where both are required. The
matrix is an allow-list, not a requirement to create every listed edge.

`run_lineage` is a higher-level composition owner, not a shared scientific
foundation. The shared owners and Preprocessing therefore do not gain a
back-edge to it.

The accepted Wave 4C decision also moves the single
`KinaseLibraryResourceProvenance` definition early and narrowly to
`phospy.kinase.references.provenance`. This does not advance Kinase reference
projection, scoring, prediction, activity or another Wave 7 responsibility.

Permission for one capability-to-capability edge does not permit imports from
the complete target package. Architecture tests must restrict those imports to
named contract/result paths.

### Architecture enforcement

Architecture tests will enforce both package direction and path-level contract
use.

The tests must:

1. map implementation modules to logical owners;
2. reject cycles between logical owners;
3. compare observed owner edges with an explicit allowed-edge table;
4. reject shared-domain imports of capability implementations;
5. restrict Signalome → Kinase imports to published result contracts;
6. reject upstream imports from datasets/preprocessing/evidence/readers into
   Differential;
7. reject general References imports from Kinase scoring/prediction/activity;
8. reject capability imports from public facades or concrete IO adapters;
9. ensure concrete configs, validators, table schemas, and result types have
   one logical owner; and
10. permit facade modules to re-export exact identities without declaring new
    domain behaviour;
11. enforce the closed `preparation_records` leaf/symbol map, its lower-only
    imports, and the exact Differential and IO reader permissions; and
12. distinguish Dataset model/durable leaves from the legitimate Dataset
    construction/planning/builder composition role.
13. enforce the complete D05B-1/D05B-2 package-condensation classification set
    as three exact source-path/module/target-module/symbol imports, reject
    aliases, wildcards, dynamic imports, registration hooks and duplicates,
    and prove the complete eager-module graph is acyclic without excluding any
    of the three edges; and
14. prove D05B-1 remains the sole Dataset-to-IO acquisition edge and that its
    eager transitive import closure cannot return to Datasets.

Signalome-to-Datasets imports are an exact `(consumer module, defining module, symbol)` allow-list containing only the seven D08-1 triples. A source-level AST guard scans every module below `phospy.signalome`, derives its consumer from the source path, resolves absolute and relative imports, and follows local `phospy` re-exports to the defining symbol without importing production modules. Only a direct, absolute, unaliased spelling of an approved triple is accepted, once; imports beneath `TYPE_CHECKING` and other nested or conditional blocks are included.

The Dataset guard rejects an approved symbol from an unlisted consumer, an unlisted module or symbol, aliases, wildcards, whole-package/module imports, relative or re-export shortcuts, duplicates, Dataset-module attribute access, `importlib`/`__import__` acquisition, reflective or container/assignment-based module indirection, private Dataset storage acquisition, and every Dataset `_internal`, construction, builder or validation import. Static symbol/callable origin and binding propagation identify those forms; unresolved dynamic imports fail closed when derived from a Dataset prefix or used at a Dataset-sensitive sink. Enforcement remains scoped to Dataset-module- or Dataset-object-derived bindings and Dataset-shaped dynamic acquisition, so unrelated Signalome relative imports, statically resolved non-Dataset dynamic imports, ordinary object-property access, scientific attribute access and the same attribute spelling on a proven non-Dataset object remain legal. The full resolver rules and positive/negative fixture matrix are implementation obligations in ARCH-MIG-08; this paragraph does not widen D08-1's permission set.

Scientific regression, parity, validation, provenance, serialization,
dataframe-ownership, public import, and constructor tests remain separate from
architecture tests.

### Wave 5A D05A-1 passive preparation-record owner

The maintainer accepted D05A-1 on 5 October 2026. The neutral
`phospy.preparation_records` owner exists to preserve one concrete runtime
definition for passive values embedded by both Preprocessing outputs and
Dataset durable records without creating a Dataset/Preprocessing cycle.

Its declaration allow-list is closed by target leaf:

| Leaf | Owned declarations |
| --- | --- |
| `policies` | `MissingDataPolicy`, `NormalisationPolicy`, `SiteSequenceConflictPolicy`, `SiteSequenceResolutionMode`, `TotalProteinCorrectionPolicy`, and `TotalProteinCorrectionIdentityMatchingPolicy` |
| `imputation` | Persisted imputation input-scale sources, operation-order values, KNN no-overlap values, their sets, and schema version |
| `report_schema` | The seven Dataset-durable preprocessing row types, their column tuples, typed dataframe conversions, `reorder_columns`, and required generic dataclass/conversion helpers |
| `batch_correction` | `BatchCorrectionReport`, `BatchCorrectionPolicy`, `BatchCorrectionDiagnostics`, status/confounding/shape aliases, passive value constants, and their structural payload helpers |
| `protein_aware` | Passive preparation eligibility, sample-alignment, transformation-state, mapping-diagnostic and site-eligibility records; `ProteinMappingStatus`; their column/value constants and structural record helpers |
| `group_evidence` | The two group-aware mask hash functions and the existing mask-normalisation helper |
| `validation.sample_groups` | `ResolvedSampleGroups`, `SampleGroupMetadataResolver`, and its structural label/group validation helpers |

The accepted source-to-target symbol inventory in the Wave 5A preflight is
normative. Dependencies required for record construction, defaults,
annotations, intrinsic validation and payload methods move with the record.
An implementation must report an additional semantic dependency instead of
silently expanding this set.

Preprocessing retains caller-intent configuration, scientific policy
selection, planning, stage output, execution and active report assembly.
Datasets retains the analysis-ready model, durable state/reports/sidecars and
binding validation. In particular, `BatchCorrectionResult` and
`ProteinAwarePreparationStageOutput` remain Preprocessing-owned.

Differential's permanent access is restricted to the exact passive leaves
`policies`, `protein_aware` and `batch_correction` for the existing readers
named by D05A-1. IO's permanent access is restricted to `policies` for the
existing processing-state bundle reader. No other Differential or IO access
is implied. Public facades gain no new names or general access to this owner.

The accepted exception-register projection is 140 ordered records: 32 exact
deletions and six in-place Signalome target substitutions. The live register
does not change until the physical movement makes those closures true. This
ownership amendment does not itself authorise implementation.

### Wave 5B D05B-1 and D05B-2 package-condensation classifications

On 6 October 2026, maintainer Shane Martin accepted D05B-1 as defined in the
corrected Wave 5B preflight and separately authorised ARCH-MIG-05B
implementation. Shane Martin subsequently accepted D05B-2 as the following
two additional exact package-condensation classifications and declared D05B-1
plus D05B-2 to be the complete approved exclusion set:

```text
source path:   src/phospy/datasets/construction/input_resolution.py
source module: phospy.datasets.construction.input_resolution
target module: phospy.validation.common.input_sources
symbol:        validate_dataframe_or_path

source path:   src/phospy/datasets/validation/inputs.py
source module: phospy.datasets.validation.inputs
target module: phospy.validation.common.input_sources
symbol:        validate_dataframe_or_path
```

D05B-1 permits exactly this permanent builder-composition acquisition edge:

```text
source path:   src/phospy/datasets/builder.py
source module: phospy.datasets.builder
target module: phospy.io.readers.dataset_inputs
symbol:        DatasetPathTableReader
```

The D05B-1 import must use
`from phospy.io.readers.dataset_inputs import DatasetPathTableReader`. No
other Dataset leaf, IO target or IO symbol is permitted, including through an
alias, wildcard, dynamic import or registration hook. This permission is
classified and removed before logical-owner or package condensation; it does
not add a broad `datasets -> io` permission. The remaining condensed graph
and the complete resolved module graph must be acyclic.

Each D05B-2 import must use
`from phospy.validation.common.input_sources import
validate_dataframe_or_path`. D05B-2 grants no broader Dataset-to-Validation
import or condensation permission. D05B-1 and D05B-2 are permanent governed
classifications, not temporary architecture exceptions. They therefore add no
record to `temporary-exceptions.json`, whose accepted Wave 5B transformation
remains exactly 140 to 114 records. The package-condensed graph is checked
after removing exactly these three edges. The complete eager-module graph is
checked with all three edges present and must remain acyclic.

The root, `phospy.io`, `phospy.io.readers`, and every eager dependency loaded
by the permitted target must remain free of Dataset imports. Dataset package,
construction and validation initializers must not eagerly import the builder.
Opposing-order reachability tests must prove both that the acquisition path
cannot return to Datasets and that IO persistence/reconstruction paths cannot
reach the builder. IO remains the sole owner of filesystem behaviour.

The accepted Wave 5B live-register delta is exactly the 26 ordered deletions
listed in the corrected preflight, leaving 114 records with no re-key,
replacement or later-wave deletion.

## Public API separation

This ADR does not itself define supported public import paths. ADR-0054 now
defines the post-reset external surfaces, removes compatibility and semi-public
routes, and maps the retained aggregate exports to their logical owners.

A domain's target internal path is not automatically a supported public path.
Conversely, retaining a supported facade path does not make that facade the
implementation owner.

This ADR does not require compatibility shims. Any retained route must be
justified by the separate API policy and must preserve exact owner-defined
object identities without owning behaviour.

## Physical migration separation

This ADR does not choose migration order, temporary bridge modules, release
boundaries, or implementation tickets.

Migration planning may begin because this ADR and ADR-0054 are accepted. The
migration must preserve scientific behaviour and must introduce or update
architecture tests as dependency boundaries become enforceable.

## Consequences

### Positive consequences

- Each major capability becomes understandable from one bounded package.
- Configuration, validation, planning, execution, results, diagnostics, and
  provenance have coherent ownership.
- Legitimate shared scientific concepts remain explicit.
- Cross-capability dependencies become narrow and testable.
- The existing workflow lifecycle is retained rather than rewritten.
- Signalome continues to use the complete Kinase result without a duplicate
  transfer object.
- Upstream evidence and preprocessing no longer depend on Differential
  implementation details.
- Kinase-specific reference semantics no longer create reverse dependencies
  from general References.
- Architecture tests can detect logical cycles that the current broad
  `phospy.science` graph node hides.

### Negative consequences

- Many current modules have correct behaviour but the wrong logical owner,
  making migration broad even without scientific changes.
- The accepted API reset requires coordinated documentation, test, and typing
  updates before physical moves.
- Some current internal-view shortcuts, especially Signalome access to Kinase
  prediction/scoring state, require a replacement that preserves performance
  and dataframe ownership.
- During migration, facade and owner packages may temporarily coexist and make
  navigation more complex.
- Logical ownership tests require a maintained owner map until the physical
  package structure fully reflects the target design.

## Alternatives considered

### Retain the horizontal package architecture

This would preserve the current broad acyclic graph but would leave lifecycle
ownership distributed across contracts, validation, workflows, and science.
It would not address the hidden logical-domain cycles.

### Make every workflow dependency-free

This would duplicate dataset, identity, quantitative, reference, statistics,
and provenance concepts. Those concepts are genuinely shared and should not be
copied into each workflow.

### Force identical workflow directory structures

This would impose unnecessary files on compact capabilities and would hide the
real complexity differences between Enrichment, Kinase, preprocessing, and
Signalome. Common responsibilities and dependency rules are sufficient.

### Introduce a dedicated Signalome input DTO

This would duplicate information already governed by `KinaseWorkflowResult`
and add another result-projection contract. Signalome always follows Kinase, so
the complete Kinase result is the clearer boundary.

### Promote SPS discovery to a peer top-level capability

SPS discovery is independently runnable, but its result is consumed only by
preprocessing/batch correction. A nested preprocessing/control-site capability
expresses both its independent lifecycle and its actual ownership.

### Keep Experimental Design as a shared domain

No independent production consumer exists. Shared placement would preserve a
hypothetical reuse case rather than represent current ownership.

### Keep all quantification-depth definitions in Differential

Readers, peptide evidence, and preprocessing already create and validate depth
before Differential runs. Keeping the base contract in Differential preserves
an upstream-to-downstream dependency inversion. Only the Differential trend
interpretation belongs there.

## Validation and review criteria

This ADR is satisfied when:

1. every significant module has one documented logical owner;
2. the logical owner graph, including the shared-domain allow-list, is acyclic;
3. cross-capability imports are limited to approved contract/result paths;
4. Signalome accepts and retains `KinaseWorkflowResult` directly;
5. preprocessing and evidence do not import Differential;
6. general References does not import Kinase implementation;
7. `phosphosites` does not obtain `Organism` through References;
8. domain-specific configs, validators, tables, results, diagnostics, and
   provenance assembly are owned by their capability;
9. public facade tests follow ADR-0054 while scientific behaviour tests remain
   unchanged unless separately authorised; and
10. architecture enforcement covers both logical owner edges and restricted
    contract paths.

## Non-goals

This ADR does not:

- alter scientific algorithms, numerical outputs, eligibility, or validation
  policy;
- change `AnalysisReadyPhosphoDataset` invariants;
- change `site_key` or `display_id` semantics;
- change public API support or deprecation policy;
- specify compatibility shims;
- select a migration sequence or create implementation tickets;
- make Activity a peer workflow; or
- make SPS discovery a peer top-level analytical capability.

## Related records

- [ADR-0002: Internal Workflow Architecture](adr_0002_internal_workflow_architecture.md)
- [ADR-0003: Analysis-Ready Dataset and Preprocessing Boundary](adr_0003-analysis_ready_dataset_and_preprocessing_boundary.md)
- [ADR-0007: Validation Domain Architecture](adr_0007_validation_domain_architecture.md)
- [ADR-0010: Internal Package and Module Layout](adr_0010_internal_package_and_module_layout.md)
- [ADR-0019: Experimental Design, Contrast, and Replicate Contract](adr_0019_experimental_design_and_contrast_contract.md)
- [ADR-0020: Peptide Evidence and Site-Level Resolution Policy](adr_0020_peptide_evidence_and_site_level_resolution_policy.md)
- [ADR-0024: Protein-Scoped Phosphosite Row Identity](adr_0024_protein_scoped_phosphosite_row_identity.md)
- [ADR-0028: Semi-Public Science Import Policy](adr_0028_semi_public_science_import_policy.md)
- [ADR-0031: Public API Stability Tiers](adr_0031_public_api_stability_tiers.md)
- [ADR-0033: Result Caveats and Scientific Warning Ownership](adr_0033_result_caveats.md)
- [ADR-0034: Quantitative State, Motif Scoring Semantics, and Reference Context](adr_0034_quantitative_state_motif_semantics_and_reference_context.md)
- [ADR-0035: Provenance Immutability and Stable Serialization](adr_0035_provenance_immutability_and_stable_serialization.md)
- [ADR-0051: SPS Reference Quantitative State and Discovery Separation](adr_0051_sps_reference_quantitative_state_and_discovery_separation.md)
- [ADR-0052: Partial-Reference SPS Consensus Ordering](adr_0052_partial_reference_sps_consensus_ordering.md)
- ADR-0054: Public API Reset and Compatibility Removal

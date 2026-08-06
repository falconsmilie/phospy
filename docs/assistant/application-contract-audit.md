# PhosPy Assistant application-contract audit

Inspection date: 2026-08-06

## Audit scope and decision rules

This audit asks whether a separate PhosPy Assistant application can perform its
initial application-layer responsibilities by consuming the current PhosPy
source through supported public PhosPy APIs. It does not design or implement the
assistant, add an LLM provider, expose private validators, or redesign the
scientific API.

The public-contract decision rule used here is stricter than Python
importability:

- Supported public API means there is evidence that external consumers are
  expected to use the symbol. Evidence includes deliberate exports from
  `phospy`, `phospy.api`, `phospy.advanced`, or another documented public
  namespace; `__all__` membership; user-facing documentation; public API tests;
  or an ADR that classifies the symbol as public.
- Merely importable implementation code is not treated as public. Validators,
  workflow interpreters, executors, private result assemblers, and internal
  dataset-processing models remain private even when Python can import them.
- Deprecated compatibility imports are not treated as the preferred public
  contract when a non-deprecated public route exists.
- Test-only helpers are not part of the application contract.

The main public-status evidence is:

- `src/phospy/__init__.py`, which curates the top-level workflow convenience
  surface.
- `src/phospy/api/__init__.py` and `src/phospy/_api_inventory.py`, which define
  the stable public API inventory.
- `src/phospy/advanced/__init__.py`, which defines advanced supported imports.
- `docs/api/guide.md`, which documents the stable, advanced, and internal API
  tiers.
- `docs/adr/adr_0001_public_api_contract.md` and
  `docs/adr/adr_0031_public_api_stability_tiers.md`, which govern namespace
  ownership and public API stability.

Verdicts use these criteria:

- YES: the requirement can be satisfied through a supported, sufficiently typed
  and documented public PhosPy contract without direct imports from internal
  implementation modules. A small assistant-owned translation layer does not
  make the PhosPy contract inadequate.
- PARTIAL: a usable public path exists, but material limitations remain. These
  include ambiguous public status, undocumented schema guarantees, weakly typed
  structured data, no supported round-trip serialization, reliance on string
  parsing, or a stable operation that requires a narrowly scoped assistant
  adapter.
- NO: the requirement cannot be completed without importing internal
  implementation modules, depending on private validation objects, relying on
  unsupported implementation details, reconstructing unavailable scientific
  state, or making assumptions that are not verifiable from the public contract.

Ownership decisions use this rule:

- PhosPy owns scientific logic, scientific validation, dataset invariants,
  workflow validation composition, request/result contracts that normal Python
  consumers need, provenance, and result-table semantics.
- PhosPy Assistant owns LLM tool schemas, conversational state, natural-language
  query adapters, presentation DTOs, provider-specific request/response models,
  and small explicit adapters over stable PhosPy workflow contracts.
- A change belongs in PhosPy only when a normal external Python consumer
  reasonably needs a supported public contract that does not currently exist.

## Summary matrix

| Assistant requirement | Existing supported public PhosPy contract | Adequate? | Required change |
| --- | --- | --- | --- |
| Construct an analysis-ready dataset | `phospy.AnalysisReadyDatasetBuilder`, `phospy.api.AnalysisReadyDatasetBuilder`, `phospy.api.DatasetBuildRequest`, `phospy.api.AnalysisReadyPhosphoDataset` | YES | None |
| Receive structured validation outcomes | `phospy.api.PhosPyValidationError`, `phospy.api.WorkflowValidationError`, `phospy.api.WorkflowBoundaryError`, `phospy.api.ContractValidationError`, `phospy.api.PhosPyInputError` | PARTIAL | Assistant adapter over public exception classes now; narrow PhosPy diagnostic contract only if machine-readable issue collections are required |
| Discover supported workflows | `phospy.DifferentialAnalysisWorkflow`, `phospy.KinaseWorkflow`, `phospy.SignalomeWorkflow`, `phospy.api.EnrichmentWorkflow`, plus request/result exports in `phospy.api` | YES | Assistant-owned workflow-capability listing; no PhosPy registry required |
| Serialise requests | `phospy.api.DatasetBuildRequest`, `phospy.api.DifferentialAnalysisRequest`, `phospy.api.ExperimentalDesign`, `phospy.api.SampleDesignRecord`, `phospy.api.Contrast`; no supported request serialization path | PARTIAL | Assistant-owned serialization for first prototype; consider narrow PhosPy request round-trip contract only for general saved-analysis use |
| Run differential analysis | `phospy.DifferentialAnalysisWorkflow.run(request: phospy.api.DifferentialAnalysisRequest) -> phospy.api.DifferentialAnalysisResult`; also `phospy.api.DifferentialAnalysisWorkflow` | YES | None |
| Inspect typed attrition | `phospy.api.DifferentialAnalysisResult.table_for(...)`, `phospy.api.DifferentialAnalysisResult.feature_eligibility`, `phospy.api.DifferentialAnalysisResult.workflow_provenance`; `phospy.provenance.RowAttritionRecord` and `phospy.provenance.RowAttritionReport` are exported from the provenance namespace but not from `phospy.api` | PARTIAL | Assistant adapter over result tables/provenance payload; optional PhosPy documentation/typing clarification for row-attrition payload status |
| Export a rerunnable recipe | `phospy.api.DifferentialAnalysisRequest` records in-memory request intent; `phospy.api.DifferentialAnalysisResult.to_payload()` serializes result payloads; `phospy.provenance.to_payload()` serializes `RunProvenance`; no versioned recipe or request+dataset round-trip contract exists | NO | Deferred narrow PhosPy recipe contract if durable rerun export is required; otherwise assistant-owned session recipe only |
| Query results by site_key | `phospy.api.DifferentialAnalysisResult.table_for(contrast_name) -> pandas.DataFrame` with `site_key` index and `site_key` column | YES | Assistant-owned lookup/presentation adapter only |

## Detailed findings

### Construct an analysis-ready dataset

#### Requirement

An external application must construct a valid `AnalysisReadyPhosphoDataset`
from user-provided phosphoproteomics tables without importing internal builders
or validators. The construction path must preserve the current invariant that
analysis-ready site metadata includes `site_sequence`.

#### Existing public contract

Supported public symbols:

- `phospy.AnalysisReadyDatasetBuilder`
- `phospy.api.AnalysisReadyDatasetBuilder`
- `phospy.api.DatasetBuildRequest`
- `phospy.api.AnalysisReadyPhosphoDataset`
- Advanced/trusted lane: `phospy.api.AnalysisReadyPhosphoDataset.from_trusted_tables(...)`

Owning source files and methods:

- `src/phospy/__init__.py` re-exports `AnalysisReadyDatasetBuilder` and
  `AnalysisReadyPhosphoDataset` from the top-level package.
- `src/phospy/api/__init__.py` includes `AnalysisReadyDatasetBuilder`,
  `DatasetBuildRequest`, and `AnalysisReadyPhosphoDataset` in the stable
  public facade.
- `src/phospy/api/builders.py` defines the public
  `AnalysisReadyDatasetBuilder.run(request: DatasetBuildRequest) ->
  AnalysisReadyPhosphoDataset` wrapper with default reader, validator,
  interpreter, executor, preprocessing, and batch-correction wiring.
- `src/phospy/contracts/dataset_build.py` defines `DatasetBuildRequest` with
  public table/path inputs, site-resolution mode, preprocessing config, organism,
  intensity-scale declaration, and peptide-evidence fields.
- `src/phospy/science/datasets/construction/analysis_ready.py` implements
  `AnalysisReadyPhosphoDataset`. Its direct constructor is sealed; the class
  docstring directs ordinary users to `AnalysisReadyDatasetBuilder.run(...)`.

Important inputs and outputs:

- `DatasetBuildRequest.phospho`: pandas `DataFrame` or file path.
- `DatasetBuildRequest.site_metadata`: pandas `DataFrame` or file path for the
  site-level-resolved lane.
- `DatasetBuildRequest.site_resolution_mode`: default
  `"site_level_resolved"`; peptide-evidence lane is also represented.
- `DatasetBuildRequest.input_intensity_scale`: public declaration such as
  `"linear"` or `"log2"`.
- Output: `AnalysisReadyPhosphoDataset`, with defensive snapshot properties such
  as `phospho`, `site_metadata`, `sample_metadata`, `preprocessing_report`, and
  `provenance`.

#### Evidence

- `docs/api/guide.md` says ordinary dataset construction should use
  `AnalysisReadyDatasetBuilder().run(DatasetBuildRequest(...))` and describes
  `AnalysisReadyPhosphoDataset` as the strict analysis-ready boundary.
- `docs/api/dataset-build-workflow.md` documents builder usage.
- `docs/adr/adr_0031_public_api_stability_tiers.md` records
  `AnalysisReadyDatasetBuilder`, `AnalysisReadyPhosphoDataset`, and
  `DatasetBuildRequest` in the stable public API inventory. It also records
  that `AnalysisReadyPhosphoDataset.from_trusted_tables(...)` is an advanced
  trusted factory.
- `tests/unit/test_public_contract_dataset.py::test_builder_exposes_only_run_request_contract`
  asserts that the public builder exposes `run` and that its type hints are
  `DatasetBuildRequest -> AnalysisReadyPhosphoDataset`.
- `tests/unit/test_public_contract_dataset.py::test_public_dataset_ingestion_story_is_builder_only`
  asserts that the public dataset ingestion story is builder-only.
- `tests/unit/test_public_contract_errors.py::test_dataset_constructor_rejects_missing_site_sequence_column`
  verifies that the dataset boundary rejects missing `site_sequence`.

#### Behavioural assessment

An external application can pass pandas tables or file paths to
`DatasetBuildRequest`, run `AnalysisReadyDatasetBuilder().run(request)`, and
receive an `AnalysisReadyPhosphoDataset` whose `phospho` and `site_metadata`
tables are keyed by `site_key`. The builder owns user-input interpretation,
preprocessing, private validation, processing-state establishment, and
provenance. The external application does not need to import dataset validators
or implementation builders.

The direct dataset constructor raises immediately. That is deliberate and does
not weaken the public path because the builder is the ordinary construction API.
The trusted-table factory exists for advanced callers with complete assertions
and provenance, but it is not the preferred assistant path.

Forced internal import: none.

#### Adequacy verdict

YES.

The public builder/request/domain contract is explicit, documented, typed, and
tested. It preserves `site_sequence` and does not require exposing private
dataset validation.

#### Ownership decision

No change required.

#### Minimal required change

None.

### Receive structured validation outcomes

#### Requirement

An external application must receive validation failures in a consumable form
when dataset construction or workflow execution fails. The relevant question is
whether outcomes are consumable, not whether private validator implementations
are public.

#### Existing public contract

Supported public symbols:

- `phospy.api.PhosPyValidationError`
- `phospy.api.WorkflowValidationError`
- `phospy.api.WorkflowBoundaryError`
- `phospy.api.ContractValidationError`
- `phospy.api.PhosPyInputError`
- `phospy.api.UnsupportedInputFormatError`

Related non-stable/ambiguous symbols:

- `phospy.errors.DatasetValidationError` is exported by `phospy.errors`, but
  ADR-0031 classifies `DatasetValidationError` as internal/experimental rather
  than stable `phospy.api`.
- Validator classes under `phospy.validation.*` and workflow validators under
  `phospy.workflows.*.validator` are not public contracts.

Owning source files and methods:

- `src/phospy/api/__init__.py` re-exports the common stable exception families.
- `src/phospy/errors/validation.py` defines validation exception classes.
- `src/phospy/errors/workflows.py` defines `WorkflowBoundaryError`, whose
  structured attributes are `seam`, `next_action`, and `details`.
- `src/phospy/errors/__init__.py` exports the broader exception taxonomy.

#### Evidence

- `docs/api/guide.md` documents catching public exception families such as
  `WorkflowValidationError` and `PhosPyValidationError`.
- `docs/adr/adr_0031_public_api_stability_tiers.md` lists common exception
  families in the stable public API and classifies validator internals and
  dataset-processing internals as internal/experimental.
- `tests/unit/test_public_contract_errors.py::test_top_level_exception_exports_match_curated_facade`
  verifies the curated public exception facade and distinguishes non-facade
  errors.
- `tests/unit/test_public_contract_boundary_honesty.py` and
  `tests/unit/test_public_contract_dataset.py` verify many failures through
  public exceptions rather than public validator objects.

#### Behavioural assessment

An external application can catch stable public exception families and route
failures by broad type:

- `ContractValidationError` for value-object construction validation.
- `PhosPyInputError` and `UnsupportedInputFormatError` for input/table contract
  failures.
- `PhosPyValidationError` for validation failures as a broad catch-all.
- `WorkflowValidationError` for workflow-level validation failures.
- `WorkflowBoundaryError` for boundary failures with structured `seam`,
  `next_action`, and `details`.

Most dataset and workflow validation failures are not exposed as structured
issue collections with stable error codes, field paths, column names, and
workflow context. Many tests assert substrings in exception messages, which is
good for human diagnostics but not a stable machine-readable issue model.
Importing private validators would provide implementation detail, not a
supported application contract.

Forced internal import: no internal import is required to catch public
exceptions. An application that requires private validator issue objects or
validator-specific return values would be forced into internal modules because
no public issue-collection contract currently exists.

Classification of any forced internal import: genuine missing general-purpose
public contract only if the consumer needs machine-readable issue collections;
otherwise it is safely handled by an assistant-owned exception adapter that
uses public exception classes and messages.

#### Adequacy verdict

PARTIAL.

Typed exception classes are public and adequate for coarse failure handling.
Only `WorkflowBoundaryError` exposes stable structured diagnostic fields. A
complete structured validation-outcome contract with issue codes and field or
column context is not currently public.

#### Ownership decision

Assistant-owned adapter for the first assistant increment; narrow PhosPy
public-contract addition only if machine-readable issue collections are required
by normal external Python consumers.

The assistant can start by catching public exceptions and presenting their
human-readable messages without importing validators. PhosPy should not expose
validator implementations merely for assistant integration.

#### Minimal required change

None for a first prototype that displays public exception messages. If durable
machine-readable validation outcomes become a general consumer requirement, add
a narrow public diagnostic payload contract to public exception classes rather
than exporting validators.

### Discover supported workflows

#### Requirement

An external application must identify the workflows PhosPy supports and their
stable public entrypoints/identifiers.

#### Existing public contract

Supported public symbols:

- `phospy.DifferentialAnalysisWorkflow`
- `phospy.KinaseWorkflow`
- `phospy.SignalomeWorkflow`
- `phospy.api.DifferentialAnalysisWorkflow`
- `phospy.api.EnrichmentWorkflow`
- `phospy.api.KinaseWorkflow`
- `phospy.api.SignalomeWorkflow`
- Request/result pairs:
  `phospy.api.DifferentialAnalysisRequest`,
  `phospy.api.DifferentialAnalysisResult`,
  `phospy.api.EnrichmentWorkflowRequest`,
  `phospy.api.EnrichmentWorkflowResult`,
  `phospy.api.KinaseWorkflowRequest`,
  `phospy.api.KinaseWorkflowResult`,
  `phospy.api.SignalomeWorkflowRequest`,
  `phospy.api.SignalomeWorkflowResult`

Owning source files and methods:

- `src/phospy/__init__.py` exports the top-level workflow convenience surface.
- `src/phospy/api/workflows.py` exports all public workflow classes.
- `src/phospy/api/requests.py` exports public request DTOs.
- `src/phospy/api/results.py` exports primary result objects.
- Each public workflow class exposes `run(...)` as the public method.

#### Evidence

- `docs/api/guide.md` contains a workflow map with Dataset builder,
  Differential analysis, Enrichment, Kinase, and Signalome pages and explicitly
  documents the request/result shape for each workflow.
- `docs/workflow_contracts.md` contains a workflow contract table listing
  workflow, public entrypoint, request, and result.
- `docs/adr/adr_0001_public_api_contract.md` defines workflow-oriented public
  API shape.
- `docs/adr/adr_0031_public_api_stability_tiers.md` lists workflow classes and
  request/result objects in the stable public API inventory.
- `tests/unit/test_public_contract_workflows.py::test_public_workflow_and_request_exports_match_contract`
  and `tests/unit/test_public_contract_workflows.py::test_public_workflows_expose_run_only`
  verify the public workflow exports and the `run`-only method surface.

#### Behavioural assessment

An external application can discover public workflows by inspecting documented
public exports (`phospy.__all__`, `phospy.api.workflows.__all__`, and
`phospy.api.__all__`) and docs. The stable identifiers are the public class and
request/result names. There is no dynamic workflow registry or enum, but none is
necessary for the assistant's first increment because the supported workflow set
is small and explicitly documented.

Forced internal import: none.

#### Adequacy verdict

YES.

The public workflow set is intentionally exported, documented, and tested. A
dynamic PhosPy registry would be an unnecessary abstraction for the current
assistant boundary.

#### Ownership decision

Assistant-owned adapter.

The assistant can maintain a small explicit mapping from its own capability
names to these stable workflow classes and request models. That mapping is
assistant presentation/orchestration, not PhosPy scientific API.

#### Minimal required change

None.

### Serialise requests

#### Requirement

An external application must serialize workflow requests, with detailed
attention to the differential-analysis request. The requirement includes
deterministic structured data, JSON compatibility, round-trip restoration,
policy handling, version identification, and rejection of arbitrary or invalid
scientific state.

#### Existing public contract

Supported public request symbols:

- `phospy.api.DatasetBuildRequest`
- `phospy.api.DifferentialAnalysisRequest`
- `phospy.api.ExperimentalDesign`
- `phospy.api.SampleDesignRecord`
- `phospy.api.Contrast`
- `phospy.api.FixedEffectCovariate`
- `phospy.api.CategoricalCovariate`
- `phospy.api.ContinuousCovariate`
- `phospy.api.BatchCovariate`
- Advanced differential config symbols when explicit non-default policy is
  needed: `phospy.advanced.DifferentialAnalysisConfig`,
  `phospy.advanced.EmpiricalBayesConfig`, and
  `phospy.advanced.MultipleTestingConfig`.

No supported PhosPy request serialization path was found for
`DifferentialAnalysisRequest` or `DatasetBuildRequest`.

Owning source files and methods:

- `src/phospy/contracts/requests.py` defines public workflow request dataclasses.
- `src/phospy/contracts/dataset_build.py` defines `DatasetBuildRequest`.
- `src/phospy/science/design/models.py` defines `ExperimentalDesign`,
  `SampleDesignRecord`, and `Contrast`.
- `src/phospy/contracts/configs/differential.py` defines differential config
  objects.

#### Evidence

- `src/phospy/contracts/requests.py` states that request dataclasses are
  lightweight command payloads and intentionally do not perform scientific
  validation during construction.
- `src/phospy/contracts/requests.py::DifferentialAnalysisRequest` has fields
  `dataset: AnalysisReadyPhosphoDataset`, `design: ExperimentalDesign`,
  `contrasts: tuple[Contrast, ...]`, and `config:
  DifferentialAnalysisConfig`.
- `tests/unit/test_public_contract_workflows.py::test_workflow_requests_keep_ingestion_outside_workflows`
  verifies these type contracts.
- Repository search found `to_payload`/`from_payload` on result and provenance
  models, but not on public workflow request models.
- `docs/api/guide.md` documents request objects as public command payloads and
  states that validation occurs when builders/workflows are run.

#### Behavioural assessment

An external application can construct typed request objects and can inspect
their dataclass fields in memory. Constructor-time validation exists for local
value-object invariants such as non-empty sample IDs and contrast names.
Workflow-level scientific validation runs at `run(...)`.

However, the public request models do not provide supported `to_payload`,
`from_payload`, JSON schema, schema version, or round-trip restoration methods.
Generic `dataclasses.asdict` is not a PhosPy-supported serialization contract
and is insufficient for request objects that contain pandas DataFrames,
datasets, reference objects, and config policy classes. For differential
analysis specifically, serializing the request also requires a supported way to
refer to or snapshot the input `AnalysisReadyPhosphoDataset`.

Forced internal import: none for in-memory request construction. No internal
import can safely solve supported request round-trip serialization; using
interpreter/executor internals would be unsupported and would not provide a
public saved-request schema.

Classification of any forced internal import: genuine missing general-purpose
public contract only if PhosPy chooses to support saved request interchange.
Assistant-specific JSON/tool formats should remain in the assistant.

#### Adequacy verdict

PARTIAL.

Typed public request models are adequate for in-memory execution. They are not
adequate as a supported request serialization and round-trip contract.

#### Ownership decision

Assistant-owned adapter for the first prototype; possible narrow PhosPy
public-contract addition for general saved-analysis use.

The assistant can serialize its own high-level inputs and rehydrate public
PhosPy request objects inside its application boundary. PhosPy should only add
request payload helpers if non-assistant consumers need a supported
request-interchange format.

#### Minimal required change

None for in-process assistant execution. For durable request interchange, add
narrow, versioned `to_payload`/`from_payload` helpers for public request
contracts and explicit dataset-reference/snapshot semantics. Do not introduce
provider-specific DTOs.

### Run differential analysis

#### Requirement

An external application must construct a dataset, construct a public
differential request, invoke the supported public workflow entrypoint, and
receive the public result type.

#### Existing public contract

Supported public symbols:

- `phospy.DifferentialAnalysisWorkflow`
- `phospy.api.DifferentialAnalysisWorkflow`
- `phospy.api.DifferentialAnalysisRequest`
- `phospy.api.ExperimentalDesign`
- `phospy.api.SampleDesignRecord`
- `phospy.api.Contrast`
- `phospy.api.DifferentialAnalysisResult`
- Optional advanced config symbols:
  `phospy.advanced.DifferentialAnalysisConfig`,
  `phospy.advanced.EmpiricalBayesConfig`,
  `phospy.advanced.MultipleTestingConfig`

Owning source files and methods:

- `src/phospy/workflows/differential/public.py` defines
  `DifferentialAnalysisWorkflow.run(request: DifferentialAnalysisRequest) ->
  DifferentialAnalysisResult`.
- `src/phospy/api/workflows.py` exports the workflow class.
- `src/phospy/api/requests.py` exports `DifferentialAnalysisRequest` and design
  primitives.
- `src/phospy/api/results.py` exports `DifferentialAnalysisResult`.

#### Evidence

- `docs/api/differential-analysis.md` documents
  `DifferentialAnalysisWorkflow.run(...)` returning `DifferentialAnalysisResult`.
- `docs/api/guide.md` includes the same workflow map.
- `docs/adr/adr_0001_public_api_contract.md` records top-level
  `DifferentialAnalysisWorkflow` compatibility expectations.
- `tests/unit/test_public_contract_workflows.py::test_workflow_run_type_contracts_are_request_to_result`
  verifies the run type hints.
- `tests/integration/test_differential_workflow_integration.py::test_differential_workflow_runs_on_builder_log2_dataset`
  demonstrates builder-created dataset execution, although that existing file
  imports internal helpers elsewhere and therefore did not satisfy the strict
  public-consumer verification requirement for this ticket.
- `tests/contract/test_public_differential_consumer.py::test_public_consumer_builds_and_runs_differential_by_site_key`
  was added by this audit as the focused public-consumer verification.

#### Behavioural assessment

An external application can create a valid analysis-ready dataset through the
public builder, construct a `DifferentialAnalysisRequest` with explicit
`ExperimentalDesign` and `Contrast` objects, call
`DifferentialAnalysisWorkflow().run(request)`, and receive a
`DifferentialAnalysisResult`. The supported convention is `run(...)`; public
workflow classes intentionally do not expose `execute`,
`run_from_analysis_ready`, or direct interpreter/executor hooks.

Forced internal import: none.

#### Adequacy verdict

YES.

The differential workflow entrypoint is public, typed, documented, and covered
by public-contract tests.

#### Ownership decision

No change required.

#### Minimal required change

None.

### Inspect typed attrition

#### Requirement

An external application must inspect row exclusions, filtering, or attrition
without parsing logs. Relevant information includes stage/reason, row counts,
affected identifiers where supported, provenance, and typed records or
well-defined structured payloads.

#### Existing public contract

Supported or externally visible symbols and surfaces:

- `phospy.api.DifferentialAnalysisResult`
- `phospy.api.DifferentialAnalysisResult.table_for(contrast_name)`
- `phospy.api.DifferentialAnalysisResult.feature_eligibility`
- `phospy.api.DifferentialAnalysisResult.workflow_provenance`
- `phospy.provenance.RowAttritionRecord`
- `phospy.provenance.RowAttritionReport`

Public-status caveat:

- `RowAttritionRecord` and `RowAttritionReport` are exported by
  `phospy.provenance.__all__`, and `docs/workflow_contracts.md` describes typed
  `row_attrition` reports. They are not included in the stable `phospy.api`
  inventory, so their status is supported by provenance-namespace evidence but
  weaker than the stable workflow/request/result API.

Owning source files and methods:

- `src/phospy/science/differential/models/results.py` defines
  `DifferentialAnalysisResult`, `table_for(...)`, `feature_eligibility`, and
  `to_payload()`.
- `src/phospy/science/differential/models/tables.py` defines result-status
  columns and validates differential result tables.
- `src/phospy/provenance/models/tables.py` defines
  `RowAttritionRecord.to_payload()`, `RowAttritionReport.from_records(...)`, and
  `RowAttritionReport.to_payload()`.
- `src/phospy/workflows/differential/provenance.py` assembles
  `row_attrition_metrics` and `row_attrition` payloads for differential
  workflow provenance.

#### Evidence

- `docs/workflow_contracts.md` states that typed `row_attrition` reports are
  causal site-row provenance and explains the semantics of `row_attrition` and
  `row_attrition_metrics`.
- `docs/api/differential-analysis.md` documents `workflow_provenance`,
  `policy_provenance`, `input_dataset_preprocessing_report`, and result table
  access.
- `tests/unit/test_workflow_stage_decomposition.py::test_differential_provenance_assembler_records_row_attrition`
  verifies differential row-attrition metrics and row-attrition records.
- `tests/unit/test_differential_analysis.py::test_differential_analysis_withholds_all_constant_site_intensities`
  verifies public result-table row status behaviour for withheld rows.
- `tests/unit/test_differential_result_contract.py` verifies the public
  differential result-table identity contract and `site_key` handling.

#### Behavioural assessment

An external application can inspect:

- per-contrast result tables through `table_for(...)`;
- `feature_eligibility` when present, with `site_key`, `result_status`, and
  `result_status_reason`;
- `workflow_provenance["row_attrition_metrics"]` and
  `workflow_provenance["row_attrition"]` when attrition occurred;
- typed provenance model classes from `phospy.provenance` if it chooses to
  validate or construct row-attrition payloads deliberately.

This is sufficient to present attrition without log parsing. The limitation is
that `DifferentialAnalysisResult.workflow_provenance` is typed as a generic
`Mapping[str, object]`; the differential result does not expose a stable
`row_attrition` property returning `RowAttritionReport`, and the row-status
constants are not stable `phospy.api` exports. The assistant can still
normalize the existing public table/provenance payload into its own display
format.

Forced internal import: none if consuming result tables and provenance payloads.
Importing `phospy.science.differential.models.tables` solely for status
constants would be a documentation/export ambiguity and should be avoided by an
assistant application.

#### Adequacy verdict

PARTIAL.

The data exists in structured public result/provenance surfaces, but the stable
API does not expose a single strongly typed differential attrition accessor.

#### Ownership decision

Assistant-owned adapter, with optional PhosPy documentation or typing
clarification.

The assistant can transform the existing table/provenance structures for
presentation. PhosPy should only add or promote a typed attrition accessor if
normal external consumers need to depend on it directly.

#### Minimal required change

No production change for the first assistant increment. A small PhosPy
documentation/typing clarification could state whether
`workflow_provenance["row_attrition"]` is a supported payload schema for
external consumers.

### Export a rerunnable recipe

#### Requirement

An external application must export a rerunnable analysis recipe sufficient to
reconstruct an equivalent workflow invocation. This must be distinguished from
ordinary provenance or result serialization.

#### Existing public contract

Available public or externally visible contracts:

- `phospy.api.DifferentialAnalysisRequest` records the in-memory request object.
- `phospy.api.DatasetBuildRequest` records in-memory dataset-build intent.
- `phospy.api.DifferentialAnalysisResult.to_payload()` serializes result tables,
  caveats, diagnostics, and workflow provenance to JSON-compatible Python data.
- `phospy.provenance.RunProvenance`
- `phospy.provenance.to_payload(provenance: RunProvenance)`
- `phospy.provenance.from_payload(payload) -> RunProvenance`

No qualifying public symbol was found for a versioned rerunnable recipe model,
recipe export, request+dataset snapshot export, or recipe import/rerun.

Owning source files and methods:

- `src/phospy/science/differential/models/results.py::DifferentialAnalysisResult.to_payload`
  serializes result payloads.
- `src/phospy/provenance/models/workflows.py::RunProvenance` defines run
  provenance.
- `src/phospy/provenance/serialization/workflows.py::to_payload` and
  `from_payload` serialize current run provenance.
- `src/phospy/contracts/requests.py` and `src/phospy/contracts/dataset_build.py`
  define in-memory request objects but no recipe or request payload methods.

#### Evidence

- `docs/workflow_contracts.md` states that workflow provenance records resolved
  request/config choices, table fingerprints, policies, diagnostics, and
  environment metadata, but also states that saved workflow bundles and
  provenance payloads are supported only for the current PhosPy schema.
- `docs/api/differential-analysis.md` documents `result.to_payload()` for
  result handoff and provenance preservation, not recipe reconstruction.
- `docs/adr/adr_0035_provenance_immutability_and_stable_serialization.md`
  governs provenance immutability and serialization semantics; it does not
  define a rerunnable recipe contract.
- Repository search did not find a public recipe model or public request
  round-trip serializer.

#### Behavioural assessment

Current PhosPy provenance and result payloads are useful for audit,
fingerprinting, diagnostics, and result handoff. They are not a recipe because
they do not provide a supported way to reconstruct:

- the original `DatasetBuildRequest`;
- the dataset tables or an immutable data snapshot;
- an equivalent `AnalysisReadyPhosphoDataset`;
- the differential `DifferentialAnalysisRequest`;
- request schema version;
- random/determinism settings beyond provenance metadata;
- all external reference/data identities needed for a rerun.

An assistant can store its own session-local inputs and rebuild public request
objects while those inputs remain available. That is an assistant recipe/log,
not a PhosPy-supported rerunnable recipe.

Forced internal import: no internal import can provide a supported recipe.
Using provenance internals or workflow interpreters to reconstruct requests
would rely on unsupported implementation details and incomplete state.

Classification of any forced internal import: genuine missing general-purpose
public PhosPy contract if durable rerunnable recipe export is a PhosPy feature;
deferred capability if the first assistant prototype can proceed without
durable recipe export.

#### Adequacy verdict

NO.

The repository supports result and provenance serialization, but not a
rerunnable analysis recipe.

#### Ownership decision

Deferred capability not required for the first assistant increment, unless the
first increment explicitly requires durable replay outside the current process.
If durable recipe export is required, it belongs in PhosPy as a narrow public
recipe contract because it must be scientifically faithful to PhosPy dataset,
request, provenance, and reference semantics.

#### Minimal required change

Do not implement a recipe in this ticket. A future change should define a
versioned recipe contract that states exactly how dataset snapshots or dataset
identities, public requests, configs, reference identities, provenance, and
determinism metadata are captured and restored.

### Query results by site_key

#### Requirement

An external application must reliably retrieve one or more differential result
rows by `site_key`, distinguish `site_key` from `display_id`, and handle missing
or duplicate display labels without importing internal result-table helpers.

#### Existing public contract

Supported public symbols:

- `phospy.api.DifferentialAnalysisResult`
- `phospy.api.DifferentialAnalysisResult.table_for(contrast_name) -> pandas.DataFrame`
- Advanced table helpers for reporting only:
  `phospy.advanced.filter_differential_results` and
  `phospy.advanced.rank_differential_results`

Owning source files and methods:

- `src/phospy/science/differential/models/results.py` defines
  `DifferentialAnalysisResult.table_for(...)` and validates contrast table
  ownership.
- `src/phospy/science/differential/models/tables.py` enforces the differential
  result table contract, including required identity columns and unique string
  indexes.
- `src/phospy/science/sites/identity_rules/result_identity.py` enforces
  result-table identity coherence.
- `src/phospy/science/tables/differential.py` defines advanced filter/rank
  helpers re-exported from `phospy.advanced`.

#### Evidence

- `docs/api/differential-analysis.md` states that each contrast result table is
  indexed by the input `site_key`, includes `site_key` and `display_id`, and
  tells users to interpret rows by `site_key`.
- `docs/workflow_contracts.md` states that site-level workflow outputs that
  materialize identity include both `site_key` and `display_id`.
- `README.md` states that `site_key` is true analysis-ready phosphosite row
  identity and `display_id` is a human-readable label.
- `tests/unit/test_differential_result_contract.py::test_result_tables_follow_public_differential_contract`
  verifies that result tables are indexed by `site_key` and include identity
  columns.
- `tests/unit/test_differential_result_contract.py::test_workflow_keeps_duplicate_display_ids_with_distinct_site_keys`
  verifies duplicate `display_id` values do not collide when `site_key` values
  differ.
- `tests/contract/test_public_differential_consumer.py::test_public_consumer_builds_and_runs_differential_by_site_key`
  verifies public result lookup by `site_key`.

#### Behavioural assessment

An external application can call `result.table_for("contrast_name")`, confirm
that the returned table index name is `site_key`, and retrieve rows with normal
pandas operations such as `table.loc[[site_key]]`. The result table includes a
`site_key` column that is required to match the index, so consumers can preserve
identity in exports. Missing keys naturally raise `KeyError`, and the index is
validated as unique. Duplicate `display_id` values are valid and do not affect
site-key lookup.

No PhosPy-specific query API is required for exact-key lookup. A restricted
natural-language result query adapter belongs in PhosPy Assistant.

Forced internal import: none.

#### Adequacy verdict

YES.

Direct filtering of the documented public result table is adequate for PhosPy's
contract. Assistant-specific lookup safeguards and UI phrasing can live outside
PhosPy.

#### Ownership decision

Assistant-owned adapter for presentation only.

#### Minimal required change

None.

## Lightweight public-consumer verification

I first inspected existing unit, integration, smoke, and architecture tests.
Relevant existing tests include:

- `tests/unit/test_public_contract_workflows.py::test_workflow_run_type_contracts_are_request_to_result`
- `tests/unit/test_public_contract_dataset.py::test_builder_exposes_only_run_request_contract`
- `tests/integration/test_differential_workflow_integration.py::test_differential_workflow_runs_on_builder_log2_dataset`
- `tests/integration/test_differential_workflow_integration.py::test_documented_two_vs_two_differential_example_contract`
- `tests/integration/test_dataset_batch_correction_validation.py::test_dataset_batch_correction_downstream_differential_uses_corrected_matrix`

No existing test satisfied all required public-consumer criteria. The runnable
differential tests either import internal modules or test-support factories at
module scope, derive expected `site_key` values through internal/test helper
routes, or do not assert the full public-consumer boundary.

Therefore this audit adds exactly one focused test:

- `tests/contract/test_public_differential_consumer.py::test_public_consumer_builds_and_runs_differential_by_site_key`

That test is self-contained and imports PhosPy only through supported public
routes:

- `phospy.AnalysisReadyDatasetBuilder`
- `phospy.DifferentialAnalysisWorkflow`
- `phospy.api.AnalysisReadyPhosphoDataset`
- `phospy.api.DatasetBuildRequest`
- `phospy.api.DifferentialAnalysisRequest`
- `phospy.api.DifferentialAnalysisResult`
- `phospy.api.ExperimentalDesign`
- `phospy.api.SampleDesignRecord`
- `phospy.api.Contrast`
- `phospy.api.Organism`

It constructs a deterministic synthetic dataset with `site_sequence`, builds an
analysis-ready dataset through the public builder, constructs the public
differential request, executes `DifferentialAnalysisWorkflow.run(...)`, verifies
the public result type, and confirms that the result table preserves and can be
looked up by `site_key`.

## Forced internal imports by requirement

| Assistant requirement | Forced internal import? | Classification |
| --- | --- | --- |
| Construct an analysis-ready dataset | No | Public builder and request are sufficient. |
| Receive structured validation outcomes | No for public exception handling; yes only if private validator issue objects are required | Genuine missing public diagnostic contract only for machine-readable issue collections; otherwise assistant adapter. |
| Discover supported workflows | No | Public exports/docs are sufficient; assistant-owned capability listing. |
| Serialise requests | No internal import helps | Genuine missing request round-trip contract if PhosPy supports saved requests; assistant-specific serialization should stay external. |
| Run differential analysis | No | Public workflow/request/result path is sufficient. |
| Inspect typed attrition | No if consuming public result/provenance payload; importing internal status constants would be avoidable | Documentation/export ambiguity; assistant adapter can normalize payloads. |
| Export a rerunnable recipe | No internal import helps | Genuine missing recipe contract if durable rerun export is required; otherwise deferred. |
| Query results by site_key | No | Public result table is sufficient. |

## Recommended follow-up tickets

Only the tickets below are recommended. Assistant-only concerns such as LLM
providers, prompts, tool-call schemas, conversational state, natural-language
querying, and presentation DTOs remain out of PhosPy.

### Ticket 1: Define public validation diagnostic payloads for boundary failures

Problem statement:

Public exception classes let external applications distinguish broad failure
families, and `WorkflowBoundaryError` exposes `seam`, `next_action`, and
`details`. Most dataset and workflow validation failures remain message-only
from a stable public-contract perspective. External consumers that need
machine-readable validation outcomes would otherwise need to parse messages or
import private validators.

Exact public contract affected:

- `phospy.api.PhosPyValidationError`
- `phospy.api.WorkflowValidationError`
- `phospy.api.ContractValidationError`
- `phospy.api.PhosPyInputError`
- `phospy.api.WorkflowBoundaryError`

Proposed ownership:

PhosPy, because validation failure semantics are owned by the scientific
package and are useful to normal Python consumers beyond the assistant.

Implementation scope:

- Define a minimal stable diagnostic payload shape on public exception types, or
  document explicitly which exception attributes are stable.
- Include code/category, message, optional field or column context, optional
  workflow stage/context, and optional next action.
- Preserve current exception hierarchy and messages.
- Do not expose validator classes or validation-domain internals.

Acceptance criteria:

- Public docs explain the stable diagnostic payload.
- Public API tests cover at least one dataset/build failure and one workflow
  failure.
- Private validators remain unexported from `phospy.api`.
- Existing exception catching remains backward compatible.

Non-goals:

- No LLM-specific validation DTOs.
- No public validator registry.
- No exposure of private validation functions/classes.
- No scientific validation behaviour changes.

### Ticket 2: Define narrow request payload round-trip helpers for public workflow requests

Problem statement:

Public request dataclasses are adequate for in-memory execution but do not
provide supported JSON-compatible serialization or restoration. External
consumers needing saved request interchange must currently write their own
schema and make unsupported assumptions for datasets, pandas tables, configs,
and policy values.

Exact public contract affected:

- `phospy.api.DatasetBuildRequest`
- `phospy.api.DifferentialAnalysisRequest`
- `phospy.api.ExperimentalDesign`
- `phospy.api.SampleDesignRecord`
- `phospy.api.Contrast`
- Advanced differential config classes as needed through `phospy.advanced`

Proposed ownership:

PhosPy, only for a general saved-analysis/request-interchange feature.
Assistant-specific tool schemas remain assistant-owned.

Implementation scope:

- Add narrow, versioned, JSON-compatible `to_payload`/`from_payload` helpers or
  equivalent functions for public request contracts.
- Explicitly define how pandas inputs, dataset references/snapshots, enums, and
  advanced config policy values are represented.
- Reject unsupported or scientifically invalid restored state through existing
  public constructor/workflow validation.

Acceptance criteria:

- Round-trip tests for a minimal differential request.
- Payload schema includes a version field.
- No provider-specific fields.
- No workflow interpreter/executor exposure.
- Existing request constructors remain usable.

Non-goals:

- No assistant facade.
- No LLM tool schema.
- No arbitrary dataframe execution.
- No change to differential calculations or validation rules.

### Ticket 3: Define a versioned rerunnable analysis recipe contract

Problem statement:

Result payloads and run provenance support audit and handoff, but they are not
a rerunnable recipe. A durable recipe must specify how to reconstruct an
equivalent dataset and workflow request, including dataset snapshot or identity,
request/config payload, reference identity, transformation state, provenance,
and deterministic settings.

Exact public contract affected:

- New narrow public recipe contract, if accepted.
- Existing related contracts:
  `phospy.api.DatasetBuildRequest`,
  `phospy.api.DifferentialAnalysisRequest`,
  `phospy.api.AnalysisReadyPhosphoDataset`,
  `phospy.api.DifferentialAnalysisResult`,
  `phospy.provenance.RunProvenance`.

Proposed ownership:

PhosPy, because recipe correctness depends on PhosPy scientific dataset,
workflow, provenance, and reference semantics.

Implementation scope:

- Define a versioned recipe model or functions for export/import.
- Specify dataset snapshot versus dataset identity semantics.
- Specify request/config payloads and compatibility policy.
- Specify reference identity/version and deterministic settings.
- Provide a rerun path that reconstructs public request objects and invokes the
  public workflow entrypoint.

Acceptance criteria:

- A recipe produced from a minimal public differential workflow can be restored
  and rerun to an equivalent public result under documented determinism limits.
- The recipe records schema version and PhosPy version.
- Unsupported or incomplete recipes fail with public exceptions.
- Tests cover missing data snapshot/reference identity failures.

Non-goals:

- No implementation in this audit ticket.
- No LLM integration.
- No prompt/tool schema generation.
- No broad workflow facade.
- No guarantee of bitwise identical outputs across environments beyond
  documented determinism limits.

### Ticket 4: Clarify differential row-attrition public payload status

Problem statement:

Differential row attrition is available through result tables and
`workflow_provenance`, and typed `RowAttritionRecord`/`RowAttritionReport`
models exist in `phospy.provenance`. The stable `phospy.api` facade does not
expose a typed differential `row_attrition` property, and the public status of
the `workflow_provenance["row_attrition"]` payload is less explicit than the
primary result-table contract.

Exact public contract affected:

- `phospy.api.DifferentialAnalysisResult.workflow_provenance`
- `phospy.api.DifferentialAnalysisResult.feature_eligibility`
- `phospy.provenance.RowAttritionRecord`
- `phospy.provenance.RowAttritionReport`

Proposed ownership:

PhosPy documentation/typing clarification.

Implementation scope:

- Document whether `workflow_provenance["row_attrition"]` is a supported payload
  schema for external consumers.
- If desired, add a narrow typed accessor on `DifferentialAnalysisResult` that
  returns a `RowAttritionReport | None`.
- Do not expose workflow provenance assemblers or validators.

Acceptance criteria:

- Docs state the supported inspection path.
- Tests verify the documented path without importing workflow internals.
- Existing provenance payloads remain backward compatible for current schema.

Non-goals:

- No assistant display DTO.
- No natural-language query API.
- No change to row filtering or scientific calculations.

## Overall conclusion

The current PhosPy public API is broadly ready for a separate assistant
application that consumes PhosPy as a scientific package rather than embedding
LLM concerns into PhosPy.

Already adequate:

- Constructing an analysis-ready dataset through the public builder.
- Discovering supported workflows through documented public exports.
- Running differential analysis through `DifferentialAnalysisWorkflow.run(...)`.
- Querying differential result rows by `site_key` through documented result
  tables.

Adequate with assistant-owned adapters:

- Workflow discovery presentation: the assistant can maintain a small explicit
  capability map over stable workflow classes.
- Validation handling for a first prototype: catch public exceptions and
  present messages without importing validators.
- Result lookup and presentation: filter documented public tables by `site_key`
  and format results in the assistant.
- Row-attrition presentation: normalize result-table status fields and
  provenance payloads without changing PhosPy.

Genuine or potential PhosPy public-contract gaps:

- Fully machine-readable validation diagnostic issue collections are not
  currently public.
- Public request objects do not have supported JSON-compatible round-trip
  serialization.
- A genuine rerunnable recipe contract does not exist.
- Differential row-attrition payload status is usable but could be documented or
  typed more explicitly.

Blocking assessment:

- No gap blocks a first assistant prototype that executes analyses in process,
  catches public exceptions, and keeps assistant-specific state/serialization in
  the assistant.
- A durable rerunnable recipe export is blocked by a real PhosPy contract gap if
  it is required for the first increment. It should be treated as a deferred
  PhosPy feature unless product scope requires it immediately.

Recommended next implementation step:

Build the first PhosPy Assistant prototype as a separate application using the
current public builder, request, workflow, and result contracts. Keep
assistant-owned adapters for workflow capability listing, validation-message
presentation, and `site_key` result lookup. Defer durable recipe export until
PhosPy defines a narrow versioned recipe/request serialization contract.

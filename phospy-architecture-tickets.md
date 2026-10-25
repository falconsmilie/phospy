# PhosPy Architecture Review — Ordered Codex Tickets

These seven tickets are listed in the recommended execution order. Each is self-contained and can be pasted into a fresh Codex session with access to the current repository. Review findings are starting evidence to reconcile with current code, not independently reverified source findings here.

ARCH-002 has been split into a confirmed reporting fix (ARCH-002A) and a separate ownership assessment (ARCH-002B). The latter is planning only; reader relocation is not presumed necessary.

## 1. ARCH-001 — Preserve complete kinase-scoring state in reloadable bundles

Priority: P1. Implement this ticket in the current PhosPy checkout.

**Context**

An architecture review at commit 2099e8d12e230ee7ffe914bfc5930e29e9b9df94 demonstrated that populated KinaseScoringResult.score_source_summary, motif_sequence_validation, and motif_library_validation become None after a kinase bundle round-trip. Other potentially omitted state includes score_source_matrix, profile_score_diagnostics, and scoring mode, source, scale, and scale metadata. Signalome reconstructs a nested scoring result separately. Existing equality helpers compare only selected scoring tables.

**Starting points**

- src/phospy/kinase/scoring/results.py
- src/phospy/io/bundles/_kinase/manifest.py and reconstruction.py
- src/phospy/io/bundles/_signalome/reconstruction.py
- tests/integration/test_kinase_bundle_integration.py
- docs/output_bundles.md
Locate equivalent current files if paths have moved.

**Working rules**

Read AGENTS.md and current architecture, API, and bundle-format policies. Record HEAD and initial dirty paths; preserve unrelated changes. The review commit is context, not a required baseline. This ticket needs no other ticket or conversation. If already resolved, verify and report that. Do not commit, push, or publish.

**Work**

1. Reproduce the loss with populated fields and add a failing regression test before changing persistence.
2. Inventory every public state-bearing scoring field. Define which state is persisted and which is faithfully derived on reconstruction. Document any justified transient exclusion; do not redefine missing durable state as optional to make the test pass.
3. Centralise the scoring-state encoding and decoding in a cohesive IO-owned component shared by kinase and Signalome bundles. Reuse existing serialization mechanisms and keep capability models independent of IO. Avoid a generic serialization framework or additional parallel field lists.
4. Persist the required tables, validation objects, diagnostics, and scoring semantics with their existing meaning, including None versus empty values. Reconstruct the original result type.
5. Version affected schemas using current policy. Preserve existing promised read support through existing mechanisms. Where no older-version support is promised, reject unsupported schemas clearly. Never invent historical state or silently replace missing semantics with authoritative defaults. Do not introduce compatibility shims.
6. Extend round-trip comparisons to the complete durable contract, including labelled table values and relevant dtypes, nested validation results, scalar semantics, and metadata. Include populated non-default fixtures so constructor defaults cannot conceal omissions.
7. Update affected format documentation and narrowly required generated artefacts.

**Acceptance and validation**

- The originally demonstrated fields survive a round-trip.
- Every supported materially distinct scoring mode is covered through kinase and nested Signalome persistence, using parametrisation rather than an unnecessary Cartesian product.
- Missing, malformed, or inconsistent newly required state is rejected through existing bundle error conventions.
- Integrity-before-reconstruction, transactional writes, public object identity, and defensive result access remain intact.
- Run focused codec/schema and bundle integration tests, affected architecture checks, and applicable Ruff/Pyright checks. Expand testing only for unresolved failures or repository requirements.
- Scientific computation and supported Python API exports remain unchanged.

**Completion report**

Give the durable-field inventory, schema versions and older-version behaviour, files changed, exact checks and results, and remaining limitations. Distinguish observed results from assumptions.

---

## 2. ARCH-004 — Repair the dependency-isolation guard and prove it detects current owners

Priority: P2. Implement this ticket in the current PhosPy checkout.

**Context**

An architecture review found that the runtime isolation test for phospy.quantitative and phospy.statistics checks retired phospy.science.activities and phospy.science.datasets prefixes. It can therefore pass while current capability packages are eagerly imported.

**Starting points**

- tests/architecture/test_package_dependency_dag.py
- Current owner classifications and logical ownership policy
- Existing subprocess import probes and removed-route tests

**Working rules**

Read AGENTS.md and current architecture policy. Record HEAD and initial dirty paths; preserve unrelated changes. Work against current code without requiring an old commit or earlier ticket. If already fixed, verify and report that. Do not commit, push, or publish.

**Work**

1. Establish the intended import-isolation rule from current policy. Derive applicable current capability roots from the maintained ownership authority and explicitly include any forbidden shared owners. Do not forbid every shared dependency indiscriminately.
2. Exercise quantitative and statistics independently in fresh Python processes so previously loaded modules cannot conceal import behaviour.
3. Check exact root names and their submodules in sys.modules. Do not rely on retired prefixes or arbitrary substring matching.
4. Make an unexpectedly empty forbidden-owner set a failure. Assert that a known current forbidden owner such as phospy.datasets is represented.
5. Add a positive control using the same guard and fresh-process mechanism: deliberately load a real forbidden current owner in a disposable probe and demonstrate that the guard detects it. Keep the intentional violation outside production files.
6. Preserve separate negative import-route tests for retired packages.
7. If the corrected guard exposes a real production import violation, trace and fix only the smallest import-boundary issue needed to satisfy the existing policy. Do not broaden permitted owners to obtain a pass.

**Acceptance and validation**

- Legitimate independent imports of quantitative and statistics pass.
- A deliberately loaded forbidden current owner is detected.
- Matching covers both the owner root and a descendant module.
- The guard remains tied to current ownership policy and cannot pass through an empty classification.
- Removed-route enforcement remains intact.
- Run the focused probes, their positive controls, and the architecture suite. Run applicable lint/type checks on changed files.
- Preserve public imports and runtime behaviour apart from correcting an actual prohibited eager dependency.

**Completion report**

Describe the isolation rule, resolved forbidden roots, positive-control evidence, any production violation discovered, and exact validation results.

---

## 3. ARCH-003 — Define deliberate internal validation contracts for run lineage

Priority: P2. Implement this ticket in the current PhosPy checkout.

**Context**

RunProvenance imports underscore-prefixed validation helpers from Preprocessing, Provenance, and References. Comments say these imports preserve validation order and error text. Those helpers currently function as undeclared cross-owner contracts.

**Starting points**

- src/phospy/run_lineage/models.py
- The preprocessing, provenance, and references helper definitions imported there
- Existing RunProvenance validation tests and architecture ownership rules

**Working rules**

Read AGENTS.md and current architecture/API policy. Record HEAD and dirty paths; preserve unrelated changes. This ticket is independent of other tickets and must adapt to the current tree. If already resolved, verify and report that. Do not commit, push, or publish.

**Work**

1. Trace every affected helper and all its callers. Identify the validation behaviour lineage actually needs, including sequencing, exception types, and currently contractual error messages.
2. Define the smallest cohesive owner-defined internal interface for that behaviour. Prefer an existing suitable contract or a focused extraction within the current owner. Keep ownership explicit and document the guarantees used by lineage.
3. Update lineage to consume those interfaces. A cosmetic rename alone is insufficient: make the intended cross-owner use, inputs, outputs, and behavioural obligations deliberate.
4. Keep these interfaces outside supported public facades and public API inventories. Do not create a generic shared-validation package, move unrelated responsibilities, or add pass-through layers solely to rename private functions.
5. Remove superseded helpers when unused, or retain them only for a demonstrated remaining internal responsibility. Avoid compatibility aliases.
6. Add focused regression coverage for the boundary and a proportionate architecture guard preventing recurrence of these private cross-owner dependencies. Do not impose a repository-wide private-import ban without assessing legitimate existing patterns.

**Acceptance and validation**

- Lineage no longer depends on private implementation helpers across the affected owners.
- Owner-defined internal interfaces have explicit responsibilities and no new dependency cycle.
- Valid inputs, invalid inputs, and multiple-invalid-field cases preserve required validation ordering, exception types, and contractual message behaviour.
- Existing owner tests and lineage tests pass.
- Public API inventories and class identities remain unchanged.
- Run focused validation tests, affected architecture checks, and applicable Ruff/Pyright checks.
- Do not change scientific algorithms or introduce new public API.

**Completion report**

List the previous dependencies, replacement contracts and their owners, preserved validation guarantees, deleted code, and exact validation results.

---

## 4. ARCH-005 — Require validated input at the Signalome interpreter boundary

Priority: P3. Implement this ticket in the current PhosPy checkout.

**Context**

SignalomeWorkflowInterpreter.run accepts both raw and validated requests. A private unchecked conversion in the validator constructs a validated DTO after only limited checks. The helper describes uses by focused tests and internal collaborators, so its full caller set must be verified before removal.

**Starting points**

- src/phospy/signalome/interpreter.py
- src/phospy/signalome/validator.py
- tests/unit/test_signalome_workflow_diagnostics.py
- Signalome workflow orchestration, validated request models, and related tests

**Working rules**

Read AGENTS.md and current architecture/API policy. Record HEAD and dirty paths; preserve unrelated changes. No other ticket or report is required. Reconcile the finding with current code and report verified prior resolution if applicable. Do not commit, push, or publish.

**Work**

1. Find all raw-request interpreter calls and unchecked-helper callers across production code, tests, tools, and relevant documentation. Include indirect uses where applicable; do not assume the seam exists solely for tests.
2. Restrict the interpreter contract to ValidatedSignalomeWorkflowRequest. A raw request reaching it must fail clearly without an unchecked conversion.
3. Route legitimate raw production requests through the existing validator before interpretation, preserving the public workflow entry point and validator → interpreter → executor sequence. Avoid duplicate validation.
4. Update focused tests to use the real validator or a test-owned fixture factory that constructs explicitly controlled validated input. Keep malformed internal fixtures deliberate and local to the tests that require them.
5. Remove the unchecked production helper and its obsolete references. Do not replace it with another production bypass or compatibility alias.

**Acceptance and validation**

- The public Signalome workflow still accepts its supported raw request and validates it.
- Interpreter calls consume validated requests exclusively.
- Raw interpreter input is rejected; no partial validation path fabricates a validated request.
- Previously legitimate internal callers have been migrated.
- Focused tests still cover interpreter diagnostics and error cases.
- Run relevant validator, interpreter, workflow, and architecture tests plus applicable Ruff/Pyright checks.
- Supported public APIs, scientific outputs, dataset-view ownership, and allocation safeguards remain intact.

**Completion report**

Provide the caller inventory and disposition, the enforced boundary, removed helper, exact test results, and any remaining uncertainty.

---

## 5. ARCH-002A — Report complete and policy-filtered dependency graphs accurately

Priority: P2. Implement the confirmed reporting/enforcement correction in the current PhosPy checkout.

**Context**

The owner graph contains datasets → io through the builder's concrete path reader, and io → datasets through reconstruction/model dependencies. A permanent classification excludes the former edge from the graph checked for acyclicity. Documentation reportedly calls the exact observed graph acyclic without that qualification.

The narrow builder-to-reader import was explicitly authorised during the refactor. An accepted exception must not be treated automatically as an unauthorised dependency.

**Starting points**

- src/phospy/datasets/builder.py
- architecture/permanent-condensation-classifications.json
- tests/architecture/test_package_dependency_dag.py
- tests/architecture/logical_ownership_policy.py and import_graph_support.py
- docs/architecture.md and current graph-reporting tools

**Working rules**

Read AGENTS.md, current ownership policy, relevant ADRs, and classification records. Record HEAD and dirty paths; preserve unrelated changes. No earlier ticket is required. Use current evidence rather than requiring the historical review baseline. Do not commit, push, or publish.

**Work**

1. Confirm how observed edges are collected and exactly where classifications are applied.
2. Represent and report two clearly named views using the existing analyser:
   - Complete statically observed owner graph before classifications, including strongly connected components.
   - Policy-filtered owner graph used by the existing DAG gate after accepted classifications.
3. Report excluded edges with their classification identifiers and rationale. Keep output deterministic and avoid building a second competing graph framework.
4. Correct active documentation and diagnostics so claims identify the graph they describe. State relevant static-analysis limitations.
5. Preserve the existing policy gate and exception precision. Do not add exclusions, silently broaden classifications, or claim that an accepted raw cycle means the DAG enforcement itself is broken.
6. Add regression coverage showing that classified edges remain visible in the complete view and are removed only from the policy-filtered view. Demonstrate that an unclassified cycle still fails enforcement.

**Acceptance and validation**

- The actual raw datasets/io component is reported if it still exists.
- A synthetic classified-cycle case proves both graph views even if current production dependencies have changed.
- Existing accepted classifications remain auditable and narrowly enforced.
- Unclassified cycles remain failures.
- Documentation makes no unqualified claim of raw-graph acyclicity when raw cycles exist.
- Run focused graph/ownership tests, the architecture suite, and applicable lint/type checks.

**Scope boundary and completion**

Do not move the reader, redesign the builder, or change the accepted ownership policy in this ticket. Report current edges/components, classification treatment, documentation changes, and exact checks. Reader relocation requires a separate ownership assessment.

---

## 6. ARCH-002B — Assess ownership of dataset path ingestion without presupposing a move

Priority: P3. Complete a bounded architecture decision assessment in the current PhosPy checkout. This is a planning ticket, not authorisation for a reader migration.

**Context**

AnalysisReadyDatasetBuilder reportedly imports and constructs phospy.io.readers.dataset_inputs.DatasetPathTableReader, while IO also depends on dataset models. The narrow datasets.builder → reader import was explicitly authorised during the refactor and later classified. A review proposed moving this reader into Datasets to remove the raw owner cycle.

An existing accepted dependency is a tradeoff to evaluate. Eliminating a cycle alone does not establish that moving IO responsibility is the best design.

**Starting points**

- src/phospy/datasets/builder.py
- src/phospy/io/readers/dataset_inputs.py
- Dataset-owned reader protocol and input-resolution code
- ADR-0053, ADR-0054, builder migration decisions, and current classification records
- Dataset importer and builder behaviour tests

**Working rules**

Read AGENTS.md and the current architecture/API policy. Record HEAD and dirty paths; preserve unrelated changes. This ticket needs no earlier graph-reporting change: inspect the current imports and tools directly. Do not commit, push, publish, move production modules, or change policy acceptance status.

**Work**

1. Establish the actual current dependency and the most recent applicable decision. Distinguish historical authorisation from current policy.
2. Inspect the reader's real responsibilities and consumers. Identify whether it primarily owns dataset semantics, general file-format IO, or a mixture.
3. Compare retaining the narrowly classified dependency with moving the reader into Datasets. Assess another option only if current code provides a concrete reason; avoid speculative redesign.
4. Evaluate responsibility cohesion, change coupling, testing, import cost, public object identity, and maintenance cost. A classified owner cycle alone is not evidence of a Python import failure.
5. Check behaviour that a future move must preserve: string and PathLike paths, blank paths, allow_none, read/error ordering and messages, reader exception propagation, injected-reader behaviour, and fallback handling for false-valued injected objects. Use current contracts to establish expected behaviour.
6. Record a concise decision assessment in the existing architecture-documentation location. Recommend retaining or changing ownership with evidence, uncertainties, and concrete consequences. Do not mark a new ADR accepted or rewrite historical acceptance records.

**Acceptance and validation**

- The assessment is based on current code and current policy.
- Both principal options receive a concrete comparison.
- The recommendation identifies a practical benefit or explains why retention is sufficient.
- If a move is recommended, provide an independently executable follow-up scope: exact modules/consumers, classification changes, behaviour checks, and public-contract constraints. Do not implement it.
- If the dependency is already resolved, document that status and avoid inventing additional restructuring.
- Run focused existing checks only when needed to resolve an uncertainty.

**Completion report**

Provide the decision-note path, current dependency status, recommendation, supporting evidence, checks performed, and any future implementation scope. Completing this assessment does not depend on obtaining approval for a migration.

---

## 7. ARCH-006 — Refresh active architecture guidance and label generated reports accurately

Priority: P3. Implement this documentation/tooling cleanup in the current PhosPy checkout.

**Context**

The active maintenance guide reportedly assigns configuration DTOs to retired phospy.contracts.configs. The validation ownership map also describes phospy.contracts as a current owner. The architecture inventory generator emits malformed introductory prose, and the orchestration candidate report contains retired phospy.science.* routes.

Historical ledgers, superseded ADRs, negative import tests, and versioned payload evidence legitimately contain old names and must be preserved.

**Starting points**

- docs/maintenance.md
- Locate validation-ownership.md in the current docs tree
- docs/architecture_inventory.md
- tools/architecture/generate_architecture_inventory.py
- docs/testing/orchestration_test_candidates.md and its generator

**Working rules**

Read AGENTS.md and current architecture/API policy. Record HEAD and dirty paths; preserve unrelated changes. This ticket is standalone and must describe the repository as it currently exists, whether other cleanup tickets have run or not. Do not commit, push, or publish.

**Work**

1. Derive current ownership from source and accepted policy. Correct active placement guidance for requests, configurations, validation, passive records, and shared contracts where the identified stale material affects them.
2. Inspect adjacent active guidance for the same obsolete ownership claims. Keep this a focused architecture cleanup, not a documentation redesign.
3. Fix malformed inventory prose in the generator, then regenerate only its intended outputs. Preserve its existing freshness/check mechanism.
4. Decide whether the orchestration candidate report is an active maintained aid or a historical snapshot:
   - If active, regenerate in an isolated disposable checkout or equivalent safe workspace, inspect all outputs, and apply only justified changes.
   - If historical, label its status and scope clearly. Include its generation revision/date only when verifiable; do not invent provenance.
5. Inspect generator side effects before running it. Do not accept unrelated test-inventory annotation churn automatically.
6. Preserve useful historical references and current compatibility/serialization contracts. Avoid global search-and-replace of retired package names.
7. Do not add a permanent CI freshness gate for a heuristic historical report without an established maintenance need.

**Acceptance and validation**

- Current placement guidance matches current owners.
- Obsolete names in active material are clearly historical or explanatory rather than instructions to use removed routes.
- Inventory regeneration produces coherent prose and passes the existing freshness check.
- The orchestration report has an explicit current or historical status.
- Historical records and negative route tests remain intact.
- Run the relevant generator check and existing documentation build/link checks. Update existing generator tests if affected; do not add tests that merely repeat prose.

**Completion report**

List corrected guidance, the generated-report disposition, intentional historical references, generator side effects handled, and exact validation results.


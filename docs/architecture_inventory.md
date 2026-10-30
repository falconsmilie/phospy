# Architecture inventory

This generated inventory records tracked production files, module ownership,
and the complete and policy-filtered dependency views. Executable architecture
tests compare it with tracked source, the one-row-per-module owner map, the
exact owner-edge policy, and the permanent classification registry.

## Counts and architecture assertions

- Tracked production files: **646**.
- Tracked Python modules: **640**.
- Owner-map entries: **640**.
- Temporary architecture exceptions: **0**.
- Retired tracked package roots (`contracts`, `errors`, `frames`, `policies`, `science`, `tables`, `validation`, `workflows`): **0**.
- Unclassified structural forwarding candidates, old-path literals, compatibility behavior, and lazy loading: **0**.
- The complete eager module graph and policy-filtered owner graph are required to be acyclic by maintained tests.

## Top-level roots

| Root | Current classification |
| --- | --- |
| `.` | `public_facade_machinery` |
| `api` | `public_facade` |
| `core` | `organisational_namespace` |
| `data` | `packaged_resource_payload` |
| `datasets` | `capability` |
| `differential` | `capability` |
| `enrichment` | `capability` |
| `io` | `adapter` |
| `kinase` | `capability` |
| `preparation_records` | `passive_shared_record_foundation` |
| `preprocessing` | `capability` |
| `run_lineage` | `composition_foundation` |
| `signalome` | `capability` |

## Complete statically observed owner graph (before classifications)

This view contains **129** owner edges before
accepted condensation classifications are applied. Every edge remains
subject to the maintained exact-rule comparison; a classified edge is an
auditable accepted dependency, not an unauthorised import.

### Strongly connected components

- `datasets, io, kinase, signalome`
- Single-owner components: **19**.

### Complete observed edges

| Source | Target | Classification | Authority |
| --- | --- | --- | --- |
| `activities` | `errors` | `shared` | ACT-10 exact activities row |
| `activities` | `frames` | `shared` | ACT-10 exact activities row |
| `activities` | `quantitative` | `shared` | ACT-10 exact activities row |
| `datasets` | `errors` | `capability` | ADR-0053 datasets direction |
| `datasets` | `evidence` | `capability` | ADR-0053 datasets direction |
| `datasets` | `frames` | `capability` | ADR-0053 datasets direction |
| `datasets` | `io` | `permanent_classification` | ADR-0053 D05B-1 dataset builder path-reader authority |
| `datasets` | `phosphosites` | `capability` | ADR-0053 datasets direction |
| `datasets` | `preparation_records` | `capability` | ADR-0053 datasets direction |
| `datasets` | `preprocessing` | `capability` | ADR-0053 datasets direction |
| `datasets` | `provenance` | `capability` | ADR-0053 datasets direction |
| `datasets` | `quantitative` | `capability` | ADR-0053 datasets direction |
| `datasets` | `references` | `capability` | ADR-0053 datasets direction |
| `datasets` | `run_lineage` | `capability` | ADR-0053 datasets direction |
| `datasets` | `taxonomy` | `capability` | ADR-0053 datasets direction |
| `datasets` | `validation` | `permanent_classification` | ADR-0053 D05B-2 dataset input-source validation authority |
| `differential` | `datasets` | `capability` | ADR-0053 differential direction |
| `differential` | `errors` | `capability` | ADR-0053 differential direction |
| `differential` | `evidence` | `capability` | ADR-0053 differential direction |
| `differential` | `frames` | `capability` | ADR-0053 differential direction |
| `differential` | `phosphosites` | `capability` | ADR-0053 differential direction |
| `differential` | `policies` | `capability` | ADR-0053 differential direction |
| `differential` | `preparation_records` | `capability` | ADR-0053 differential direction |
| `differential` | `provenance` | `capability` | ADR-0053 differential direction |
| `differential` | `quantitative` | `capability` | ADR-0053 differential direction |
| `differential` | `references` | `capability` | ADR-0053 differential direction |
| `differential` | `run_lineage` | `capability` | ADR-0053 differential direction |
| `differential` | `statistics` | `capability` | ADR-0053 differential direction |
| `enrichment` | `errors` | `capability` | ADR-0053 enrichment direction |
| `enrichment` | `frames` | `capability` | ADR-0053 enrichment direction |
| `enrichment` | `provenance` | `capability` | ADR-0053 enrichment direction |
| `enrichment` | `quantitative` | `capability` | ADR-0053 enrichment direction |
| `enrichment` | `run_lineage` | `capability` | ADR-0053 enrichment direction |
| `enrichment` | `statistics` | `capability` | ADR-0053 enrichment direction |
| `evidence` | `errors` | `shared` | ADR-0053 evidence row |
| `evidence` | `frames` | `shared` | ADR-0053 evidence row |
| `evidence` | `phosphosites` | `shared` | ADR-0053 evidence row |
| `evidence` | `policies` | `shared` | ADR-0053 evidence row |
| `evidence` | `quantitative` | `shared` | ADR-0053 evidence row |
| `frames` | `errors` | `shared` | ADR-0053 frames row |
| `io` | `activities` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `datasets` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `enrichment` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `errors` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `evidence` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `frames` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `kinase` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `phosphosites` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `preparation_records` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `preprocessing` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `provenance` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `quantitative` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `references` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `run_lineage` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `signalome` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `taxonomy` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `validation` | `adapter` | ADR-0053 inward IO adapter |
| `kinase` | `activities` | `capability` | ADR-0053 kinase direction |
| `kinase` | `datasets` | `capability` | ADR-0053 kinase direction |
| `kinase` | `errors` | `capability` | ADR-0053 kinase direction |
| `kinase` | `evidence` | `capability` | ADR-0053 kinase direction |
| `kinase` | `frames` | `capability` | ADR-0053 kinase direction |
| `kinase` | `phosphosites` | `capability` | ADR-0053 kinase direction |
| `kinase` | `policies` | `capability` | ADR-0053 kinase direction |
| `kinase` | `provenance` | `capability` | ADR-0053 kinase direction |
| `kinase` | `quantitative` | `capability` | ADR-0053 kinase direction |
| `kinase` | `references` | `capability` | ADR-0053 kinase direction |
| `kinase` | `run_lineage` | `capability` | ADR-0053 kinase direction |
| `kinase` | `statistics` | `capability` | ADR-0053 kinase direction |
| `kinase` | `taxonomy` | `capability` | ADR-0053 kinase direction |
| `phosphosites` | `errors` | `shared` | ADR-0053 phosphosites row |
| `phosphosites` | `frames` | `shared` | ADR-0053 phosphosites row |
| `phosphosites` | `taxonomy` | `shared` | ADR-0053 phosphosites row |
| `policies` | `errors` | `shared` | ADR-0053 policies row |
| `preparation_records` | `errors` | `shared` | ADR-0053 preparation_records row |
| `preparation_records` | `frames` | `shared` | ADR-0053 preparation_records row |
| `preparation_records` | `policies` | `shared` | ADR-0053 preparation_records row |
| `preparation_records` | `provenance` | `shared` | ADR-0053 preparation_records row |
| `preprocessing` | `errors` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `evidence` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `frames` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `phosphosites` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `policies` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `preparation_records` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `provenance` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `quantitative` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `references` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `taxonomy` | `capability` | ADR-0053 preprocessing direction |
| `provenance` | `errors` | `shared` | ADR-0053 provenance row |
| `public_facades` | `datasets` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `differential` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `enrichment` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `errors` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `evidence` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `io` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `kinase` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `preprocessing` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `provenance` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `quantitative` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `references` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `signalome` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `statistics` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `taxonomy` | `public_facade` | ADR-0054 aggregate facade |
| `quantitative` | `errors` | `shared` | ADR-0053 quantitative row |
| `quantitative` | `provenance` | `shared` | ADR-0053 quantitative row |
| `references` | `errors` | `shared` | ADR-0053 references row |
| `references` | `frames` | `shared` | ADR-0053 references row |
| `references` | `phosphosites` | `shared` | ADR-0053 references row |
| `references` | `policies` | `shared` | ADR-0053 references row |
| `references` | `provenance` | `shared` | ADR-0053 references row |
| `references` | `taxonomy` | `shared` | ADR-0053 references row |
| `run_lineage` | `errors` | `composition` | ADR-0053 run-lineage composition row |
| `run_lineage` | `evidence` | `composition` | ADR-0053 run-lineage composition row |
| `run_lineage` | `preprocessing` | `composition` | ADR-0053 run-lineage composition row |
| `run_lineage` | `provenance` | `composition` | ADR-0053 run-lineage composition row |
| `run_lineage` | `references` | `composition` | ADR-0053 run-lineage composition row |
| `signalome` | `datasets` | `capability` | ADR-0053 signalome direction |
| `signalome` | `errors` | `capability` | ADR-0053 signalome direction |
| `signalome` | `evidence` | `capability` | ADR-0053 signalome direction |
| `signalome` | `frames` | `capability` | ADR-0053 signalome direction |
| `signalome` | `kinase` | `capability` | ADR-0053 signalome direction |
| `signalome` | `phosphosites` | `capability` | ADR-0053 signalome direction |
| `signalome` | `policies` | `capability` | ADR-0053 signalome direction |
| `signalome` | `provenance` | `capability` | ADR-0053 signalome direction |
| `signalome` | `quantitative` | `capability` | ADR-0053 signalome direction |
| `signalome` | `references` | `capability` | ADR-0053 signalome direction |
| `signalome` | `run_lineage` | `capability` | ADR-0053 signalome direction |
| `signalome` | `validation` | `capability` | ADR-0053 signalome direction |
| `validation` | `errors` | `shared` | ADR-0053 validation row |

## Excluded imports from accepted classifications

These exact module imports remain visible in the complete view and are
removed only before policy condensation and the DAG gate.

| Classification ID | Source path | Source owner | Target owner | Exact import | Rationale |
| --- | --- | --- | --- | --- | --- |
| `dataset-builder-path-reader` | `src/phospy/datasets/builder.py` | `datasets` | `io` | `phospy.datasets.builder -> phospy.io.readers.dataset_inputs.DatasetPathTableReader` | ADR-0053 D05B-1 authorises this exact Dataset builder acquisition of the IO-owned path reader while IO retains filesystem ownership. |
| `dataset-input-source-validation` | `src/phospy/datasets/construction/input_resolution.py` | `datasets` | `validation` | `phospy.datasets.construction.input_resolution -> phospy.core.validation.common.input_sources.validate_dataframe_or_path` | ADR-0053 D05B-2 authorises this exact Dataset input-resolution use of the shared DataFrame-or-path validation helper. |
| `dataset-request-source-validation` | `src/phospy/datasets/validation/inputs.py` | `datasets` | `validation` | `phospy.datasets.validation.inputs -> phospy.core.validation.common.input_sources.validate_dataframe_or_path` | ADR-0053 D05B-2 authorises this exact Dataset request-validation use of the shared DataFrame-or-path validation helper. |

## Policy-filtered owner graph used by the DAG gate

This view contains **127** owner edges
after the exact accepted imports above are removed. It is the graph checked
for owner-level acyclicity.

### Strongly connected components

- No multi-owner strongly connected components; the view is acyclic.
- Single-owner components: **23**.

### Policy-filtered edges

| Source | Target | Classification | Authority |
| --- | --- | --- | --- |
| `activities` | `errors` | `shared` | ACT-10 exact activities row |
| `activities` | `frames` | `shared` | ACT-10 exact activities row |
| `activities` | `quantitative` | `shared` | ACT-10 exact activities row |
| `datasets` | `errors` | `capability` | ADR-0053 datasets direction |
| `datasets` | `evidence` | `capability` | ADR-0053 datasets direction |
| `datasets` | `frames` | `capability` | ADR-0053 datasets direction |
| `datasets` | `phosphosites` | `capability` | ADR-0053 datasets direction |
| `datasets` | `preparation_records` | `capability` | ADR-0053 datasets direction |
| `datasets` | `preprocessing` | `capability` | ADR-0053 datasets direction |
| `datasets` | `provenance` | `capability` | ADR-0053 datasets direction |
| `datasets` | `quantitative` | `capability` | ADR-0053 datasets direction |
| `datasets` | `references` | `capability` | ADR-0053 datasets direction |
| `datasets` | `run_lineage` | `capability` | ADR-0053 datasets direction |
| `datasets` | `taxonomy` | `capability` | ADR-0053 datasets direction |
| `differential` | `datasets` | `capability` | ADR-0053 differential direction |
| `differential` | `errors` | `capability` | ADR-0053 differential direction |
| `differential` | `evidence` | `capability` | ADR-0053 differential direction |
| `differential` | `frames` | `capability` | ADR-0053 differential direction |
| `differential` | `phosphosites` | `capability` | ADR-0053 differential direction |
| `differential` | `policies` | `capability` | ADR-0053 differential direction |
| `differential` | `preparation_records` | `capability` | ADR-0053 differential direction |
| `differential` | `provenance` | `capability` | ADR-0053 differential direction |
| `differential` | `quantitative` | `capability` | ADR-0053 differential direction |
| `differential` | `references` | `capability` | ADR-0053 differential direction |
| `differential` | `run_lineage` | `capability` | ADR-0053 differential direction |
| `differential` | `statistics` | `capability` | ADR-0053 differential direction |
| `enrichment` | `errors` | `capability` | ADR-0053 enrichment direction |
| `enrichment` | `frames` | `capability` | ADR-0053 enrichment direction |
| `enrichment` | `provenance` | `capability` | ADR-0053 enrichment direction |
| `enrichment` | `quantitative` | `capability` | ADR-0053 enrichment direction |
| `enrichment` | `run_lineage` | `capability` | ADR-0053 enrichment direction |
| `enrichment` | `statistics` | `capability` | ADR-0053 enrichment direction |
| `evidence` | `errors` | `shared` | ADR-0053 evidence row |
| `evidence` | `frames` | `shared` | ADR-0053 evidence row |
| `evidence` | `phosphosites` | `shared` | ADR-0053 evidence row |
| `evidence` | `policies` | `shared` | ADR-0053 evidence row |
| `evidence` | `quantitative` | `shared` | ADR-0053 evidence row |
| `frames` | `errors` | `shared` | ADR-0053 frames row |
| `io` | `activities` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `datasets` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `enrichment` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `errors` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `evidence` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `frames` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `kinase` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `phosphosites` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `preparation_records` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `preprocessing` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `provenance` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `quantitative` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `references` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `run_lineage` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `signalome` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `taxonomy` | `adapter` | ADR-0053 inward IO adapter |
| `io` | `validation` | `adapter` | ADR-0053 inward IO adapter |
| `kinase` | `activities` | `capability` | ADR-0053 kinase direction |
| `kinase` | `datasets` | `capability` | ADR-0053 kinase direction |
| `kinase` | `errors` | `capability` | ADR-0053 kinase direction |
| `kinase` | `evidence` | `capability` | ADR-0053 kinase direction |
| `kinase` | `frames` | `capability` | ADR-0053 kinase direction |
| `kinase` | `phosphosites` | `capability` | ADR-0053 kinase direction |
| `kinase` | `policies` | `capability` | ADR-0053 kinase direction |
| `kinase` | `provenance` | `capability` | ADR-0053 kinase direction |
| `kinase` | `quantitative` | `capability` | ADR-0053 kinase direction |
| `kinase` | `references` | `capability` | ADR-0053 kinase direction |
| `kinase` | `run_lineage` | `capability` | ADR-0053 kinase direction |
| `kinase` | `statistics` | `capability` | ADR-0053 kinase direction |
| `kinase` | `taxonomy` | `capability` | ADR-0053 kinase direction |
| `phosphosites` | `errors` | `shared` | ADR-0053 phosphosites row |
| `phosphosites` | `frames` | `shared` | ADR-0053 phosphosites row |
| `phosphosites` | `taxonomy` | `shared` | ADR-0053 phosphosites row |
| `policies` | `errors` | `shared` | ADR-0053 policies row |
| `preparation_records` | `errors` | `shared` | ADR-0053 preparation_records row |
| `preparation_records` | `frames` | `shared` | ADR-0053 preparation_records row |
| `preparation_records` | `policies` | `shared` | ADR-0053 preparation_records row |
| `preparation_records` | `provenance` | `shared` | ADR-0053 preparation_records row |
| `preprocessing` | `errors` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `evidence` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `frames` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `phosphosites` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `policies` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `preparation_records` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `provenance` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `quantitative` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `references` | `capability` | ADR-0053 preprocessing direction |
| `preprocessing` | `taxonomy` | `capability` | ADR-0053 preprocessing direction |
| `provenance` | `errors` | `shared` | ADR-0053 provenance row |
| `public_facades` | `datasets` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `differential` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `enrichment` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `errors` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `evidence` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `io` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `kinase` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `preprocessing` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `provenance` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `quantitative` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `references` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `signalome` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `statistics` | `public_facade` | ADR-0054 aggregate facade |
| `public_facades` | `taxonomy` | `public_facade` | ADR-0054 aggregate facade |
| `quantitative` | `errors` | `shared` | ADR-0053 quantitative row |
| `quantitative` | `provenance` | `shared` | ADR-0053 quantitative row |
| `references` | `errors` | `shared` | ADR-0053 references row |
| `references` | `frames` | `shared` | ADR-0053 references row |
| `references` | `phosphosites` | `shared` | ADR-0053 references row |
| `references` | `policies` | `shared` | ADR-0053 references row |
| `references` | `provenance` | `shared` | ADR-0053 references row |
| `references` | `taxonomy` | `shared` | ADR-0053 references row |
| `run_lineage` | `errors` | `composition` | ADR-0053 run-lineage composition row |
| `run_lineage` | `evidence` | `composition` | ADR-0053 run-lineage composition row |
| `run_lineage` | `preprocessing` | `composition` | ADR-0053 run-lineage composition row |
| `run_lineage` | `provenance` | `composition` | ADR-0053 run-lineage composition row |
| `run_lineage` | `references` | `composition` | ADR-0053 run-lineage composition row |
| `signalome` | `datasets` | `capability` | ADR-0053 signalome direction |
| `signalome` | `errors` | `capability` | ADR-0053 signalome direction |
| `signalome` | `evidence` | `capability` | ADR-0053 signalome direction |
| `signalome` | `frames` | `capability` | ADR-0053 signalome direction |
| `signalome` | `kinase` | `capability` | ADR-0053 signalome direction |
| `signalome` | `phosphosites` | `capability` | ADR-0053 signalome direction |
| `signalome` | `policies` | `capability` | ADR-0053 signalome direction |
| `signalome` | `provenance` | `capability` | ADR-0053 signalome direction |
| `signalome` | `quantitative` | `capability` | ADR-0053 signalome direction |
| `signalome` | `references` | `capability` | ADR-0053 signalome direction |
| `signalome` | `run_lineage` | `capability` | ADR-0053 signalome direction |
| `signalome` | `validation` | `capability` | ADR-0053 signalome direction |
| `validation` | `errors` | `shared` | ADR-0053 validation row |

## Static-analysis scope

The analyser includes absolute and relative imports, imports nested under
`TYPE_CHECKING` or function bodies, and literal `importlib.import_module` or
`__import__` calls. It does not resolve computed module names, bindings created
through assignment or re-export, plugin discovery, import hooks, or other
runtime-generated imports. The complete owner view therefore reports static
dependency potential; the separately enforced eager-module graph describes
imports reachable while modules and class bodies execute.

## Public facades, adapters, and foundations

The root table is authoritative for every top-level owner. The public
supported facades are `phospy.api` and `phospy.api.advanced`; the
`phospy` root provides metadata only, and `phospy.io`
is the inward adapter. Shared scientific, record, composition, and
technical foundations are the roots carrying those classifications above.

Exact permanent facade/safety modules:

- `src/phospy/api/_inventory.py`: authoritative aggregate inventory data; explicitly free of compatibility routes
- `src/phospy/differential/aggregation/experimental.py`: permanent fail-closed scientific safety seam governed by ADR-0041; it is not a supported compatibility import route

Exact permanent lazy loaders:

- `src/phospy/io/__init__.py`: IO adapter composition
- `src/phospy/io/bundles/__init__.py`: IO bundle adapter composition
- `src/phospy/io/readers/__init__.py`: IO reader adapter composition
- `src/phospy/preprocessing/contracts/plan.py`: cycle-free owner-internal plan/interpreter composition

## Structurally detected forwarding candidates

The AST detector identifies modules whose executable top-level body is
limited to imports, `__all__` maintenance, and direct aliases of imported
names. A docstring and `__future__` imports do not affect detection. The
detector establishes forwarding-shaped structure only; each candidate
still requires the explicit policy disposition below.

The declared scope does not interpret computed exports, forwarding inside
functions or classes, module `__getattr__`, import hooks, or arbitrary
runtime binding, and it does not attempt to model arbitrary Python behavior.

| Source path | Policy disposition | Rationale |
| --- | --- | --- |
| `src/phospy/api/__init__.py` | `supported_aggregate_facade` | ADR-0054 ordinary aggregate API; exact owner identities only. |
| `src/phospy/api/advanced/__init__.py` | `supported_aggregate_facade` | ADR-0055 specialist aggregate API; exact owner identities only. |
| `src/phospy/api/builders.py` | `deliberate_internal_facade` | Internal aggregate composition used by the supported API; unsupported as an external route. |
| `src/phospy/core/errors/__init__.py` | `deliberate_internal_facade` | Owner-package aggregation for internal imports. |
| `src/phospy/core/evidence/__init__.py` | `deliberate_internal_facade` | Owner-package aggregation for evidence records. |
| `src/phospy/core/frames/__init__.py` | `deliberate_internal_facade` | Owner-package aggregation for neutral frame primitives. |
| `src/phospy/core/phosphosites/__init__.py` | `deliberate_internal_facade` | Phosphosite owner-package composition. |
| `src/phospy/core/policies/__init__.py` | `deliberate_internal_facade` | Shared policy owner-package composition. |
| `src/phospy/core/provenance/__init__.py` | `deliberate_internal_facade` | Provenance owner-package composition. |
| `src/phospy/core/provenance/models/__init__.py` | `deliberate_internal_facade` | Provenance model package composition. |
| `src/phospy/core/provenance/serialization/__init__.py` | `deliberate_internal_facade` | Provenance serialization package composition. |
| `src/phospy/core/quantitative/__init__.py` | `deliberate_internal_facade` | Quantitative owner-package composition. |
| `src/phospy/core/quantitative/models.py` | `deliberate_internal_facade` | Quantitative model aggregation retained for owner-internal composition. |
| `src/phospy/core/quantitative/transformers/__init__.py` | `deliberate_internal_facade` | Quantitative transformer package composition. |
| `src/phospy/core/references/manifest.py` | `deliberate_internal_facade` | Reference manifest contract composition. |
| `src/phospy/core/references/sequences/__init__.py` | `deliberate_internal_facade` | Reference sequence package composition. |
| `src/phospy/core/references/validation/__init__.py` | `deliberate_internal_facade` | Reference validation package composition. |
| `src/phospy/core/statistics/__init__.py` | `deliberate_internal_facade` | Statistics owner-package composition. |
| `src/phospy/core/taxonomy/__init__.py` | `deliberate_internal_facade` | Taxonomy owner-package composition. |
| `src/phospy/io/bundles/kinase.py` | `deliberate_internal_facade` | IO-owned Kinase bundle composition. |
| `src/phospy/io/bundles/signalome.py` | `deliberate_internal_facade` | IO-owned Signalome bundle composition. |
| `src/phospy/io/publishers/__init__.py` | `deliberate_internal_facade` | IO publisher package composition. |
| `src/phospy/io/readers/fragpipe/__init__.py` | `deliberate_internal_facade` | FragPipe reader package composition. |
| `src/phospy/kinase/scoring/motif/__init__.py` | `deliberate_internal_facade` | Kinase motif implementation package composition. |
| `src/phospy/preprocessing/batch_correction/__init__.py` | `deliberate_internal_facade` | Batch-correction capability package composition. |
| `src/phospy/preprocessing/batch_correction/adapters/__init__.py` | `deliberate_internal_facade` | Batch-correction adapter package composition. |
| `src/phospy/preprocessing/batch_correction/contracts/__init__.py` | `deliberate_internal_facade` | Batch-correction contract package composition. |
| `src/phospy/preprocessing/batch_correction/execution/__init__.py` | `deliberate_internal_facade` | Batch-correction execution package composition. |
| `src/phospy/preprocessing/batch_correction/provenance/__init__.py` | `deliberate_internal_facade` | Batch-correction provenance package composition. |
| `src/phospy/preprocessing/batch_correction/validation/__init__.py` | `deliberate_internal_facade` | Batch-correction validation package composition. |
| `src/phospy/preprocessing/contracts/__init__.py` | `deliberate_internal_facade` | Preprocessing contract package composition. |
| `src/phospy/preprocessing/contracts/configuration/__init__.py` | `deliberate_internal_facade` | Preprocessing configuration package composition. |
| `src/phospy/preprocessing/sequence/__init__.py` | `deliberate_internal_facade` | Preprocessing sequence package composition. |
| `src/phospy/preprocessing/stages/__init__.py` | `deliberate_internal_facade` | Preprocessing stage package composition. |
| `src/phospy/preprocessing/stages/missing_data/__init__.py` | `deliberate_internal_facade` | Missing-data stage package composition. |
| `src/phospy/preprocessing/stages/site_matrix_components/__init__.py` | `deliberate_internal_facade` | Site-matrix stage component composition. |
| `src/phospy/preprocessing/workflow.py` | `deliberate_internal_facade` | Preprocessing-owned workflow composition route. |
| `src/phospy/run_lineage/__init__.py` | `deliberate_internal_facade` | Run-lineage composition foundation facade. |
| `src/phospy/signalome/clustering/backends/__init__.py` | `deliberate_internal_facade` | Signalome clustering backend package composition. |
| `src/phospy/signalome/tables/__init__.py` | `deliberate_internal_facade` | Signalome table package composition. |

## Classified versioned and historical behavior

These records distinguish legitimate versioned schema or legacy-data
identifiers from live compatibility routes. Their exact sets are enforced.

### Retired-route string literals

- `src/phospy/core/quantitative/authority.py`
  - `phospy.science.batch_correction.sps_discovery`
  - `phospy.science.datasets.builders.transformation_resolver`
  - `phospy.science.datasets.preprocessing.state_builder`
  - `phospy.science.transformations.transformers.identity`
- `src/phospy/core/quantitative/policy.py`
  - `phospy.science.batch_correction.sps_discovery`
  - `phospy.science.datasets.preprocessing.state_builder`
  - `phospy.science.transformations.transformers.identity`
- `src/phospy/core/quantitative/quantitative_contracts.py`
  - `phospy.science.transformations.quantitative_contracts`
- `src/phospy/core/quantitative/scale_values.py`
  - `phospy.science.transformations.transformers.log2`
- `src/phospy/core/quantitative/transformers/identity.py`
  - `phospy.science.transformations.transformers.identity.IdentityTransformer`
- `src/phospy/core/quantitative/transformers/log2.py`
  - `phospy.science.transformations.transformers.log2`
  - `phospy.science.transformations.transformers.log2.Log2Transformer`
- `src/phospy/datasets/construction/provenance.py`
  - `phospy.science.datasets.builders.transformation_state.DatasetTransformationStateResolver.run`
- `src/phospy/datasets/construction/state_assembly.py`
  - `phospy.science.datasets.preprocessing.state_builder`
- `src/phospy/datasets/construction/transformation_state.py`
  - `phospy.science.datasets.builders.executor.input_intensity_scale`
- `src/phospy/kinase/references/projection_contracts.py`
  - `phospy.workflows.kinase.reference_projector.v1`
- `src/phospy/kinase/site_sequence_support.py`
  - `phospy.workflows.kinase.site_sequence_support_builder.v1`
- `src/phospy/preprocessing/control_sites/discovery_contracts.py`
  - `phospy.science.transformations.transformers.identity.IdentityTransformer`
- `src/phospy/preprocessing/stages/sequence.py`
  - `phospy.science.sequences`
  - `phospy.science.sequences.resolver.v1`

### Legacy behavior markers

- `scientific_policy_identifiers`
  - `src/phospy/differential/moderation/empirical_bayes.py`
  - `src/phospy/kinase/activity/execution.py`
  - `src/phospy/kinase/activity/membership.py`
  - `src/phospy/kinase/activity/methods/ssgsea_substrate_enrichment.py`
  - `src/phospy/kinase/activity/models.py`
  - `src/phospy/kinase/activity/result_validation.py`
  - `src/phospy/kinase/activity/results.py`
  - `src/phospy/kinase/prediction/results.py`
  - `src/phospy/kinase/scoring/config.py`
  - `src/phospy/core/phosphosites/display_identity.py`
  - `src/phospy/preprocessing/control_sites/contracts.py`
  - `src/phospy/preprocessing/control_sites/discovery_contracts.py`
  - `src/phospy/preprocessing/control_sites/discovery_provenance.py`
  - `src/phospy/preprocessing/control_sites/discovery_results.py`
  - `src/phospy/preprocessing/sequence/reference.py`
  - `src/phospy/core/provenance/environment.py`
- `signalome_versioned_input_and_diagnostic_aliases`
  - `src/phospy/signalome/assignments.py`
  - `src/phospy/signalome/caveats.py`
  - `src/phospy/signalome/clustering/models.py`
  - `src/phospy/signalome/clustering/protein_modules.py`
  - `src/phospy/signalome/constants.py`
  - `src/phospy/signalome/contracts.py`
  - `src/phospy/signalome/models.py`
  - `src/phospy/signalome/network.py`
  - `src/phospy/signalome/protein_resolution.py`
  - `src/phospy/signalome/provenance.py`
  - `src/phospy/signalome/requests.py`
  - `src/phospy/signalome/row_attrition.py`
  - `src/phospy/signalome/validator.py`
- `versioned_schema_and_payload_compatibility`
  - `src/phospy/datasets/construction/input_resolution.py`
  - `src/phospy/datasets/trusted_assertions.py`
  - `src/phospy/enrichment/contracts/collections.py`
  - `src/phospy/core/evidence/__init__.py`
  - `src/phospy/core/evidence/dataset_resolution/__init__.py`
  - `src/phospy/core/evidence/dataset_resolution/models.py`
  - `src/phospy/core/evidence/serialization/historical_peptide.py`
  - `src/phospy/io/bundles/_kinase/manifest.py`
  - `src/phospy/io/bundles/_kinase/reconstruction.py`
  - `src/phospy/io/bundles/_kinase/snapshots.py`
  - `src/phospy/io/bundles/_shared/processing_state.py`
  - `src/phospy/io/bundles/_signalome/config.py`
  - `src/phospy/io/bundles/_signalome/diagnostics.py`
  - `src/phospy/io/bundles/_signalome/manifest.py`
  - `src/phospy/io/bundles/_signalome/reconstruction.py`
  - `src/phospy/io/bundles/_signalome/snapshots.py`
  - `src/phospy/io/bundles/_signalome/tables.py`
  - `src/phospy/preprocessing/provenance/serialization.py`
  - `src/phospy/preprocessing/results/quantitative_evidence.py`
  - `src/phospy/core/provenance/serialization/payload.py`
  - `src/phospy/core/provenance/serialization/environment.py`
  - `src/phospy/core/provenance/serialization/tables.py`
  - `src/phospy/core/quantitative/models.py`
  - `src/phospy/core/quantitative/policy.py`
  - `src/phospy/core/quantitative/scale_state.py`
  - `src/phospy/core/quantitative/serialization.py`

## Complete tracked production file inventory

- `src/phospy/__init__.py`
- `src/phospy/api/__init__.py`
- `src/phospy/api/_inventory.py`
- `src/phospy/api/advanced/__init__.py`
- `src/phospy/api/builders.py`
- `src/phospy/core/__init__.py`
- `src/phospy/core/activities/__init__.py`
- `src/phospy/core/activities/contracts.py`
- `src/phospy/core/activities/semantics.py`
- `src/phospy/core/activities/validation.py`
- `src/phospy/core/errors/__init__.py`
- `src/phospy/core/errors/base.py`
- `src/phospy/core/errors/build.py`
- `src/phospy/core/errors/input.py`
- `src/phospy/core/errors/provenance.py`
- `src/phospy/core/errors/references.py`
- `src/phospy/core/errors/transformations.py`
- `src/phospy/core/errors/validation.py`
- `src/phospy/core/errors/workflows.py`
- `src/phospy/core/evidence/__init__.py`
- `src/phospy/core/evidence/caveats.py`
- `src/phospy/core/evidence/config.py`
- `src/phospy/core/evidence/dataset_resolution/__init__.py`
- `src/phospy/core/evidence/dataset_resolution/allocation.py`
- `src/phospy/core/evidence/dataset_resolution/mapping.py`
- `src/phospy/core/evidence/dataset_resolution/models.py`
- `src/phospy/core/evidence/dataset_resolution/policies.py`
- `src/phospy/core/evidence/dataset_resolution/resolver.py`
- `src/phospy/core/evidence/dataset_resolution/site_metadata.py`
- `src/phospy/core/evidence/dataset_resolution/site_sequence.py`
- `src/phospy/core/evidence/dataset_resolution/summary.py`
- `src/phospy/core/evidence/localisation.py`
- `src/phospy/core/evidence/models.py`
- `src/phospy/core/evidence/modified_peptides.py`
- `src/phospy/core/evidence/multi_site.py`
- `src/phospy/core/evidence/quantification_depth.py`
- `src/phospy/core/evidence/sequence_resolution.py`
- `src/phospy/core/evidence/serialization/historical_peptide.py`
- `src/phospy/core/evidence/site_metadata.py`
- `src/phospy/core/frames/__init__.py`
- `src/phospy/core/frames/comparison.py`
- `src/phospy/core/frames/missing_values.py`
- `src/phospy/core/frames/numeric.py`
- `src/phospy/core/frames/ownership.py`
- `src/phospy/core/frames/pandas_typing.py`
- `src/phospy/core/frames/table_formats.py`
- `src/phospy/core/frames/table_schema.py`
- `src/phospy/core/frames/validation.py`
- `src/phospy/core/phosphosites/__init__.py`
- `src/phospy/core/phosphosites/contracts.py`
- `src/phospy/core/phosphosites/dataset_identity.py`
- `src/phospy/core/phosphosites/display_identity.py`
- `src/phospy/core/phosphosites/identifiers.py`
- `src/phospy/core/phosphosites/identity.py`
- `src/phospy/core/phosphosites/messages.py`
- `src/phospy/core/phosphosites/metadata_validation.py`
- `src/phospy/core/phosphosites/parsing.py`
- `src/phospy/core/phosphosites/result_identity.py`
- `src/phospy/core/phosphosites/sequence_context.py`
- `src/phospy/core/phosphosites/site_keys.py`
- `src/phospy/core/phosphosites/validation.py`
- `src/phospy/core/phosphosites/workflow_identity.py`
- `src/phospy/core/policies/__init__.py`
- `src/phospy/core/policies/config_values.py`
- `src/phospy/core/policies/policy_base.py`
- `src/phospy/core/provenance/__init__.py`
- `src/phospy/core/provenance/environment.py`
- `src/phospy/core/provenance/hashing.py`
- `src/phospy/core/provenance/immutability.py`
- `src/phospy/core/provenance/models/__init__.py`
- `src/phospy/core/provenance/models/batch_correction.py`
- `src/phospy/core/provenance/models/environment.py`
- `src/phospy/core/provenance/models/tables.py`
- `src/phospy/core/provenance/result_caveats.py`
- `src/phospy/core/provenance/row_attrition.py`
- `src/phospy/core/provenance/scalar_validation.py`
- `src/phospy/core/provenance/scientific_policy_models.py`
- `src/phospy/core/provenance/serialization/__init__.py`
- `src/phospy/core/provenance/serialization/batch_correction.py`
- `src/phospy/core/provenance/serialization/environment.py`
- `src/phospy/core/provenance/serialization/payload.py`
- `src/phospy/core/provenance/serialization/tables.py`
- `src/phospy/core/quantitative/__init__.py`
- `src/phospy/core/quantitative/_payloads.py`
- `src/phospy/core/quantitative/assertions.py`
- `src/phospy/core/quantitative/authority.py`
- `src/phospy/core/quantitative/contracts.py`
- `src/phospy/core/quantitative/declared_scale_diagnostics.py`
- `src/phospy/core/quantitative/evidence.py`
- `src/phospy/core/quantitative/models.py`
- `src/phospy/core/quantitative/policy.py`
- `src/phospy/core/quantitative/provenance.py`
- `src/phospy/core/quantitative/quantitative_contracts.py`
- `src/phospy/core/quantitative/scale_state.py`
- `src/phospy/core/quantitative/scale_values.py`
- `src/phospy/core/quantitative/serialization.py`
- `src/phospy/core/quantitative/state_coherence.py`
- `src/phospy/core/quantitative/transformers/__init__.py`
- `src/phospy/core/quantitative/transformers/identifiers.py`
- `src/phospy/core/quantitative/transformers/identity.py`
- `src/phospy/core/quantitative/transformers/log2.py`
- `src/phospy/core/quantitative/validation.py`
- `src/phospy/core/references/__init__.py`
- `src/phospy/core/references/acquisition.py`
- `src/phospy/core/references/builder.py`
- `src/phospy/core/references/builder_validation.py`
- `src/phospy/core/references/compatibility.py`
- `src/phospy/core/references/config.py`
- `src/phospy/core/references/context_compatibility.py`
- `src/phospy/core/references/errors.py`
- `src/phospy/core/references/fingerprints.py`
- `src/phospy/core/references/identifiers.py`
- `src/phospy/core/references/manifest.py`
- `src/phospy/core/references/manifest_common.py`
- `src/phospy/core/references/manifest_files.py`
- `src/phospy/core/references/manifest_model.py`
- `src/phospy/core/references/manifest_policy.py`
- `src/phospy/core/references/models.py`
- `src/phospy/core/references/provenance/context.py`
- `src/phospy/core/references/provenance/identifiers.py`
- `src/phospy/core/references/provenance/models.py`
- `src/phospy/core/references/provenance/serialization.py`
- `src/phospy/core/references/redistribution.py`
- `src/phospy/core/references/resolution.py`
- `src/phospy/core/references/resources.py`
- `src/phospy/core/references/sequences/__init__.py`
- `src/phospy/core/references/sequences/fasta.py`
- `src/phospy/core/references/sequences/models.py`
- `src/phospy/core/references/sequences/repository.py`
- `src/phospy/core/references/sequences/validation.py`
- `src/phospy/core/references/tables.py`
- `src/phospy/core/references/validation/__init__.py`
- `src/phospy/core/references/validation/_diagnostics.py`
- `src/phospy/core/references/validation/bundle.py`
- `src/phospy/core/references/validation/bundle_semantics.py`
- `src/phospy/core/references/validation/manifest_schema.py`
- `src/phospy/core/references/validation/redistribution_policy.py`
- `src/phospy/core/references/validation/resource_integrity.py`
- `src/phospy/core/statistics/__init__.py`
- `src/phospy/core/statistics/multiple_testing.py`
- `src/phospy/core/taxonomy/__init__.py`
- `src/phospy/core/taxonomy/organisms.py`
- `src/phospy/core/validation/__init__.py`
- `src/phospy/core/validation/common/__init__.py`
- `src/phospy/core/validation/common/config_values.py`
- `src/phospy/core/validation/common/input_sources.py`
- `src/phospy/core/validation/common/numbers.py`
- `src/phospy/core/validation/common/paths.py`
- `src/phospy/core/validation/frame_columns.py`
- `src/phospy/core/validation/ownership.py`
- `src/phospy/data/reference_bundles/rat/l6_native/ATTRIBUTION.md`
- `src/phospy/data/reference_bundles/rat/l6_native/manifest.json`
- `src/phospy/data/reference_bundles/rat/l6_native/motif_scores.csv`
- `src/phospy/data/reference_bundles/rat/l6_native/motif_sizes.csv`
- `src/phospy/data/reference_bundles/rat/l6_native/site_sequences.csv`
- `src/phospy/data/reference_bundles/rat/l6_native/substrate_map.csv`
- `src/phospy/datasets/__init__.py`
- `src/phospy/datasets/_construction/__init__.py`
- `src/phospy/datasets/_construction/fingerprints.py`
- `src/phospy/datasets/_construction/organism_coherence.py`
- `src/phospy/datasets/_construction/protein_aware_binding.py`
- `src/phospy/datasets/_construction/service.py`
- `src/phospy/datasets/_construction/trusted_assertions.py`
- `src/phospy/datasets/_construction/validation.py`
- `src/phospy/datasets/_internal/__init__.py`
- `src/phospy/datasets/_internal/frame_store.py`
- `src/phospy/datasets/_processing_state/__init__.py`
- `src/phospy/datasets/_processing_state/json_contracts.py`
- `src/phospy/datasets/_processing_state/missing_data.py`
- `src/phospy/datasets/_processing_state/models.py`
- `src/phospy/datasets/_processing_state/total_protein.py`
- `src/phospy/datasets/builder.py`
- `src/phospy/datasets/construction/__init__.py`
- `src/phospy/datasets/construction/contracts.py`
- `src/phospy/datasets/construction/executor.py`
- `src/phospy/datasets/construction/input_resolution.py`
- `src/phospy/datasets/construction/intensity_scale.py`
- `src/phospy/datasets/construction/interpreter.py`
- `src/phospy/datasets/construction/normalization.py`
- `src/phospy/datasets/construction/normalization_reporting.py`
- `src/phospy/datasets/construction/planning.py`
- `src/phospy/datasets/construction/preprocessing.py`
- `src/phospy/datasets/construction/protein_aware.py`
- `src/phospy/datasets/construction/provenance.py`
- `src/phospy/datasets/construction/report_assembly.py`
- `src/phospy/datasets/construction/sample_metadata.py`
- `src/phospy/datasets/construction/sequence_derivation.py`
- `src/phospy/datasets/construction/site_identity.py`
- `src/phospy/datasets/construction/site_metadata.py`
- `src/phospy/datasets/construction/site_sequence_validation.py`
- `src/phospy/datasets/construction/state_assembly.py`
- `src/phospy/datasets/construction/total_matrix.py`
- `src/phospy/datasets/construction/transformation_state.py`
- `src/phospy/datasets/contracts.py`
- `src/phospy/datasets/derived_quantitative.py`
- `src/phospy/datasets/direct_construction.py`
- `src/phospy/datasets/imputation_metadata.py`
- `src/phospy/datasets/internal_view.py`
- `src/phospy/datasets/models.py`
- `src/phospy/datasets/processing_state.py`
- `src/phospy/datasets/protein_aware.py`
- `src/phospy/datasets/schemas.py`
- `src/phospy/datasets/trusted_assertions.py`
- `src/phospy/datasets/trusted_construction_caveats.py`
- `src/phospy/datasets/validation/__init__.py`
- `src/phospy/datasets/validation/inputs.py`
- `src/phospy/datasets/validation/model_boundary.py`
- `src/phospy/datasets/validation/request.py`
- `src/phospy/differential/__init__.py`
- `src/phospy/differential/_payload_codec.py`
- `src/phospy/differential/aggregation/__init__.py`
- `src/phospy/differential/aggregation/executor.py`
- `src/phospy/differential/aggregation/experimental.py`
- `src/phospy/differential/aggregation/models.py`
- `src/phospy/differential/aggregation/scientific_policies.py`
- `src/phospy/differential/caveats.py`
- `src/phospy/differential/config.py`
- `src/phospy/differential/contracts.py`
- `src/phospy/differential/design/__init__.py`
- `src/phospy/differential/design/contrast_helpers.py`
- `src/phospy/differential/design/matrix_builder.py`
- `src/phospy/differential/design/models.py`
- `src/phospy/differential/diagnostics.py`
- `src/phospy/differential/execution/__init__.py`
- `src/phospy/differential/execution/compound_symmetry_gls.py`
- `src/phospy/differential/execution/duplicate_correlation.py`
- `src/phospy/differential/execution/fitting.py`
- `src/phospy/differential/execution/linear_model.py`
- `src/phospy/differential/execution/model_executor.py`
- `src/phospy/differential/execution/protein_covariate_adjusted.py`
- `src/phospy/differential/execution/workflow.py`
- `src/phospy/differential/internal_view.py`
- `src/phospy/differential/models/__init__.py`
- `src/phospy/differential/models/design.py`
- `src/phospy/differential/models/duplicate_correlation.py`
- `src/phospy/differential/models/fit.py`
- `src/phospy/differential/models/protein_aware.py`
- `src/phospy/differential/moderation/__init__.py`
- `src/phospy/differential/moderation/config.py`
- `src/phospy/differential/moderation/empirical_bayes.py`
- `src/phospy/differential/planning/__init__.py`
- `src/phospy/differential/planning/design.py`
- `src/phospy/differential/planning/eligibility.py`
- `src/phospy/differential/planning/imputation.py`
- `src/phospy/differential/planning/interpreter.py`
- `src/phospy/differential/planning/models.py`
- `src/phospy/differential/planning/protein_aware_eligibility.py`
- `src/phospy/differential/planning/protein_aware_inputs.py`
- `src/phospy/differential/planning/replicates.py`
- `src/phospy/differential/policies.py`
- `src/phospy/differential/provenance/__init__.py`
- `src/phospy/differential/provenance/assembly.py`
- `src/phospy/differential/provenance/lineage.py`
- `src/phospy/differential/provenance/models.py`
- `src/phospy/differential/public.py`
- `src/phospy/differential/quantification_depth.py`
- `src/phospy/differential/result_assembly.py`
- `src/phospy/differential/results.py`
- `src/phospy/differential/tables.py`
- `src/phospy/differential/validation/__init__.py`
- `src/phospy/differential/validation/computation.py`
- `src/phospy/differential/validation/dataset.py`
- `src/phospy/differential/validation/design.py`
- `src/phospy/differential/validation/reliability.py`
- `src/phospy/differential/validation/result_tables.py`
- `src/phospy/differential/validation/workflow.py`
- `src/phospy/enrichment/__init__.py`
- `src/phospy/enrichment/assembly.py`
- `src/phospy/enrichment/caveats.py`
- `src/phospy/enrichment/contracts/__init__.py`
- `src/phospy/enrichment/contracts/collections.py`
- `src/phospy/enrichment/contracts/config.py`
- `src/phospy/enrichment/contracts/identifier_sets.py`
- `src/phospy/enrichment/contracts/requests.py`
- `src/phospy/enrichment/contracts/results.py`
- `src/phospy/enrichment/contracts/vocabulary.py`
- `src/phospy/enrichment/execution.py`
- `src/phospy/enrichment/filtering.py`
- `src/phospy/enrichment/interpreter.py`
- `src/phospy/enrichment/ora.py`
- `src/phospy/enrichment/ora_execution.py`
- `src/phospy/enrichment/planning.py`
- `src/phospy/enrichment/provenance.py`
- `src/phospy/enrichment/validation.py`
- `src/phospy/enrichment/workflow.py`
- `src/phospy/io/__init__.py`
- `src/phospy/io/bundles/__init__.py`
- `src/phospy/io/bundles/_kinase/__init__.py`
- `src/phospy/io/bundles/_kinase/constants.py`
- `src/phospy/io/bundles/_kinase/loader.py`
- `src/phospy/io/bundles/_kinase/manifest.py`
- `src/phospy/io/bundles/_kinase/models.py`
- `src/phospy/io/bundles/_kinase/reconstruction.py`
- `src/phospy/io/bundles/_kinase/snapshots.py`
- `src/phospy/io/bundles/_kinase/writer.py`
- `src/phospy/io/bundles/_shared/__init__.py`
- `src/phospy/io/bundles/_shared/integrity.py`
- `src/phospy/io/bundles/_shared/json_files.py`
- `src/phospy/io/bundles/_shared/kinase_scoring_state.py`
- `src/phospy/io/bundles/_shared/organisms.py`
- `src/phospy/io/bundles/_shared/paths.py`
- `src/phospy/io/bundles/_shared/primitives.py`
- `src/phospy/io/bundles/_shared/processing_state.py`
- `src/phospy/io/bundles/_shared/tables.py`
- `src/phospy/io/bundles/_shared/transactions.py`
- `src/phospy/io/bundles/_shared/trusted_dataset_assertions.py`
- `src/phospy/io/bundles/_signalome/__init__.py`
- `src/phospy/io/bundles/_signalome/config.py`
- `src/phospy/io/bundles/_signalome/constants.py`
- `src/phospy/io/bundles/_signalome/diagnostics.py`
- `src/phospy/io/bundles/_signalome/loader.py`
- `src/phospy/io/bundles/_signalome/manifest.py`
- `src/phospy/io/bundles/_signalome/models.py`
- `src/phospy/io/bundles/_signalome/primitives.py`
- `src/phospy/io/bundles/_signalome/reconstruction.py`
- `src/phospy/io/bundles/_signalome/snapshots.py`
- `src/phospy/io/bundles/_signalome/tables.py`
- `src/phospy/io/bundles/_signalome/writer.py`
- `src/phospy/io/bundles/kinase.py`
- `src/phospy/io/bundles/kinase_library.py`
- `src/phospy/io/bundles/signalome.py`
- `src/phospy/io/importers/__init__.py`
- `src/phospy/io/importers/contracts.py`
- `src/phospy/io/importers/fragpipe_validation.py`
- `src/phospy/io/importers/maxquant_validation.py`
- `src/phospy/io/importers/validation.py`
- `src/phospy/io/publishers/__init__.py`
- `src/phospy/io/publishers/workflows.py`
- `src/phospy/io/readers/__init__.py`
- `src/phospy/io/readers/_table_parsing.py`
- `src/phospy/io/readers/dataset_inputs.py`
- `src/phospy/io/readers/enrichment_sets.py`
- `src/phospy/io/readers/fragpipe/__init__.py`
- `src/phospy/io/readers/fragpipe/columns.py`
- `src/phospy/io/readers/fragpipe/constants.py`
- `src/phospy/io/readers/fragpipe/conversion.py`
- `src/phospy/io/readers/fragpipe/filtering.py`
- `src/phospy/io/readers/fragpipe/importer.py`
- `src/phospy/io/readers/fragpipe/models.py`
- `src/phospy/io/readers/fragpipe/normalization.py`
- `src/phospy/io/readers/fragpipe/raw.py`
- `src/phospy/io/readers/fragpipe/reporting.py`
- `src/phospy/io/readers/importers.py`
- `src/phospy/io/readers/maxquant.py`
- `src/phospy/io/readers/tables.py`
- `src/phospy/kinase/__init__.py`
- `src/phospy/kinase/_matrix_validation.py`
- `src/phospy/kinase/activity/__init__.py`
- `src/phospy/kinase/activity/config.py`
- `src/phospy/kinase/activity/contracts.py`
- `src/phospy/kinase/activity/diagnostics.py`
- `src/phospy/kinase/activity/execution.py`
- `src/phospy/kinase/activity/inputs.py`
- `src/phospy/kinase/activity/membership.py`
- `src/phospy/kinase/activity/membership_selection.py`
- `src/phospy/kinase/activity/methods/__init__.py`
- `src/phospy/kinase/activity/methods/ksea_zscore.py`
- `src/phospy/kinase/activity/methods/ssgsea_substrate_enrichment.py`
- `src/phospy/kinase/activity/methods/weighted_substrate_activity.py`
- `src/phospy/kinase/activity/models.py`
- `src/phospy/kinase/activity/policies.py`
- `src/phospy/kinase/activity/result_validation.py`
- `src/phospy/kinase/activity/results.py`
- `src/phospy/kinase/activity/runner.py`
- `src/phospy/kinase/activity/stage_contracts.py`
- `src/phospy/kinase/activity/statistics.py`
- `src/phospy/kinase/activity/tables.py`
- `src/phospy/kinase/activity/threshold_membership.py`
- `src/phospy/kinase/activity/validation.py`
- `src/phospy/kinase/attrition_metrics.py`
- `src/phospy/kinase/caveats.py`
- `src/phospy/kinase/contracts.py`
- `src/phospy/kinase/eligibility.py`
- `src/phospy/kinase/executor.py`
- `src/phospy/kinase/interpreter.py`
- `src/phospy/kinase/policies.py`
- `src/phospy/kinase/prediction/__init__.py`
- `src/phospy/kinase/prediction/candidates.py`
- `src/phospy/kinase/prediction/config.py`
- `src/phospy/kinase/prediction/diagnostics.py`
- `src/phospy/kinase/prediction/execution.py`
- `src/phospy/kinase/prediction/internal_view.py`
- `src/phospy/kinase/prediction/outputs.py`
- `src/phospy/kinase/prediction/policies.py`
- `src/phospy/kinase/prediction/results.py`
- `src/phospy/kinase/prediction/runner.py`
- `src/phospy/kinase/prediction/sampling_core.py`
- `src/phospy/kinase/prediction/sampling_runtime.py`
- `src/phospy/kinase/prediction/svm.py`
- `src/phospy/kinase/prediction/tables.py`
- `src/phospy/kinase/provenance.py`
- `src/phospy/kinase/public.py`
- `src/phospy/kinase/references/__init__.py`
- `src/phospy/kinase/references/loading.py`
- `src/phospy/kinase/references/models.py`
- `src/phospy/kinase/references/projection.py`
- `src/phospy/kinase/references/projection_contracts.py`
- `src/phospy/kinase/references/provenance.py`
- `src/phospy/kinase/references/validation.py`
- `src/phospy/kinase/requests.py`
- `src/phospy/kinase/resolved_validator.py`
- `src/phospy/kinase/result_assembly.py`
- `src/phospy/kinase/result_boundary.py`
- `src/phospy/kinase/results.py`
- `src/phospy/kinase/row_attrition.py`
- `src/phospy/kinase/scoring/__init__.py`
- `src/phospy/kinase/scoring/config.py`
- `src/phospy/kinase/scoring/contracts.py`
- `src/phospy/kinase/scoring/contributions.py`
- `src/phospy/kinase/scoring/execution.py`
- `src/phospy/kinase/scoring/kinase_library.py`
- `src/phospy/kinase/scoring/motif/__init__.py`
- `src/phospy/kinase/scoring/motif/frequency_matrices.py`
- `src/phospy/kinase/scoring/motif/kinase_library.py`
- `src/phospy/kinase/scoring/motif/library_validation.py`
- `src/phospy/kinase/scoring/motif/models.py`
- `src/phospy/kinase/scoring/motif/scaling.py`
- `src/phospy/kinase/scoring/motif/scoring.py`
- `src/phospy/kinase/scoring/policies.py`
- `src/phospy/kinase/scoring/profile.py`
- `src/phospy/kinase/scoring/rank_fusion.py`
- `src/phospy/kinase/scoring/results.py`
- `src/phospy/kinase/scoring/selection.py`
- `src/phospy/kinase/scoring/sequence_validation.py`
- `src/phospy/kinase/scoring/tables.py`
- `src/phospy/kinase/scoring/transforms.py`
- `src/phospy/kinase/sequence_contracts.py`
- `src/phospy/kinase/site_attrition.py`
- `src/phospy/kinase/site_sequence_policy.py`
- `src/phospy/kinase/site_sequence_support.py`
- `src/phospy/kinase/validator.py`
- `src/phospy/preparation_records/__init__.py`
- `src/phospy/preparation_records/batch_correction.py`
- `src/phospy/preparation_records/group_evidence.py`
- `src/phospy/preparation_records/imputation.py`
- `src/phospy/preparation_records/policies.py`
- `src/phospy/preparation_records/protein_aware.py`
- `src/phospy/preparation_records/report_schema.py`
- `src/phospy/preparation_records/validation/__init__.py`
- `src/phospy/preparation_records/validation/sample_groups.py`
- `src/phospy/preprocessing/__init__.py`
- `src/phospy/preprocessing/batch_correction/__init__.py`
- `src/phospy/preprocessing/batch_correction/adapters/__init__.py`
- `src/phospy/preprocessing/batch_correction/adapters/preprocessing.py`
- `src/phospy/preprocessing/batch_correction/contracts/__init__.py`
- `src/phospy/preprocessing/batch_correction/contracts/_validation.py`
- `src/phospy/preprocessing/batch_correction/contracts/configuration.py`
- `src/phospy/preprocessing/batch_correction/contracts/missingness.py`
- `src/phospy/preprocessing/batch_correction/contracts/output.py`
- `src/phospy/preprocessing/batch_correction/contracts/request.py`
- `src/phospy/preprocessing/batch_correction/contracts/results.py`
- `src/phospy/preprocessing/batch_correction/contracts/stage.py`
- `src/phospy/preprocessing/batch_correction/contracts/workflow.py`
- `src/phospy/preprocessing/batch_correction/environment.py`
- `src/phospy/preprocessing/batch_correction/execution/__init__.py`
- `src/phospy/preprocessing/batch_correction/execution/linear.py`
- `src/phospy/preprocessing/batch_correction/execution/missingness.py`
- `src/phospy/preprocessing/batch_correction/execution/ruv_iii.py`
- `src/phospy/preprocessing/batch_correction/execution/ruv_iii_executor.py`
- `src/phospy/preprocessing/batch_correction/execution/sps_ruv.py`
- `src/phospy/preprocessing/batch_correction/execution/stage.py`
- `src/phospy/preprocessing/batch_correction/metadata.py`
- `src/phospy/preprocessing/batch_correction/planning.py`
- `src/phospy/preprocessing/batch_correction/provenance/__init__.py`
- `src/phospy/preprocessing/batch_correction/provenance/_payloads.py`
- `src/phospy/preprocessing/batch_correction/provenance/assembly.py`
- `src/phospy/preprocessing/batch_correction/provenance/validation.py`
- `src/phospy/preprocessing/batch_correction/validation/__init__.py`
- `src/phospy/preprocessing/batch_correction/validation/_helpers.py`
- `src/phospy/preprocessing/batch_correction/validation/controls.py`
- `src/phospy/preprocessing/batch_correction/validation/design.py`
- `src/phospy/preprocessing/batch_correction/validation/feasibility.py`
- `src/phospy/preprocessing/batch_correction/validation/missingness.py`
- `src/phospy/preprocessing/batch_correction/validation/output.py`
- `src/phospy/preprocessing/batch_correction/validation/protocols.py`
- `src/phospy/preprocessing/batch_correction/validation/request.py`
- `src/phospy/preprocessing/batch_correction/validation/stage_order.py`
- `src/phospy/preprocessing/batch_correction/workflow.py`
- `src/phospy/preprocessing/contracts/__init__.py`
- `src/phospy/preprocessing/contracts/configuration/__init__.py`
- `src/phospy/preprocessing/contracts/configuration/comparisons.py`
- `src/phospy/preprocessing/contracts/configuration/coverage_filter.py`
- `src/phospy/preprocessing/contracts/configuration/dataset.py`
- `src/phospy/preprocessing/contracts/configuration/intensity_transform.py`
- `src/phospy/preprocessing/contracts/configuration/localisation.py`
- `src/phospy/preprocessing/contracts/configuration/missing_data.py`
- `src/phospy/preprocessing/contracts/configuration/normalisation.py`
- `src/phospy/preprocessing/contracts/configuration/site_matrix.py`
- `src/phospy/preprocessing/contracts/configuration/site_sequence.py`
- `src/phospy/preprocessing/contracts/configuration/total_protein.py`
- `src/phospy/preprocessing/contracts/imputation_scale.py`
- `src/phospy/preprocessing/contracts/plan.py`
- `src/phospy/preprocessing/contracts/plan_resolution.py`
- `src/phospy/preprocessing/contracts/policies.py`
- `src/phospy/preprocessing/contracts/stage.py`
- `src/phospy/preprocessing/contracts/stage_keys.py`
- `src/phospy/preprocessing/contracts/state.py`
- `src/phospy/preprocessing/control_sites/__init__.py`
- `src/phospy/preprocessing/control_sites/_discovery_payload.py`
- `src/phospy/preprocessing/control_sites/contracts.py`
- `src/phospy/preprocessing/control_sites/discovery_contracts.py`
- `src/phospy/preprocessing/control_sites/discovery_execution.py`
- `src/phospy/preprocessing/control_sites/discovery_planning.py`
- `src/phospy/preprocessing/control_sites/discovery_provenance.py`
- `src/phospy/preprocessing/control_sites/discovery_ranking.py`
- `src/phospy/preprocessing/control_sites/discovery_results.py`
- `src/phospy/preprocessing/control_sites/discovery_validation.py`
- `src/phospy/preprocessing/control_sites/validation.py`
- `src/phospy/preprocessing/control_sites/workflow.py`
- `src/phospy/preprocessing/dataset_construction.py`
- `src/phospy/preprocessing/evidence/missing_data.py`
- `src/phospy/preprocessing/execution/corrected_output.py`
- `src/phospy/preprocessing/execution/diagnostics_normalization.py`
- `src/phospy/preprocessing/execution/fingerprints.py`
- `src/phospy/preprocessing/execution/pipeline.py`
- `src/phospy/preprocessing/execution/trace_builder.py`
- `src/phospy/preprocessing/planning/assembly.py`
- `src/phospy/preprocessing/planning/boundary.py`
- `src/phospy/preprocessing/planning/configuration.py`
- `src/phospy/preprocessing/planning/ordering.py`
- `src/phospy/preprocessing/planning/registry.py`
- `src/phospy/preprocessing/planning/rules.py`
- `src/phospy/preprocessing/planning/validation.py`
- `src/phospy/preprocessing/protein_aware/__init__.py`
- `src/phospy/preprocessing/protein_aware/alignment.py`
- `src/phospy/preprocessing/protein_aware/identity.py`
- `src/phospy/preprocessing/protein_aware/mapping.py`
- `src/phospy/preprocessing/provenance/composition.py`
- `src/phospy/preprocessing/provenance/models.py`
- `src/phospy/preprocessing/provenance/policies.py`
- `src/phospy/preprocessing/provenance/serialization.py`
- `src/phospy/preprocessing/results/diagnostics.py`
- `src/phospy/preprocessing/results/execution.py`
- `src/phospy/preprocessing/results/missing_data.py`
- `src/phospy/preprocessing/results/protein_aware.py`
- `src/phospy/preprocessing/results/quantitative_evidence.py`
- `src/phospy/preprocessing/results/report_rows.py`
- `src/phospy/preprocessing/results/report_schema.py`
- `src/phospy/preprocessing/sequence/__init__.py`
- `src/phospy/preprocessing/sequence/conflicts.py`
- `src/phospy/preprocessing/sequence/diagnostics.py`
- `src/phospy/preprocessing/sequence/metadata.py`
- `src/phospy/preprocessing/sequence/models.py`
- `src/phospy/preprocessing/sequence/reference.py`
- `src/phospy/preprocessing/sequence/request.py`
- `src/phospy/preprocessing/stages/__init__.py`
- `src/phospy/preprocessing/stages/batch_correction.py`
- `src/phospy/preprocessing/stages/comparisons.py`
- `src/phospy/preprocessing/stages/group_coverage.py`
- `src/phospy/preprocessing/stages/intensity.py`
- `src/phospy/preprocessing/stages/localisation.py`
- `src/phospy/preprocessing/stages/missing_data/__init__.py`
- `src/phospy/preprocessing/stages/missing_data/audit.py`
- `src/phospy/preprocessing/stages/missing_data/diagnostics.py`
- `src/phospy/preprocessing/stages/missing_data/forbid.py`
- `src/phospy/preprocessing/stages/missing_data/group_aware_routing.py`
- `src/phospy/preprocessing/stages/missing_data/knn.py`
- `src/phospy/preprocessing/stages/missing_data/minprob.py`
- `src/phospy/preprocessing/stages/missing_data/models.py`
- `src/phospy/preprocessing/stages/missing_data/row_median.py`
- `src/phospy/preprocessing/stages/missing_data/stage.py`
- `src/phospy/preprocessing/stages/normalisation.py`
- `src/phospy/preprocessing/stages/protein_aware.py`
- `src/phospy/preprocessing/stages/sequence.py`
- `src/phospy/preprocessing/stages/site_matrix.py`
- `src/phospy/preprocessing/stages/site_matrix_components/__init__.py`
- `src/phospy/preprocessing/stages/site_matrix_components/audit.py`
- `src/phospy/preprocessing/stages/site_matrix_components/duplicates.py`
- `src/phospy/preprocessing/stages/site_matrix_components/identity.py`
- `src/phospy/preprocessing/stages/site_matrix_components/metadata.py`
- `src/phospy/preprocessing/stages/total_protein.py`
- `src/phospy/preprocessing/validation/configuration.py`
- `src/phospy/preprocessing/validation/configuration_fields.py`
- `src/phospy/preprocessing/validation/group_coverage.py`
- `src/phospy/preprocessing/validation/quantitative_evidence.py`
- `src/phospy/preprocessing/validation/quantitative_transition.py`
- `src/phospy/preprocessing/workflow.py`
- `src/phospy/run_lineage/__init__.py`
- `src/phospy/run_lineage/models.py`
- `src/phospy/run_lineage/serialization.py`
- `src/phospy/signalome/__init__.py`
- `src/phospy/signalome/alignment_diagnostics.py`
- `src/phospy/signalome/assignments.py`
- `src/phospy/signalome/boundary_errors.py`
- `src/phospy/signalome/caveats.py`
- `src/phospy/signalome/clustering/__init__.py`
- `src/phospy/signalome/clustering/backend_dispatch.py`
- `src/phospy/signalome/clustering/backends/__init__.py`
- `src/phospy/signalome/clustering/backends/exact_python.py`
- `src/phospy/signalome/clustering/backends/scipy_hierarchical.py`
- `src/phospy/signalome/clustering/candidate_scoring.py`
- `src/phospy/signalome/clustering/candidate_selection.py`
- `src/phospy/signalome/clustering/contracts.py`
- `src/phospy/signalome/clustering/diagnostic_schemas.py`
- `src/phospy/signalome/clustering/diagnostics.py`
- `src/phospy/signalome/clustering/exact_python.py`
- `src/phospy/signalome/clustering/models.py`
- `src/phospy/signalome/clustering/module_selection.py`
- `src/phospy/signalome/clustering/orchestration.py`
- `src/phospy/signalome/clustering/policies.py`
- `src/phospy/signalome/clustering/protein_modules.py`
- `src/phospy/signalome/clustering/scale_guards.py`
- `src/phospy/signalome/clustering/scientific_policies.py`
- `src/phospy/signalome/clustering/scipy_hierarchical.py`
- `src/phospy/signalome/clustering/scoring.py`
- `src/phospy/signalome/clustering/selection.py`
- `src/phospy/signalome/clustering/stability.py`
- `src/phospy/signalome/clustering/tree_building.py`
- `src/phospy/signalome/clustering/tree_engine_adapter.py`
- `src/phospy/signalome/clustering/validation.py`
- `src/phospy/signalome/clustering_runner.py`
- `src/phospy/signalome/component_helpers.py`
- `src/phospy/signalome/component_models.py`
- `src/phospy/signalome/config.py`
- `src/phospy/signalome/constants.py`
- `src/phospy/signalome/context.py`
- `src/phospy/signalome/context_tables.py`
- `src/phospy/signalome/contracts.py`
- `src/phospy/signalome/errors.py`
- `src/phospy/signalome/executor.py`
- `src/phospy/signalome/expanded.py`
- `src/phospy/signalome/interpreter.py`
- `src/phospy/signalome/matrix_alignment.py`
- `src/phospy/signalome/models.py`
- `src/phospy/signalome/module_tables.py`
- `src/phospy/signalome/modules.py`
- `src/phospy/signalome/network.py`
- `src/phospy/signalome/network_builder.py`
- `src/phospy/signalome/protein_resolution.py`
- `src/phospy/signalome/provenance.py`
- `src/phospy/signalome/public.py`
- `src/phospy/signalome/requests.py`
- `src/phospy/signalome/result_assembly.py`
- `src/phospy/signalome/result_validation.py`
- `src/phospy/signalome/results.py`
- `src/phospy/signalome/row_attrition.py`
- `src/phospy/signalome/scientific_policies.py`
- `src/phospy/signalome/score_matrix_selection.py`
- `src/phospy/signalome/score_preconditioning.py`
- `src/phospy/signalome/tables/__init__.py`
- `src/phospy/signalome/tables/assignments.py`
- `src/phospy/signalome/tables/common.py`
- `src/phospy/signalome/tables/context.py`
- `src/phospy/signalome/tables/modules.py`
- `src/phospy/signalome/tables/network.py`
- `src/phospy/signalome/validator.py`

## Complete Python module ownership inventory

| Source path | Module | Owner | Rationale |
| --- | --- | --- | --- |
| `src/phospy/__init__.py` | `phospy` | `public_facades` | Metadata-only package root |
| `src/phospy/api/__init__.py` | `phospy.api` | `public_facades` | Supported or compatibility aggregate adapter |
| `src/phospy/api/_inventory.py` | `phospy.api._inventory` | `public_facades` | Supported aggregate inventory |
| `src/phospy/api/advanced/__init__.py` | `phospy.api.advanced` | `public_facades` | Supported specialist aggregate adapter |
| `src/phospy/api/builders.py` | `phospy.api.builders` | `public_facades` | Internal API composition route; unsupported for direct consumer import |
| `src/phospy/core/__init__.py` | `phospy.core` | `namespace_containers` | Inert organisational namespace initializer |
| `src/phospy/core/activities/__init__.py` | `phospy.core.activities` | `activities` | Passive activity quantitative contract package |
| `src/phospy/core/activities/contracts.py` | `phospy.core.activities.contracts` | `activities` | Generic activity quantitative method contract |
| `src/phospy/core/activities/semantics.py` | `phospy.core.activities.semantics` | `activities` | Generic activity profile and quantitative semantics |
| `src/phospy/core/activities/validation.py` | `phospy.core.activities.validation` | `activities` | Generic activity quantitative contract validation |
| `src/phospy/core/errors/__init__.py` | `phospy.core.errors` | `errors` | Shared exception taxonomy |
| `src/phospy/core/errors/base.py` | `phospy.core.errors.base` | `errors` | Shared exception taxonomy |
| `src/phospy/core/errors/build.py` | `phospy.core.errors.build` | `errors` | Shared exception taxonomy |
| `src/phospy/core/errors/input.py` | `phospy.core.errors.input` | `errors` | Shared exception taxonomy |
| `src/phospy/core/errors/provenance.py` | `phospy.core.errors.provenance` | `errors` | Shared exception taxonomy |
| `src/phospy/core/errors/references.py` | `phospy.core.errors.references` | `errors` | Shared exception taxonomy |
| `src/phospy/core/errors/transformations.py` | `phospy.core.errors.transformations` | `errors` | Shared exception taxonomy |
| `src/phospy/core/errors/validation.py` | `phospy.core.errors.validation` | `errors` | Shared exception taxonomy |
| `src/phospy/core/errors/workflows.py` | `phospy.core.errors.workflows` | `errors` | Shared exception taxonomy |
| `src/phospy/core/evidence/__init__.py` | `phospy.core.evidence` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/caveats.py` | `phospy.core.evidence.caveats` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/config.py` | `phospy.core.evidence.config` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/dataset_resolution/__init__.py` | `phospy.core.evidence.dataset_resolution` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/dataset_resolution/allocation.py` | `phospy.core.evidence.dataset_resolution.allocation` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/dataset_resolution/mapping.py` | `phospy.core.evidence.dataset_resolution.mapping` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/dataset_resolution/models.py` | `phospy.core.evidence.dataset_resolution.models` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/dataset_resolution/policies.py` | `phospy.core.evidence.dataset_resolution.policies` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/dataset_resolution/resolver.py` | `phospy.core.evidence.dataset_resolution.resolver` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/dataset_resolution/site_metadata.py` | `phospy.core.evidence.dataset_resolution.site_metadata` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/dataset_resolution/site_sequence.py` | `phospy.core.evidence.dataset_resolution.site_sequence` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/dataset_resolution/summary.py` | `phospy.core.evidence.dataset_resolution.summary` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/localisation.py` | `phospy.core.evidence.localisation` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/models.py` | `phospy.core.evidence.models` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/modified_peptides.py` | `phospy.core.evidence.modified_peptides` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/multi_site.py` | `phospy.core.evidence.multi_site` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/quantification_depth.py` | `phospy.core.evidence.quantification_depth` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/sequence_resolution.py` | `phospy.core.evidence.sequence_resolution` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/evidence/serialization/historical_peptide.py` | `phospy.core.evidence.serialization.historical_peptide` | `evidence` | Historical peptide evidence compatibility codec |
| `src/phospy/core/evidence/site_metadata.py` | `phospy.core.evidence.site_metadata` | `evidence` | Evidence policy, records, allocation, or structural resolution |
| `src/phospy/core/frames/__init__.py` | `phospy.core.frames` | `frames` | Generic dataframe ownership or schema mechanics |
| `src/phospy/core/frames/comparison.py` | `phospy.core.frames.comparison` | `frames` | Generic dataframe ownership or schema mechanics |
| `src/phospy/core/frames/missing_values.py` | `phospy.core.frames.missing_values` | `frames` | Generic dataframe ownership or schema mechanics |
| `src/phospy/core/frames/numeric.py` | `phospy.core.frames.numeric` | `frames` | Generic dataframe ownership or schema mechanics |
| `src/phospy/core/frames/ownership.py` | `phospy.core.frames.ownership` | `frames` | Generic dataframe ownership or schema mechanics |
| `src/phospy/core/frames/pandas_typing.py` | `phospy.core.frames.pandas_typing` | `frames` | Generic dataframe ownership or schema mechanics |
| `src/phospy/core/frames/table_formats.py` | `phospy.core.frames.table_formats` | `frames` | Generic dataframe ownership or schema mechanics |
| `src/phospy/core/frames/table_schema.py` | `phospy.core.frames.table_schema` | `frames` | Generic dataframe ownership or schema mechanics |
| `src/phospy/core/frames/validation.py` | `phospy.core.frames.validation` | `frames` | Generic dataframe ownership or schema mechanics |
| `src/phospy/core/phosphosites/__init__.py` | `phospy.core.phosphosites` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/phosphosites/contracts.py` | `phospy.core.phosphosites.contracts` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/phosphosites/dataset_identity.py` | `phospy.core.phosphosites.dataset_identity` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/phosphosites/display_identity.py` | `phospy.core.phosphosites.display_identity` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/phosphosites/identifiers.py` | `phospy.core.phosphosites.identifiers` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/phosphosites/identity.py` | `phospy.core.phosphosites.identity` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/phosphosites/messages.py` | `phospy.core.phosphosites.messages` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/phosphosites/metadata_validation.py` | `phospy.core.phosphosites.metadata_validation` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/phosphosites/parsing.py` | `phospy.core.phosphosites.parsing` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/phosphosites/result_identity.py` | `phospy.core.phosphosites.result_identity` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/phosphosites/sequence_context.py` | `phospy.core.phosphosites.sequence_context` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/phosphosites/site_keys.py` | `phospy.core.phosphosites.site_keys` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/phosphosites/validation.py` | `phospy.core.phosphosites.validation` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/phosphosites/workflow_identity.py` | `phospy.core.phosphosites.workflow_identity` | `phosphosites` | Shared phosphosite identity or supplied-window rule |
| `src/phospy/core/policies/__init__.py` | `phospy.core.policies` | `policies` | Generic parsing and scalar policy |
| `src/phospy/core/policies/config_values.py` | `phospy.core.policies.config_values` | `policies` | Generic parsing and scalar policy |
| `src/phospy/core/policies/policy_base.py` | `phospy.core.policies.policy_base` | `policies` | Generic parsing and scalar policy |
| `src/phospy/core/provenance/__init__.py` | `phospy.core.provenance` | `provenance` | Generic immutable provenance infrastructure |
| `src/phospy/core/provenance/environment.py` | `phospy.core.provenance.environment` | `provenance` | Generic immutable provenance infrastructure |
| `src/phospy/core/provenance/hashing.py` | `phospy.core.provenance.hashing` | `provenance` | Generic immutable provenance infrastructure |
| `src/phospy/core/provenance/immutability.py` | `phospy.core.provenance.immutability` | `provenance` | Generic immutable provenance infrastructure |
| `src/phospy/core/provenance/models/__init__.py` | `phospy.core.provenance.models` | `provenance` | Generic immutable provenance infrastructure |
| `src/phospy/core/provenance/models/batch_correction.py` | `phospy.core.provenance.models.batch_correction` | `provenance` | Immutable batch-correction audit wire contract |
| `src/phospy/core/provenance/models/environment.py` | `phospy.core.provenance.models.environment` | `provenance` | Generic immutable provenance infrastructure |
| `src/phospy/core/provenance/models/tables.py` | `phospy.core.provenance.models.tables` | `provenance` | Generic immutable provenance infrastructure |
| `src/phospy/core/provenance/result_caveats.py` | `phospy.core.provenance.result_caveats` | `provenance` | Generic immutable provenance infrastructure |
| `src/phospy/core/provenance/row_attrition.py` | `phospy.core.provenance.row_attrition` | `provenance` | Generic immutable provenance infrastructure |
| `src/phospy/core/provenance/scalar_validation.py` | `phospy.core.provenance.scalar_validation` | `provenance` | Owner-internal provenance scalar-validation contract |
| `src/phospy/core/provenance/scientific_policy_models.py` | `phospy.core.provenance.scientific_policy_models` | `provenance` | Generic immutable provenance infrastructure |
| `src/phospy/core/provenance/serialization/__init__.py` | `phospy.core.provenance.serialization` | `provenance` | Generic immutable provenance infrastructure |
| `src/phospy/core/provenance/serialization/batch_correction.py` | `phospy.core.provenance.serialization.batch_correction` | `provenance` | Direct historical batch-correction provenance codec |
| `src/phospy/core/provenance/serialization/environment.py` | `phospy.core.provenance.serialization.environment` | `provenance` | Generic immutable provenance infrastructure |
| `src/phospy/core/provenance/serialization/payload.py` | `phospy.core.provenance.serialization.payload` | `provenance` | Owner-internal provenance payload coercion contract |
| `src/phospy/core/provenance/serialization/tables.py` | `phospy.core.provenance.serialization.tables` | `provenance` | Generic immutable provenance infrastructure |
| `src/phospy/core/quantitative/__init__.py` | `phospy.core.quantitative` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/_payloads.py` | `phospy.core.quantitative._payloads` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/assertions.py` | `phospy.core.quantitative.assertions` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/authority.py` | `phospy.core.quantitative.authority` | `quantitative` | Owner-internal quantitative authority contract |
| `src/phospy/core/quantitative/contracts.py` | `phospy.core.quantitative.contracts` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/declared_scale_diagnostics.py` | `phospy.core.quantitative.declared_scale_diagnostics` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/evidence.py` | `phospy.core.quantitative.evidence` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/models.py` | `phospy.core.quantitative.models` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/policy.py` | `phospy.core.quantitative.policy` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/provenance.py` | `phospy.core.quantitative.provenance` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/quantitative_contracts.py` | `phospy.core.quantitative.quantitative_contracts` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/scale_state.py` | `phospy.core.quantitative.scale_state` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/scale_values.py` | `phospy.core.quantitative.scale_values` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/serialization.py` | `phospy.core.quantitative.serialization` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/state_coherence.py` | `phospy.core.quantitative.state_coherence` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/transformers/__init__.py` | `phospy.core.quantitative.transformers` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/transformers/identifiers.py` | `phospy.core.quantitative.transformers.identifiers` | `quantitative` | Owner-internal quantitative transformer identifier contract |
| `src/phospy/core/quantitative/transformers/identity.py` | `phospy.core.quantitative.transformers.identity` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/transformers/log2.py` | `phospy.core.quantitative.transformers.log2` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/quantitative/validation.py` | `phospy.core.quantitative.validation` | `quantitative` | Shared quantitative state, transformation, validation, or serialization contract |
| `src/phospy/core/references/__init__.py` | `phospy.core.references` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/acquisition.py` | `phospy.core.references.acquisition` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/builder.py` | `phospy.core.references.builder` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/builder_validation.py` | `phospy.core.references.builder_validation` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/compatibility.py` | `phospy.core.references.compatibility` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/config.py` | `phospy.core.references.config` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/context_compatibility.py` | `phospy.core.references.context_compatibility` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/errors.py` | `phospy.core.references.errors` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/fingerprints.py` | `phospy.core.references.fingerprints` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/identifiers.py` | `phospy.core.references.identifiers` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/manifest.py` | `phospy.core.references.manifest` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/manifest_common.py` | `phospy.core.references.manifest_common` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/manifest_files.py` | `phospy.core.references.manifest_files` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/manifest_model.py` | `phospy.core.references.manifest_model` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/manifest_policy.py` | `phospy.core.references.manifest_policy` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/models.py` | `phospy.core.references.models` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/provenance/context.py` | `phospy.core.references.provenance.context` | `references` | Reference provenance record or direct codec |
| `src/phospy/core/references/provenance/identifiers.py` | `phospy.core.references.provenance.identifiers` | `references` | Reference provenance record or direct codec |
| `src/phospy/core/references/provenance/models.py` | `phospy.core.references.provenance.models` | `references` | Reference provenance record or direct codec |
| `src/phospy/core/references/provenance/serialization.py` | `phospy.core.references.provenance.serialization` | `references` | Reference provenance record or direct codec |
| `src/phospy/core/references/redistribution.py` | `phospy.core.references.redistribution` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/resolution.py` | `phospy.core.references.resolution` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/resources.py` | `phospy.core.references.resources` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/sequences/__init__.py` | `phospy.core.references.sequences` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/sequences/fasta.py` | `phospy.core.references.sequences.fasta` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/sequences/models.py` | `phospy.core.references.sequences.models` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/sequences/repository.py` | `phospy.core.references.sequences.repository` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/sequences/validation.py` | `phospy.core.references.sequences.validation` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/tables.py` | `phospy.core.references.tables` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/validation/__init__.py` | `phospy.core.references.validation` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/validation/_diagnostics.py` | `phospy.core.references.validation._diagnostics` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/validation/bundle.py` | `phospy.core.references.validation.bundle` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/validation/bundle_semantics.py` | `phospy.core.references.validation.bundle_semantics` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/validation/manifest_schema.py` | `phospy.core.references.validation.manifest_schema` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/validation/redistribution_policy.py` | `phospy.core.references.validation.redistribution_policy` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/references/validation/resource_integrity.py` | `phospy.core.references.validation.resource_integrity` | `references` | General reference identity, resources, sequences, or construction |
| `src/phospy/core/statistics/__init__.py` | `phospy.core.statistics` | `statistics` | Reusable statistical kernel or generic correction value |
| `src/phospy/core/statistics/multiple_testing.py` | `phospy.core.statistics.multiple_testing` | `statistics` | Reusable statistical kernel or generic correction value |
| `src/phospy/core/taxonomy/__init__.py` | `phospy.core.taxonomy` | `taxonomy` | Typed organism vocabulary and normalisation |
| `src/phospy/core/taxonomy/organisms.py` | `phospy.core.taxonomy.organisms` | `taxonomy` | Typed organism vocabulary and normalisation |
| `src/phospy/core/validation/__init__.py` | `phospy.core.validation` | `validation` | Reusable validation mechanics |
| `src/phospy/core/validation/common/__init__.py` | `phospy.core.validation.common` | `validation` | Reusable validation mechanics |
| `src/phospy/core/validation/common/config_values.py` | `phospy.core.validation.common.config_values` | `validation` | Reusable validation mechanics |
| `src/phospy/core/validation/common/input_sources.py` | `phospy.core.validation.common.input_sources` | `validation` | Policy-free DataFrame-or-path input validation mechanics |
| `src/phospy/core/validation/common/numbers.py` | `phospy.core.validation.common.numbers` | `validation` | Reusable validation mechanics |
| `src/phospy/core/validation/common/paths.py` | `phospy.core.validation.common.paths` | `validation` | Reusable validation mechanics |
| `src/phospy/core/validation/frame_columns.py` | `phospy.core.validation.frame_columns` | `validation` | Reusable validation mechanics |
| `src/phospy/core/validation/ownership.py` | `phospy.core.validation.ownership` | `validation` | Generic validation governance |
| `src/phospy/datasets/__init__.py` | `phospy.datasets` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_construction/__init__.py` | `phospy.datasets._construction` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_construction/fingerprints.py` | `phospy.datasets._construction.fingerprints` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_construction/organism_coherence.py` | `phospy.datasets._construction.organism_coherence` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_construction/protein_aware_binding.py` | `phospy.datasets._construction.protein_aware_binding` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_construction/service.py` | `phospy.datasets._construction.service` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_construction/trusted_assertions.py` | `phospy.datasets._construction.trusted_assertions` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_construction/validation.py` | `phospy.datasets._construction.validation` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_internal/__init__.py` | `phospy.datasets._internal` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_internal/frame_store.py` | `phospy.datasets._internal.frame_store` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_processing_state/__init__.py` | `phospy.datasets._processing_state` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_processing_state/json_contracts.py` | `phospy.datasets._processing_state.json_contracts` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_processing_state/missing_data.py` | `phospy.datasets._processing_state.missing_data` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_processing_state/models.py` | `phospy.datasets._processing_state.models` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/_processing_state/total_protein.py` | `phospy.datasets._processing_state.total_protein` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/builder.py` | `phospy.datasets.builder` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/__init__.py` | `phospy.datasets.construction` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/contracts.py` | `phospy.datasets.construction.contracts` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/executor.py` | `phospy.datasets.construction.executor` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/input_resolution.py` | `phospy.datasets.construction.input_resolution` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/intensity_scale.py` | `phospy.datasets.construction.intensity_scale` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/interpreter.py` | `phospy.datasets.construction.interpreter` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/normalization.py` | `phospy.datasets.construction.normalization` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/normalization_reporting.py` | `phospy.datasets.construction.normalization_reporting` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/planning.py` | `phospy.datasets.construction.planning` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/preprocessing.py` | `phospy.datasets.construction.preprocessing` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/protein_aware.py` | `phospy.datasets.construction.protein_aware` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/provenance.py` | `phospy.datasets.construction.provenance` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/report_assembly.py` | `phospy.datasets.construction.report_assembly` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/sample_metadata.py` | `phospy.datasets.construction.sample_metadata` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/sequence_derivation.py` | `phospy.datasets.construction.sequence_derivation` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/site_identity.py` | `phospy.datasets.construction.site_identity` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/site_metadata.py` | `phospy.datasets.construction.site_metadata` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/site_sequence_validation.py` | `phospy.datasets.construction.site_sequence_validation` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/state_assembly.py` | `phospy.datasets.construction.state_assembly` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/total_matrix.py` | `phospy.datasets.construction.total_matrix` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/construction/transformation_state.py` | `phospy.datasets.construction.transformation_state` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/contracts.py` | `phospy.datasets.contracts` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/derived_quantitative.py` | `phospy.datasets.derived_quantitative` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/direct_construction.py` | `phospy.datasets.direct_construction` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/imputation_metadata.py` | `phospy.datasets.imputation_metadata` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/internal_view.py` | `phospy.datasets.internal_view` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/models.py` | `phospy.datasets.models` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/processing_state.py` | `phospy.datasets.processing_state` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/protein_aware.py` | `phospy.datasets.protein_aware` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/schemas.py` | `phospy.datasets.schemas` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/trusted_assertions.py` | `phospy.datasets.trusted_assertions` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/trusted_construction_caveats.py` | `phospy.datasets.trusted_construction_caveats` | `datasets` | Dataset model construction invariant or durable attached record |
| `src/phospy/datasets/validation/__init__.py` | `phospy.datasets.validation` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/validation/inputs.py` | `phospy.datasets.validation.inputs` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/validation/model_boundary.py` | `phospy.datasets.validation.model_boundary` | `datasets` | Dataset construction lifecycle |
| `src/phospy/datasets/validation/request.py` | `phospy.datasets.validation.request` | `datasets` | Dataset construction lifecycle |
| `src/phospy/differential/__init__.py` | `phospy.differential` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/_payload_codec.py` | `phospy.differential._payload_codec` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/aggregation/__init__.py` | `phospy.differential.aggregation` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/aggregation/executor.py` | `phospy.differential.aggregation.executor` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/aggregation/experimental.py` | `phospy.differential.aggregation.experimental` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/aggregation/models.py` | `phospy.differential.aggregation.models` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/aggregation/scientific_policies.py` | `phospy.differential.aggregation.scientific_policies` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/caveats.py` | `phospy.differential.caveats` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/config.py` | `phospy.differential.config` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/contracts.py` | `phospy.differential.contracts` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/design/__init__.py` | `phospy.differential.design` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/design/contrast_helpers.py` | `phospy.differential.design.contrast_helpers` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/design/matrix_builder.py` | `phospy.differential.design.matrix_builder` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/design/models.py` | `phospy.differential.design.models` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/diagnostics.py` | `phospy.differential.diagnostics` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/execution/__init__.py` | `phospy.differential.execution` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/execution/compound_symmetry_gls.py` | `phospy.differential.execution.compound_symmetry_gls` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/execution/duplicate_correlation.py` | `phospy.differential.execution.duplicate_correlation` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/execution/fitting.py` | `phospy.differential.execution.fitting` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/execution/linear_model.py` | `phospy.differential.execution.linear_model` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/execution/model_executor.py` | `phospy.differential.execution.model_executor` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/execution/protein_covariate_adjusted.py` | `phospy.differential.execution.protein_covariate_adjusted` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/execution/workflow.py` | `phospy.differential.execution.workflow` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/internal_view.py` | `phospy.differential.internal_view` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/models/__init__.py` | `phospy.differential.models` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/models/design.py` | `phospy.differential.models.design` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/models/duplicate_correlation.py` | `phospy.differential.models.duplicate_correlation` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/models/fit.py` | `phospy.differential.models.fit` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/models/protein_aware.py` | `phospy.differential.models.protein_aware` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/moderation/__init__.py` | `phospy.differential.moderation` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/moderation/config.py` | `phospy.differential.moderation.config` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/moderation/empirical_bayes.py` | `phospy.differential.moderation.empirical_bayes` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/planning/__init__.py` | `phospy.differential.planning` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/planning/design.py` | `phospy.differential.planning.design` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/planning/eligibility.py` | `phospy.differential.planning.eligibility` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/planning/imputation.py` | `phospy.differential.planning.imputation` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/planning/interpreter.py` | `phospy.differential.planning.interpreter` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/planning/models.py` | `phospy.differential.planning.models` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/planning/protein_aware_eligibility.py` | `phospy.differential.planning.protein_aware_eligibility` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/planning/protein_aware_inputs.py` | `phospy.differential.planning.protein_aware_inputs` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/planning/replicates.py` | `phospy.differential.planning.replicates` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/policies.py` | `phospy.differential.policies` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/provenance/__init__.py` | `phospy.differential.provenance` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/provenance/assembly.py` | `phospy.differential.provenance.assembly` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/provenance/lineage.py` | `phospy.differential.provenance.lineage` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/provenance/models.py` | `phospy.differential.provenance.models` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/public.py` | `phospy.differential.public` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/quantification_depth.py` | `phospy.differential.quantification_depth` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/result_assembly.py` | `phospy.differential.result_assembly` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/results.py` | `phospy.differential.results` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/tables.py` | `phospy.differential.tables` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/validation/__init__.py` | `phospy.differential.validation` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/validation/computation.py` | `phospy.differential.validation.computation` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/validation/dataset.py` | `phospy.differential.validation.dataset` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/validation/design.py` | `phospy.differential.validation.design` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/validation/reliability.py` | `phospy.differential.validation.reliability` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/validation/result_tables.py` | `phospy.differential.validation.result_tables` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/differential/validation/workflow.py` | `phospy.differential.validation.workflow` | `differential` | Complete Differential capability lifecycle |
| `src/phospy/enrichment/__init__.py` | `phospy.enrichment` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/assembly.py` | `phospy.enrichment.assembly` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/caveats.py` | `phospy.enrichment.caveats` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/contracts/__init__.py` | `phospy.enrichment.contracts` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/contracts/collections.py` | `phospy.enrichment.contracts.collections` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/contracts/config.py` | `phospy.enrichment.contracts.config` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/contracts/identifier_sets.py` | `phospy.enrichment.contracts.identifier_sets` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/contracts/requests.py` | `phospy.enrichment.contracts.requests` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/contracts/results.py` | `phospy.enrichment.contracts.results` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/contracts/vocabulary.py` | `phospy.enrichment.contracts.vocabulary` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/execution.py` | `phospy.enrichment.execution` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/filtering.py` | `phospy.enrichment.filtering` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/interpreter.py` | `phospy.enrichment.interpreter` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/ora.py` | `phospy.enrichment.ora` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/ora_execution.py` | `phospy.enrichment.ora_execution` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/planning.py` | `phospy.enrichment.planning` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/provenance.py` | `phospy.enrichment.provenance` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/validation.py` | `phospy.enrichment.validation` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/enrichment/workflow.py` | `phospy.enrichment.workflow` | `enrichment` | Complete Enrichment capability lifecycle |
| `src/phospy/io/__init__.py` | `phospy.io` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/__init__.py` | `phospy.io.bundles` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_kinase/__init__.py` | `phospy.io.bundles._kinase` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_kinase/constants.py` | `phospy.io.bundles._kinase.constants` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_kinase/loader.py` | `phospy.io.bundles._kinase.loader` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_kinase/manifest.py` | `phospy.io.bundles._kinase.manifest` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_kinase/models.py` | `phospy.io.bundles._kinase.models` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_kinase/reconstruction.py` | `phospy.io.bundles._kinase.reconstruction` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_kinase/snapshots.py` | `phospy.io.bundles._kinase.snapshots` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_kinase/writer.py` | `phospy.io.bundles._kinase.writer` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_shared/__init__.py` | `phospy.io.bundles._shared` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_shared/integrity.py` | `phospy.io.bundles._shared.integrity` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_shared/json_files.py` | `phospy.io.bundles._shared.json_files` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_shared/kinase_scoring_state.py` | `phospy.io.bundles._shared.kinase_scoring_state` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_shared/organisms.py` | `phospy.io.bundles._shared.organisms` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_shared/paths.py` | `phospy.io.bundles._shared.paths` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_shared/primitives.py` | `phospy.io.bundles._shared.primitives` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_shared/processing_state.py` | `phospy.io.bundles._shared.processing_state` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_shared/tables.py` | `phospy.io.bundles._shared.tables` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_shared/transactions.py` | `phospy.io.bundles._shared.transactions` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_shared/trusted_dataset_assertions.py` | `phospy.io.bundles._shared.trusted_dataset_assertions` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_signalome/__init__.py` | `phospy.io.bundles._signalome` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_signalome/config.py` | `phospy.io.bundles._signalome.config` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_signalome/constants.py` | `phospy.io.bundles._signalome.constants` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_signalome/diagnostics.py` | `phospy.io.bundles._signalome.diagnostics` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_signalome/loader.py` | `phospy.io.bundles._signalome.loader` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_signalome/manifest.py` | `phospy.io.bundles._signalome.manifest` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_signalome/models.py` | `phospy.io.bundles._signalome.models` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_signalome/primitives.py` | `phospy.io.bundles._signalome.primitives` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_signalome/reconstruction.py` | `phospy.io.bundles._signalome.reconstruction` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_signalome/snapshots.py` | `phospy.io.bundles._signalome.snapshots` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_signalome/tables.py` | `phospy.io.bundles._signalome.tables` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/_signalome/writer.py` | `phospy.io.bundles._signalome.writer` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/kinase.py` | `phospy.io.bundles.kinase` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/kinase_library.py` | `phospy.io.bundles.kinase_library` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/bundles/signalome.py` | `phospy.io.bundles.signalome` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/importers/__init__.py` | `phospy.io.importers` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/importers/contracts.py` | `phospy.io.importers.contracts` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/importers/fragpipe_validation.py` | `phospy.io.importers.fragpipe_validation` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/importers/maxquant_validation.py` | `phospy.io.importers.maxquant_validation` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/importers/validation.py` | `phospy.io.importers.validation` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/publishers/__init__.py` | `phospy.io.publishers` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/publishers/workflows.py` | `phospy.io.publishers.workflows` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/__init__.py` | `phospy.io.readers` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/_table_parsing.py` | `phospy.io.readers._table_parsing` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/dataset_inputs.py` | `phospy.io.readers.dataset_inputs` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/enrichment_sets.py` | `phospy.io.readers.enrichment_sets` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/fragpipe/__init__.py` | `phospy.io.readers.fragpipe` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/fragpipe/columns.py` | `phospy.io.readers.fragpipe.columns` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/fragpipe/constants.py` | `phospy.io.readers.fragpipe.constants` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/fragpipe/conversion.py` | `phospy.io.readers.fragpipe.conversion` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/fragpipe/filtering.py` | `phospy.io.readers.fragpipe.filtering` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/fragpipe/importer.py` | `phospy.io.readers.fragpipe.importer` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/fragpipe/models.py` | `phospy.io.readers.fragpipe.models` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/fragpipe/normalization.py` | `phospy.io.readers.fragpipe.normalization` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/fragpipe/raw.py` | `phospy.io.readers.fragpipe.raw` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/fragpipe/reporting.py` | `phospy.io.readers.fragpipe.reporting` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/importers.py` | `phospy.io.readers.importers` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/maxquant.py` | `phospy.io.readers.maxquant` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/io/readers/tables.py` | `phospy.io.readers.tables` | `io` | Reader, persistence, or publisher adapter |
| `src/phospy/kinase/__init__.py` | `phospy.kinase` | `kinase` | Kinase capability package |
| `src/phospy/kinase/_matrix_validation.py` | `phospy.kinase._matrix_validation` | `kinase` | Kinase shared policy or mixed-owner correction |
| `src/phospy/kinase/activity/__init__.py` | `phospy.kinase.activity` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/config.py` | `phospy.kinase.activity.config` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/contracts.py` | `phospy.kinase.activity.contracts` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/diagnostics.py` | `phospy.kinase.activity.diagnostics` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/execution.py` | `phospy.kinase.activity.execution` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/inputs.py` | `phospy.kinase.activity.inputs` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/membership.py` | `phospy.kinase.activity.membership` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/membership_selection.py` | `phospy.kinase.activity.membership_selection` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/methods/__init__.py` | `phospy.kinase.activity.methods` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/methods/ksea_zscore.py` | `phospy.kinase.activity.methods.ksea_zscore` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/methods/ssgsea_substrate_enrichment.py` | `phospy.kinase.activity.methods.ssgsea_substrate_enrichment` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/methods/weighted_substrate_activity.py` | `phospy.kinase.activity.methods.weighted_substrate_activity` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/models.py` | `phospy.kinase.activity.models` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/policies.py` | `phospy.kinase.activity.policies` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/result_validation.py` | `phospy.kinase.activity.result_validation` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/results.py` | `phospy.kinase.activity.results` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/runner.py` | `phospy.kinase.activity.runner` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/stage_contracts.py` | `phospy.kinase.activity.stage_contracts` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/statistics.py` | `phospy.kinase.activity.statistics` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/tables.py` | `phospy.kinase.activity.tables` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/threshold_membership.py` | `phospy.kinase.activity.threshold_membership` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/activity/validation.py` | `phospy.kinase.activity.validation` | `kinase` | Kinase activity capability implementation |
| `src/phospy/kinase/attrition_metrics.py` | `phospy.kinase.attrition_metrics` | `kinase` | Kinase capability package |
| `src/phospy/kinase/caveats.py` | `phospy.kinase.caveats` | `kinase` | Kinase capability package |
| `src/phospy/kinase/contracts.py` | `phospy.kinase.contracts` | `kinase` | Kinase capability package |
| `src/phospy/kinase/eligibility.py` | `phospy.kinase.eligibility` | `kinase` | Kinase capability package |
| `src/phospy/kinase/executor.py` | `phospy.kinase.executor` | `kinase` | Kinase capability package |
| `src/phospy/kinase/interpreter.py` | `phospy.kinase.interpreter` | `kinase` | Kinase capability package |
| `src/phospy/kinase/policies.py` | `phospy.kinase.policies` | `kinase` | Kinase shared policy or mixed-owner correction |
| `src/phospy/kinase/prediction/__init__.py` | `phospy.kinase.prediction` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/prediction/candidates.py` | `phospy.kinase.prediction.candidates` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/prediction/config.py` | `phospy.kinase.prediction.config` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/prediction/diagnostics.py` | `phospy.kinase.prediction.diagnostics` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/prediction/execution.py` | `phospy.kinase.prediction.execution` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/prediction/internal_view.py` | `phospy.kinase.prediction.internal_view` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/prediction/outputs.py` | `phospy.kinase.prediction.outputs` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/prediction/policies.py` | `phospy.kinase.prediction.policies` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/prediction/results.py` | `phospy.kinase.prediction.results` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/prediction/runner.py` | `phospy.kinase.prediction.runner` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/prediction/sampling_core.py` | `phospy.kinase.prediction.sampling_core` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/prediction/sampling_runtime.py` | `phospy.kinase.prediction.sampling_runtime` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/prediction/svm.py` | `phospy.kinase.prediction.svm` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/prediction/tables.py` | `phospy.kinase.prediction.tables` | `kinase` | Kinase prediction contracts and scientific execution |
| `src/phospy/kinase/provenance.py` | `phospy.kinase.provenance` | `kinase` | Kinase capability package |
| `src/phospy/kinase/public.py` | `phospy.kinase.public` | `kinase` | Kinase capability package |
| `src/phospy/kinase/references/__init__.py` | `phospy.kinase.references` | `kinase` | Kinase reference projection and resource semantics |
| `src/phospy/kinase/references/loading.py` | `phospy.kinase.references.loading` | `kinase` | Kinase reference projection and resource semantics |
| `src/phospy/kinase/references/models.py` | `phospy.kinase.references.models` | `kinase` | Kinase reference projection and resource semantics |
| `src/phospy/kinase/references/projection.py` | `phospy.kinase.references.projection` | `kinase` | Kinase reference projection and resource semantics |
| `src/phospy/kinase/references/projection_contracts.py` | `phospy.kinase.references.projection_contracts` | `kinase` | Kinase reference projection and resource semantics |
| `src/phospy/kinase/references/provenance.py` | `phospy.kinase.references.provenance` | `kinase` | Kinase reference resource provenance |
| `src/phospy/kinase/references/validation.py` | `phospy.kinase.references.validation` | `kinase` | Kinase reference projection and resource semantics |
| `src/phospy/kinase/requests.py` | `phospy.kinase.requests` | `kinase` | Kinase capability package |
| `src/phospy/kinase/resolved_validator.py` | `phospy.kinase.resolved_validator` | `kinase` | Kinase capability package |
| `src/phospy/kinase/result_assembly.py` | `phospy.kinase.result_assembly` | `kinase` | Kinase capability package |
| `src/phospy/kinase/result_boundary.py` | `phospy.kinase.result_boundary` | `kinase` | Kinase capability package |
| `src/phospy/kinase/results.py` | `phospy.kinase.results` | `kinase` | Kinase capability package |
| `src/phospy/kinase/row_attrition.py` | `phospy.kinase.row_attrition` | `kinase` | Kinase capability package |
| `src/phospy/kinase/scoring/__init__.py` | `phospy.kinase.scoring` | `kinase` | Kinase scoring scientific domain |
| `src/phospy/kinase/scoring/config.py` | `phospy.kinase.scoring.config` | `kinase` | Kinase scoring scientific domain |
| `src/phospy/kinase/scoring/contracts.py` | `phospy.kinase.scoring.contracts` | `kinase` | Kinase scoring scientific domain |
| `src/phospy/kinase/scoring/contributions.py` | `phospy.kinase.scoring.contributions` | `kinase` | Kinase scoring scientific domain |
| `src/phospy/kinase/scoring/execution.py` | `phospy.kinase.scoring.execution` | `kinase` | Kinase scoring scientific domain |
| `src/phospy/kinase/scoring/kinase_library.py` | `phospy.kinase.scoring.kinase_library` | `kinase` | Kinase scoring scientific domain |
| `src/phospy/kinase/scoring/motif/__init__.py` | `phospy.kinase.scoring.motif` | `kinase` | Kinase motif-scoring scientific domain |
| `src/phospy/kinase/scoring/motif/frequency_matrices.py` | `phospy.kinase.scoring.motif.frequency_matrices` | `kinase` | Kinase motif-scoring scientific domain |
| `src/phospy/kinase/scoring/motif/kinase_library.py` | `phospy.kinase.scoring.motif.kinase_library` | `kinase` | Kinase motif-scoring scientific domain |
| `src/phospy/kinase/scoring/motif/library_validation.py` | `phospy.kinase.scoring.motif.library_validation` | `kinase` | Kinase motif-scoring scientific domain |
| `src/phospy/kinase/scoring/motif/models.py` | `phospy.kinase.scoring.motif.models` | `kinase` | Kinase motif-scoring scientific domain |
| `src/phospy/kinase/scoring/motif/scaling.py` | `phospy.kinase.scoring.motif.scaling` | `kinase` | Kinase motif-scoring scientific domain |
| `src/phospy/kinase/scoring/motif/scoring.py` | `phospy.kinase.scoring.motif.scoring` | `kinase` | Kinase motif-scoring scientific domain |
| `src/phospy/kinase/scoring/policies.py` | `phospy.kinase.scoring.policies` | `kinase` | Kinase scoring scientific domain |
| `src/phospy/kinase/scoring/profile.py` | `phospy.kinase.scoring.profile` | `kinase` | Kinase scoring scientific domain |
| `src/phospy/kinase/scoring/rank_fusion.py` | `phospy.kinase.scoring.rank_fusion` | `kinase` | Kinase scoring scientific domain |
| `src/phospy/kinase/scoring/results.py` | `phospy.kinase.scoring.results` | `kinase` | Kinase scoring scientific domain |
| `src/phospy/kinase/scoring/selection.py` | `phospy.kinase.scoring.selection` | `kinase` | Kinase scoring scientific domain |
| `src/phospy/kinase/scoring/sequence_validation.py` | `phospy.kinase.scoring.sequence_validation` | `kinase` | Kinase shared policy or mixed-owner correction |
| `src/phospy/kinase/scoring/tables.py` | `phospy.kinase.scoring.tables` | `kinase` | Kinase scoring scientific domain |
| `src/phospy/kinase/scoring/transforms.py` | `phospy.kinase.scoring.transforms` | `kinase` | Kinase scoring scientific domain |
| `src/phospy/kinase/sequence_contracts.py` | `phospy.kinase.sequence_contracts` | `kinase` | Kinase capability package |
| `src/phospy/kinase/site_attrition.py` | `phospy.kinase.site_attrition` | `kinase` | Kinase capability package |
| `src/phospy/kinase/site_sequence_policy.py` | `phospy.kinase.site_sequence_policy` | `kinase` | Kinase capability package |
| `src/phospy/kinase/site_sequence_support.py` | `phospy.kinase.site_sequence_support` | `kinase` | Kinase capability package |
| `src/phospy/kinase/validator.py` | `phospy.kinase.validator` | `kinase` | Kinase capability package |
| `src/phospy/preparation_records/__init__.py` | `phospy.preparation_records` | `preparation_records` | Passive preparation record or structural validation |
| `src/phospy/preparation_records/batch_correction.py` | `phospy.preparation_records.batch_correction` | `preparation_records` | Passive preparation record or structural validation |
| `src/phospy/preparation_records/group_evidence.py` | `phospy.preparation_records.group_evidence` | `preparation_records` | Passive preparation record or structural validation |
| `src/phospy/preparation_records/imputation.py` | `phospy.preparation_records.imputation` | `preparation_records` | Passive preparation record or structural validation |
| `src/phospy/preparation_records/policies.py` | `phospy.preparation_records.policies` | `preparation_records` | Passive preparation record or structural validation |
| `src/phospy/preparation_records/protein_aware.py` | `phospy.preparation_records.protein_aware` | `preparation_records` | Passive preparation record or structural validation |
| `src/phospy/preparation_records/report_schema.py` | `phospy.preparation_records.report_schema` | `preparation_records` | Passive preparation record or structural validation |
| `src/phospy/preparation_records/validation/__init__.py` | `phospy.preparation_records.validation` | `preparation_records` | Passive preparation record or structural validation |
| `src/phospy/preparation_records/validation/sample_groups.py` | `phospy.preparation_records.validation.sample_groups` | `preparation_records` | Passive preparation record or structural validation |
| `src/phospy/preprocessing/__init__.py` | `phospy.preprocessing` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/batch_correction/__init__.py` | `phospy.preprocessing.batch_correction` | `preprocessing` | Nested preprocessing batch-correction capability |
| `src/phospy/preprocessing/batch_correction/adapters/__init__.py` | `phospy.preprocessing.batch_correction.adapters` | `preprocessing` | Batch-correction preprocessing adaptation |
| `src/phospy/preprocessing/batch_correction/adapters/preprocessing.py` | `phospy.preprocessing.batch_correction.adapters.preprocessing` | `preprocessing` | Batch-correction preprocessing adaptation |
| `src/phospy/preprocessing/batch_correction/contracts/__init__.py` | `phospy.preprocessing.batch_correction.contracts` | `preprocessing` | Passive batch-correction capability contract |
| `src/phospy/preprocessing/batch_correction/contracts/_validation.py` | `phospy.preprocessing.batch_correction.contracts._validation` | `preprocessing` | Passive batch-correction capability contract |
| `src/phospy/preprocessing/batch_correction/contracts/configuration.py` | `phospy.preprocessing.batch_correction.contracts.configuration` | `preprocessing` | Passive batch-correction capability contract |
| `src/phospy/preprocessing/batch_correction/contracts/missingness.py` | `phospy.preprocessing.batch_correction.contracts.missingness` | `preprocessing` | Passive batch-correction capability contract |
| `src/phospy/preprocessing/batch_correction/contracts/output.py` | `phospy.preprocessing.batch_correction.contracts.output` | `preprocessing` | Passive batch-correction capability contract |
| `src/phospy/preprocessing/batch_correction/contracts/request.py` | `phospy.preprocessing.batch_correction.contracts.request` | `preprocessing` | Passive batch-correction capability contract |
| `src/phospy/preprocessing/batch_correction/contracts/results.py` | `phospy.preprocessing.batch_correction.contracts.results` | `preprocessing` | Passive batch-correction capability contract |
| `src/phospy/preprocessing/batch_correction/contracts/stage.py` | `phospy.preprocessing.batch_correction.contracts.stage` | `preprocessing` | Passive batch-correction capability contract |
| `src/phospy/preprocessing/batch_correction/contracts/workflow.py` | `phospy.preprocessing.batch_correction.contracts.workflow` | `preprocessing` | Passive batch-correction capability contract |
| `src/phospy/preprocessing/batch_correction/environment.py` | `phospy.preprocessing.batch_correction.environment` | `preprocessing` | Batch-correction environment dependency selection |
| `src/phospy/preprocessing/batch_correction/execution/__init__.py` | `phospy.preprocessing.batch_correction.execution` | `preprocessing` | Batch-correction numerical execution |
| `src/phospy/preprocessing/batch_correction/execution/linear.py` | `phospy.preprocessing.batch_correction.execution.linear` | `preprocessing` | Batch-correction numerical execution |
| `src/phospy/preprocessing/batch_correction/execution/missingness.py` | `phospy.preprocessing.batch_correction.execution.missingness` | `preprocessing` | Batch-correction numerical execution |
| `src/phospy/preprocessing/batch_correction/execution/ruv_iii.py` | `phospy.preprocessing.batch_correction.execution.ruv_iii` | `preprocessing` | Batch-correction numerical execution |
| `src/phospy/preprocessing/batch_correction/execution/ruv_iii_executor.py` | `phospy.preprocessing.batch_correction.execution.ruv_iii_executor` | `preprocessing` | Batch-correction numerical execution |
| `src/phospy/preprocessing/batch_correction/execution/sps_ruv.py` | `phospy.preprocessing.batch_correction.execution.sps_ruv` | `preprocessing` | Batch-correction numerical execution |
| `src/phospy/preprocessing/batch_correction/execution/stage.py` | `phospy.preprocessing.batch_correction.execution.stage` | `preprocessing` | Batch-correction numerical execution |
| `src/phospy/preprocessing/batch_correction/metadata.py` | `phospy.preprocessing.batch_correction.metadata` | `preprocessing` | Batch-correction metadata resolution |
| `src/phospy/preprocessing/batch_correction/planning.py` | `phospy.preprocessing.batch_correction.planning` | `preprocessing` | Batch-correction interpretation and planning |
| `src/phospy/preprocessing/batch_correction/provenance/__init__.py` | `phospy.preprocessing.batch_correction.provenance` | `preprocessing` | Batch-correction provenance assembly and scientific validation |
| `src/phospy/preprocessing/batch_correction/provenance/_payloads.py` | `phospy.preprocessing.batch_correction.provenance._payloads` | `preprocessing` | Batch-correction provenance assembly and scientific validation |
| `src/phospy/preprocessing/batch_correction/provenance/assembly.py` | `phospy.preprocessing.batch_correction.provenance.assembly` | `preprocessing` | Batch-correction provenance assembly and scientific validation |
| `src/phospy/preprocessing/batch_correction/provenance/validation.py` | `phospy.preprocessing.batch_correction.provenance.validation` | `preprocessing` | Batch-correction provenance assembly and scientific validation |
| `src/phospy/preprocessing/batch_correction/validation/__init__.py` | `phospy.preprocessing.batch_correction.validation` | `preprocessing` | Batch-correction scientific and request validation |
| `src/phospy/preprocessing/batch_correction/validation/_helpers.py` | `phospy.preprocessing.batch_correction.validation._helpers` | `preprocessing` | Batch-correction scientific and request validation |
| `src/phospy/preprocessing/batch_correction/validation/controls.py` | `phospy.preprocessing.batch_correction.validation.controls` | `preprocessing` | Batch-correction scientific and request validation |
| `src/phospy/preprocessing/batch_correction/validation/design.py` | `phospy.preprocessing.batch_correction.validation.design` | `preprocessing` | Batch-correction scientific and request validation |
| `src/phospy/preprocessing/batch_correction/validation/feasibility.py` | `phospy.preprocessing.batch_correction.validation.feasibility` | `preprocessing` | Batch-correction scientific and request validation |
| `src/phospy/preprocessing/batch_correction/validation/missingness.py` | `phospy.preprocessing.batch_correction.validation.missingness` | `preprocessing` | Batch-correction scientific and request validation |
| `src/phospy/preprocessing/batch_correction/validation/output.py` | `phospy.preprocessing.batch_correction.validation.output` | `preprocessing` | Batch-correction scientific and request validation |
| `src/phospy/preprocessing/batch_correction/validation/protocols.py` | `phospy.preprocessing.batch_correction.validation.protocols` | `preprocessing` | Batch-correction scientific and request validation |
| `src/phospy/preprocessing/batch_correction/validation/request.py` | `phospy.preprocessing.batch_correction.validation.request` | `preprocessing` | Batch-correction scientific and request validation |
| `src/phospy/preprocessing/batch_correction/validation/stage_order.py` | `phospy.preprocessing.batch_correction.validation.stage_order` | `preprocessing` | Batch-correction scientific and request validation |
| `src/phospy/preprocessing/batch_correction/workflow.py` | `phospy.preprocessing.batch_correction.workflow` | `preprocessing` | Batch-correction workflow composition |
| `src/phospy/preprocessing/contracts/__init__.py` | `phospy.preprocessing.contracts` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/configuration/__init__.py` | `phospy.preprocessing.contracts.configuration` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/configuration/comparisons.py` | `phospy.preprocessing.contracts.configuration.comparisons` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/configuration/coverage_filter.py` | `phospy.preprocessing.contracts.configuration.coverage_filter` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/configuration/dataset.py` | `phospy.preprocessing.contracts.configuration.dataset` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/configuration/intensity_transform.py` | `phospy.preprocessing.contracts.configuration.intensity_transform` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/configuration/localisation.py` | `phospy.preprocessing.contracts.configuration.localisation` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/configuration/missing_data.py` | `phospy.preprocessing.contracts.configuration.missing_data` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/configuration/normalisation.py` | `phospy.preprocessing.contracts.configuration.normalisation` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/configuration/site_matrix.py` | `phospy.preprocessing.contracts.configuration.site_matrix` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/configuration/site_sequence.py` | `phospy.preprocessing.contracts.configuration.site_sequence` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/configuration/total_protein.py` | `phospy.preprocessing.contracts.configuration.total_protein` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/imputation_scale.py` | `phospy.preprocessing.contracts.imputation_scale` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/plan.py` | `phospy.preprocessing.contracts.plan` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/plan_resolution.py` | `phospy.preprocessing.contracts.plan_resolution` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/policies.py` | `phospy.preprocessing.contracts.policies` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/stage.py` | `phospy.preprocessing.contracts.stage` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/stage_keys.py` | `phospy.preprocessing.contracts.stage_keys` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/contracts/state.py` | `phospy.preprocessing.contracts.state` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/control_sites/__init__.py` | `phospy.preprocessing.control_sites` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/control_sites/_discovery_payload.py` | `phospy.preprocessing.control_sites._discovery_payload` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/control_sites/contracts.py` | `phospy.preprocessing.control_sites.contracts` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/control_sites/discovery_contracts.py` | `phospy.preprocessing.control_sites.discovery_contracts` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/control_sites/discovery_execution.py` | `phospy.preprocessing.control_sites.discovery_execution` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/control_sites/discovery_planning.py` | `phospy.preprocessing.control_sites.discovery_planning` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/control_sites/discovery_provenance.py` | `phospy.preprocessing.control_sites.discovery_provenance` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/control_sites/discovery_ranking.py` | `phospy.preprocessing.control_sites.discovery_ranking` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/control_sites/discovery_results.py` | `phospy.preprocessing.control_sites.discovery_results` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/control_sites/discovery_validation.py` | `phospy.preprocessing.control_sites.discovery_validation` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/control_sites/validation.py` | `phospy.preprocessing.control_sites.validation` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/control_sites/workflow.py` | `phospy.preprocessing.control_sites.workflow` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/dataset_construction.py` | `phospy.preprocessing.dataset_construction` | `preprocessing` | Published Preprocessing composition boundary for Dataset construction |
| `src/phospy/preprocessing/evidence/missing_data.py` | `phospy.preprocessing.evidence.missing_data` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/execution/corrected_output.py` | `phospy.preprocessing.execution.corrected_output` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/execution/diagnostics_normalization.py` | `phospy.preprocessing.execution.diagnostics_normalization` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/execution/fingerprints.py` | `phospy.preprocessing.execution.fingerprints` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/execution/pipeline.py` | `phospy.preprocessing.execution.pipeline` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/execution/trace_builder.py` | `phospy.preprocessing.execution.trace_builder` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/planning/assembly.py` | `phospy.preprocessing.planning.assembly` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/planning/boundary.py` | `phospy.preprocessing.planning.boundary` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/planning/configuration.py` | `phospy.preprocessing.planning.configuration` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/planning/ordering.py` | `phospy.preprocessing.planning.ordering` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/planning/registry.py` | `phospy.preprocessing.planning.registry` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/planning/rules.py` | `phospy.preprocessing.planning.rules` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/planning/validation.py` | `phospy.preprocessing.planning.validation` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/protein_aware/__init__.py` | `phospy.preprocessing.protein_aware` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/protein_aware/alignment.py` | `phospy.preprocessing.protein_aware.alignment` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/protein_aware/identity.py` | `phospy.preprocessing.protein_aware.identity` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/protein_aware/mapping.py` | `phospy.preprocessing.protein_aware.mapping` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/provenance/composition.py` | `phospy.preprocessing.provenance.composition` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/provenance/models.py` | `phospy.preprocessing.provenance.models` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/provenance/policies.py` | `phospy.preprocessing.provenance.policies` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/provenance/serialization.py` | `phospy.preprocessing.provenance.serialization` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/results/diagnostics.py` | `phospy.preprocessing.results.diagnostics` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/results/execution.py` | `phospy.preprocessing.results.execution` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/results/missing_data.py` | `phospy.preprocessing.results.missing_data` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/results/protein_aware.py` | `phospy.preprocessing.results.protein_aware` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/results/quantitative_evidence.py` | `phospy.preprocessing.results.quantitative_evidence` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/results/report_rows.py` | `phospy.preprocessing.results.report_rows` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/results/report_schema.py` | `phospy.preprocessing.results.report_schema` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/sequence/__init__.py` | `phospy.preprocessing.sequence` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/sequence/conflicts.py` | `phospy.preprocessing.sequence.conflicts` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/sequence/diagnostics.py` | `phospy.preprocessing.sequence.diagnostics` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/sequence/metadata.py` | `phospy.preprocessing.sequence.metadata` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/sequence/models.py` | `phospy.preprocessing.sequence.models` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/sequence/reference.py` | `phospy.preprocessing.sequence.reference` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/sequence/request.py` | `phospy.preprocessing.sequence.request` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/__init__.py` | `phospy.preprocessing.stages` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/batch_correction.py` | `phospy.preprocessing.stages.batch_correction` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/comparisons.py` | `phospy.preprocessing.stages.comparisons` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/group_coverage.py` | `phospy.preprocessing.stages.group_coverage` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/intensity.py` | `phospy.preprocessing.stages.intensity` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/localisation.py` | `phospy.preprocessing.stages.localisation` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/missing_data/__init__.py` | `phospy.preprocessing.stages.missing_data` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/missing_data/audit.py` | `phospy.preprocessing.stages.missing_data.audit` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/missing_data/diagnostics.py` | `phospy.preprocessing.stages.missing_data.diagnostics` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/missing_data/forbid.py` | `phospy.preprocessing.stages.missing_data.forbid` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/missing_data/group_aware_routing.py` | `phospy.preprocessing.stages.missing_data.group_aware_routing` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/missing_data/knn.py` | `phospy.preprocessing.stages.missing_data.knn` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/missing_data/minprob.py` | `phospy.preprocessing.stages.missing_data.minprob` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/missing_data/models.py` | `phospy.preprocessing.stages.missing_data.models` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/missing_data/row_median.py` | `phospy.preprocessing.stages.missing_data.row_median` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/missing_data/stage.py` | `phospy.preprocessing.stages.missing_data.stage` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/normalisation.py` | `phospy.preprocessing.stages.normalisation` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/protein_aware.py` | `phospy.preprocessing.stages.protein_aware` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/sequence.py` | `phospy.preprocessing.stages.sequence` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/site_matrix.py` | `phospy.preprocessing.stages.site_matrix` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/site_matrix_components/__init__.py` | `phospy.preprocessing.stages.site_matrix_components` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/site_matrix_components/audit.py` | `phospy.preprocessing.stages.site_matrix_components.audit` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/site_matrix_components/duplicates.py` | `phospy.preprocessing.stages.site_matrix_components.duplicates` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/site_matrix_components/identity.py` | `phospy.preprocessing.stages.site_matrix_components.identity` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/site_matrix_components/metadata.py` | `phospy.preprocessing.stages.site_matrix_components.metadata` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/stages/total_protein.py` | `phospy.preprocessing.stages.total_protein` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/validation/configuration.py` | `phospy.preprocessing.validation.configuration` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/validation/configuration_fields.py` | `phospy.preprocessing.validation.configuration_fields` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/validation/group_coverage.py` | `phospy.preprocessing.validation.group_coverage` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/validation/quantitative_evidence.py` | `phospy.preprocessing.validation.quantitative_evidence` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/validation/quantitative_transition.py` | `phospy.preprocessing.validation.quantitative_transition` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/preprocessing/workflow.py` | `phospy.preprocessing.workflow` | `preprocessing` | Preprocessing control-site capability lifecycle |
| `src/phospy/run_lineage/__init__.py` | `phospy.run_lineage` | `run_lineage` | Cross-capability run-lineage composition owner |
| `src/phospy/run_lineage/models.py` | `phospy.run_lineage.models` | `run_lineage` | Cross-capability run-lineage composition owner |
| `src/phospy/run_lineage/serialization.py` | `phospy.run_lineage.serialization` | `run_lineage` | Cross-capability run-lineage composition owner |
| `src/phospy/signalome/__init__.py` | `phospy.signalome` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/alignment_diagnostics.py` | `phospy.signalome.alignment_diagnostics` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/assignments.py` | `phospy.signalome.assignments` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/boundary_errors.py` | `phospy.signalome.boundary_errors` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/caveats.py` | `phospy.signalome.caveats` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/__init__.py` | `phospy.signalome.clustering` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/backend_dispatch.py` | `phospy.signalome.clustering.backend_dispatch` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/backends/__init__.py` | `phospy.signalome.clustering.backends` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/backends/exact_python.py` | `phospy.signalome.clustering.backends.exact_python` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/backends/scipy_hierarchical.py` | `phospy.signalome.clustering.backends.scipy_hierarchical` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/candidate_scoring.py` | `phospy.signalome.clustering.candidate_scoring` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/candidate_selection.py` | `phospy.signalome.clustering.candidate_selection` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/contracts.py` | `phospy.signalome.clustering.contracts` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/diagnostic_schemas.py` | `phospy.signalome.clustering.diagnostic_schemas` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/diagnostics.py` | `phospy.signalome.clustering.diagnostics` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/exact_python.py` | `phospy.signalome.clustering.exact_python` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/models.py` | `phospy.signalome.clustering.models` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/module_selection.py` | `phospy.signalome.clustering.module_selection` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/orchestration.py` | `phospy.signalome.clustering.orchestration` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/policies.py` | `phospy.signalome.clustering.policies` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/protein_modules.py` | `phospy.signalome.clustering.protein_modules` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/scale_guards.py` | `phospy.signalome.clustering.scale_guards` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/scientific_policies.py` | `phospy.signalome.clustering.scientific_policies` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/scipy_hierarchical.py` | `phospy.signalome.clustering.scipy_hierarchical` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/scoring.py` | `phospy.signalome.clustering.scoring` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/selection.py` | `phospy.signalome.clustering.selection` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/stability.py` | `phospy.signalome.clustering.stability` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/tree_building.py` | `phospy.signalome.clustering.tree_building` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/tree_engine_adapter.py` | `phospy.signalome.clustering.tree_engine_adapter` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering/validation.py` | `phospy.signalome.clustering.validation` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/clustering_runner.py` | `phospy.signalome.clustering_runner` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/component_helpers.py` | `phospy.signalome.component_helpers` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/component_models.py` | `phospy.signalome.component_models` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/config.py` | `phospy.signalome.config` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/constants.py` | `phospy.signalome.constants` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/context.py` | `phospy.signalome.context` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/context_tables.py` | `phospy.signalome.context_tables` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/contracts.py` | `phospy.signalome.contracts` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/errors.py` | `phospy.signalome.errors` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/executor.py` | `phospy.signalome.executor` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/expanded.py` | `phospy.signalome.expanded` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/interpreter.py` | `phospy.signalome.interpreter` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/matrix_alignment.py` | `phospy.signalome.matrix_alignment` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/models.py` | `phospy.signalome.models` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/module_tables.py` | `phospy.signalome.module_tables` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/modules.py` | `phospy.signalome.modules` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/network.py` | `phospy.signalome.network` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/network_builder.py` | `phospy.signalome.network_builder` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/protein_resolution.py` | `phospy.signalome.protein_resolution` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/provenance.py` | `phospy.signalome.provenance` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/public.py` | `phospy.signalome.public` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/requests.py` | `phospy.signalome.requests` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/result_assembly.py` | `phospy.signalome.result_assembly` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/result_validation.py` | `phospy.signalome.result_validation` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/results.py` | `phospy.signalome.results` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/row_attrition.py` | `phospy.signalome.row_attrition` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/scientific_policies.py` | `phospy.signalome.scientific_policies` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/score_matrix_selection.py` | `phospy.signalome.score_matrix_selection` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/score_preconditioning.py` | `phospy.signalome.score_preconditioning` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/tables/__init__.py` | `phospy.signalome.tables` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/tables/assignments.py` | `phospy.signalome.tables.assignments` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/tables/common.py` | `phospy.signalome.tables.common` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/tables/context.py` | `phospy.signalome.tables.context` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/tables/modules.py` | `phospy.signalome.tables.modules` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/tables/network.py` | `phospy.signalome.tables.network` | `signalome` | Complete Signalome capability lifecycle |
| `src/phospy/signalome/validator.py` | `phospy.signalome.validator` | `signalome` | Complete Signalome capability lifecycle |

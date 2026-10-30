# Architecture

PhosPy uses a domain-first package layout. Capability roots own application
behaviour: `datasets`, `preprocessing`, `differential`, `enrichment`, `kinase`,
and `signalome`. Shared scientific roots own reusable meaning and data contracts:
`activities`, `evidence`, `phosphosites`, `quantitative`, `references`,
`statistics`, `taxonomy`, and `provenance`. `preparation_records` owns passive
preparation records, while `run_lineage` composes provenance across capabilities.
`errors`, `frames`, `policies`, and `validation` are technical foundations.

Dependencies point from capabilities toward shared foundations and published
contracts. The `io` root is an inward adapter and can depend only on reviewed
contract, model, provenance, and serialization paths. Architecture reporting
distinguishes two owner-graph views:

- The **complete statically observed owner graph (before classifications)**
  includes every owner edge found by the maintained import analyser. It is
  checked against an explicit exact-rule table. The current graph contains the
  non-trivial strongly connected component `datasets`, `io`, `kinase`, and
  `signalome`; in particular, `datasets -> io` is the accepted builder-to-path
  reader dependency and `io -> datasets` covers reconstruction/model use.
- The **policy-filtered owner graph used by the DAG gate** removes only the
  three exact imports recorded in
  `architecture/permanent-condensation-classifications.json` before owner
  condensation. The current policy-filtered graph is acyclic. Each excluded
  import remains reported with its classification identifier and rationale.

An accepted classification is an auditable dependency and is not reported as
an unauthorised edge. The classifications do not create broader owner or
package permissions, and temporary architecture exceptions are not supported.
The current non-policy
[dataset path ingestion ownership assessment](adr/assessment_arch_002b_dataset_path_ingestion_ownership.md)
recommends retaining the exact builder-to-reader classification; it records the
responsibility and behaviour evidence without authorising a migration.

The analyser is static. It includes imports below `TYPE_CHECKING` and function
bodies, plus literal `importlib.import_module` and `__import__` calls. It cannot
resolve computed module names, assignment/re-export bindings, plugin discovery,
import hooks, or other runtime-generated imports. The separately enforced eager
module graph covers imports reachable during module and class-body execution.

Capability imports from shared scientific and technical-foundation owners are
also checked at the resolved target. A target module component or imported
symbol beginning with `_` is rejected unless an exact reviewed
`(source module, target module, target symbol)` permission exists. Import aliases
do not change that resolved identity. ARCH-016 retains no such private-import
permission; the three ADR-0053 package-condensation classifications remain a
separate policy and are not private-import exemptions.

The reviewed cross-owner seams now use these owner-internal contracts. None of
the superseded imports had an exact authorization; they passed only through the
former broad capability-to-shared/foundation rule.

| Superseded private seam | Owner-internal contract | Consumers | Purpose |
| --- | --- | --- | --- |
| `phospy.quantitative._authority` | `phospy.quantitative.authority` | Dataset intensity-scale/state construction and preprocessing SPS discovery | Opaque scale-establishment and quantitative-meaning authority tokens; historical serialized `source` values remain unchanged. |
| `phospy.quantitative.transformers._identifiers` | `phospy.quantitative.transformers.identifiers` | Dataset intensity-scale resolution | Stable built-in/custom transformer provenance identifiers. |
| `_run_declared_scale_sanity_diagnostics` | `phospy.quantitative.declared_scale_diagnostics.run_declared_scale_sanity_diagnostics` | Dataset intensity-scale resolution and SPS reference declarations | Shared declared-scale sanity diagnostics. |
| `phospy.provenance.models._shared` and its private scalar guards | `phospy.provenance.scalar_validation` | Dataset assertions, preprocessing provenance/batch contracts, and Kinase reference provenance | JSON value typing and field-qualified scalar guards. |
| `phospy.provenance.serialization._payload` | `phospy.provenance.serialization.payload` | Preprocessing provenance serialization | Historical provenance payload validation/coercion. |
| `phospy.frames.ownership._borrow_dataframe` and `_borrow_optional_dataframe` | Existing `borrow_dataframe` and `borrow_optional_dataframe` owner contract | Kinase scoring and prediction results | Immutable internal dataframe snapshots without public mutable borrowing. |

These modules are not exported by `phospy`, `phospy.api`, or
`phospy.advanced`; direct owner-module imports remain unsupported for external
callers.

The supported public surfaces are the `phospy`, `phospy.api`, and
`phospy.advanced` facades. Public exports preserve the exact object identities of
their domain owners. Domain code cannot depend on those aggregate facades, and
removed compatibility routes are tested as negative imports.

The machine-readable authority lives in the top-level `architecture/` directory:

- `module-owners.csv` assigns every Python source module to exactly one owner;
- `public-api-baseline.json` fixes the supported public inventory and identities;
- `final-state-classifications.json` classifies permanent lazy loading, safety
  seams, and versioned legacy-data behaviour;
- `permanent-condensation-classifications.json` records the exact exceptional
  package edges that are part of the permanent architecture.

The generated [architecture inventory](architecture_inventory.md) lists every
tracked production file, module owner, top-level root, and permitted owner edge.
Its source is the policy code and machine-readable inputs above. The governing
decisions are [ADR-0053](adr/adr_0053_domain_first_capability_architecture.md) and
[ADR-0054](adr/adr_0054_public_api_reset_and_compatibility_removal.md).

[ADR-0055](adr/adr_0055_nested_owner_namespaces_and_api_consolidation.md)
approves the next directory-only refinement: the reusable owners will become
separate children of the inert `phospy.core` namespace, the specialist facade
will become `phospy.api.advanced`, the private inventory will move beneath
`phospy.api`, and the three root conveniences will be retired while remaining
available through `phospy.api`. Until the corresponding migration tickets move
those files atomically, the live routes and the maintained API baseline above
remain authoritative.

Architecture enforcement resolves owner-pure package units from the owner map,
so nested siblings such as future `phospy.core.frames` and
`phospy.core.quantitative` remain distinct. Physical top-level roots are
reported separately and do not define owner permissions or owner-level cycles.

Reusable graph analysis, ownership, final-policy, and inventory rendering code
is maintained in `tools/architecture`. Architecture tests import that maintainer
implementation; the runtime `phospy` package does not depend on it. From the
repository root, generate or check the inventory with:

```console
python -m tools.architecture.generate_architecture_inventory --repo-root .
python -m tools.architecture.generate_architecture_inventory --repo-root . --check
```

Direct execution of `tools/architecture/generate_architecture_inventory.py`
remains supported for existing maintainer workflows.

## Migration Status

The domain-first architecture migration, including its final Wave 10 gate, was
accepted on 19 October 2026 at commit
`a4b9929fb41f6507cab3011dafc1012e44c5ce19`, tree
`ecaba2700d978f31e170f02efe38d2a83f11ee6d`, with sole parent
`f2c00c9290187dd6b73a97114240808eb0f3abc0`. The independent review returned
`INDEPENDENT_REVIEW_PASS`, and authoritative hosted run
[`39161306082`](https://github.com/falconsmilie/phospy/actions/runs/39161306082)
succeeded in all 43 jobs. The accepted qualifications and complete hosted
matrix are recorded in [Release validation](release_validation.md).

This acceptance completes the architecture migration. It does not create or
authorise a tag, package publication, deployment, or release publication.

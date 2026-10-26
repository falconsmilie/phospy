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

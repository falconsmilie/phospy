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
contract, model, provenance, and serialization paths. The exact observed owner
graph is acyclic and is checked against an explicit allow-list. Three narrowly
scoped permanent condensation edges are recorded in
`architecture/permanent-condensation-classifications.json`.
Temporary architecture exceptions are not supported.

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

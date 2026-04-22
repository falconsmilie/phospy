# Troubleshooting

Use this page when a supported workflow fails and you need a fast recovery path.

## First Triage

1. Confirm you are following the supported chain:
   dataset build -> kinase -> signalome (optional).
2. Confirm required inputs exist:
   `phospho`, `site_metadata`, and `organism` when using `ReferencePreset.AUTO`.
3. Match the error to one of the sections below.

## Dataset Build and Input-Shape Issues

Common causes:

- missing required `site_metadata` fields (`gene_symbol`, `site`)
- non-canonical or mismatched site IDs across `phospho.index` and metadata
- unsupported input source or file type

Next steps:

- review [Quickstart](getting-started/quickstart-first-workflow.md) input examples
- verify strict rules in [Validation Guide](validation.md)
- confirm builder contract details in [API Guide](api.md)

## Reference and Organism Resolution Issues

Common causes:

- `ReferencePreset.AUTO` used without `dataset.organism`
- expecting bundled human/mouse references (bundled runtime is rat-only)
- dataset organism and reference preset/bundle mismatch

Next steps:

- set `organism` during dataset build
- for non-rat lanes, provide an explicit `ReferenceBundle`
- review [Validation Guide](validation.md#reference-validation)

## Signalome Preconditions and Runtime Issues

Common causes:

- missing or empty `site_metadata.protein_id`
- low-support rows causing downstream score limitations
- invalid values (for example infinite values) in upstream matrices

Next steps:

- ensure non-empty `protein_id` values for all interpreted sites
- inspect policy and diagnostics in [API Guide](api.md)
- review workflow seam checks in [Validation Guide](validation.md#workflow-validation)

## CLI-Specific Failures

Common causes:

- wrong file paths or unsupported file extension
- parquet used without optional parquet dependencies
- incompatible CLI flags for the selected command

Next steps:

- check command examples in [CLI Guide](cli.md)
- confirm bundle/output expectations in [Output Bundles](output_bundles.md)

## If the Error Persists

1. Capture the full error text and the command/request used.
2. Confirm whether failure is at builder, reference resolution, kinase, or signalome boundary.
3. Cross-check the relevant section in [Validation Guide](validation.md) and
   [API Guide](api.md) before opening an issue.

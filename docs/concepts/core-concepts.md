# Core Concepts

If you keep these five boundaries in mind, PhosPy becomes much easier to read.

## 1. Dataset boundary

`AnalysisReadyPhosphoDataset` is the strict workflow input.

## 2. Builder boundary

`AnalysisReadyDatasetBuilder` is the public ingestion path from DataFrames or files.

## 3. Workflow boundary

Each workflow takes one request object and returns one typed result object.

## 4. Reference boundary

Workflows use either a `ReferencePreset` or an explicit `ReferenceBundle`.

## 5. Signalome protein identity boundary

Signalome needs explicit `protein_id`. Site IDs are not a fallback.

For day-to-day use, return to the [Quickstart](../getting-started/quickstart-first-workflow.md) or [API Guide](../api.md).

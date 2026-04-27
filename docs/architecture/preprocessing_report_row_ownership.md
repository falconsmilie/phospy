# Preprocessing Report Row Ownership

This page defines the ownership contract for preprocessing report construction.

## Ownership model

- Preprocessing stages own scientific report-row emission for the sections they implement.
- The preprocessing pipeline validates and aggregates stage-emitted rows.
- The dataset builder composes those rows into `DatasetPreprocessingReport`.
- Cross-stage summary sections remain centrally assembled.

## Stage-owned sections

Supported stage-owned report tables are:

- `row_audit` (`PreprocessingRowAuditRow`)
- `duplicate_site_resolution` (`DuplicateSiteResolutionRow`)
- `metadata_conflicts` (`MetadataConflictRow`)
- `comparison_group_stats` (`ComparisonGroupStatsRow`)
- `comparison_pair_stats` (`ComparisonPairStatsRow`)

Current stage ownership:

- `missing_data`: `row_audit`
- `site_matrix`: `row_audit`, `duplicate_site_resolution`, `metadata_conflicts`
- `comparisons`: `comparison_group_stats`, `comparison_pair_stats`

## Centrally assembled sections

These sections are intentionally central (cross-stage summaries):

- `row_counts`
- `operations`
- final boundary rows for `final_dataset_construction`

## Contract behaviour

- `PreprocessingStageResult.report_rows` must use supported table names and typed row payloads.
- Unsupported table names or mismatched payload types fail fast with `DatasetBuildError`.
- There is no unbounded generic report-row bucket.

## Adding a new stage contribution

1. Emit typed rows in the stage via `PreprocessingStageResult.report_rows`.
2. Use existing report-schema row dataclasses where available.
3. If introducing a new report table, add:
   - typed row schema in `report_schema.py`
   - composition/validation wiring in `preprocessing/report_rows.py`
   - coverage in unit tests for emission, composition, and final report exposure.

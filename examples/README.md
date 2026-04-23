# Public Examples

These scripts are release-facing examples for the supported 1.5.0 public lane.
They intentionally mirror the recommended quickstart story instead of trying to
show every supported option.

## Preferred example order

1. `dataset_builder_demo.py`
   - builds an `AnalysisReadyPhosphoDataset` from in-memory `DataFrame` inputs
   - sets `organism=Organism.RAT` so the dataset is ready for bundled-reference
     execution in the recommended first run
   - keeps builder preprocessing at defaults rather than introducing advanced
     preprocessing or site-matrix branches
2. `kinase_workflow_demo.py`
   - starts from the same dataset-builder lane
   - runs `KinaseWorkflowRequest(..., references=ReferencePreset.AUTO)`
   - keeps the demo focused on the bundled-reference scoring/prediction lane
     that the quickstart introduces first
3. `signalome_workflow_demo.py`
   - continues from the same dataset -> kinase lineage
   - keeps the signalome happy path explicit by including non-empty
     `site_metadata.protein_id` values for every site

## What these examples intentionally do not show

To keep onboarding clear, the examples do not mix in alternate or advanced
lanes such as:

- caller-supplied `ReferenceBundle` wiring for non-rat execution
- advanced builder preprocessing combinations
- extra signalome tuning beyond the default supported configuration

Those routes remain available in the full API and documentation, but these
scripts are meant to feel like official first-run guidance.

## Smoke coverage

The example smoke tests treat these scripts as public artefacts. The tests check
that they run and that their console output continues to advertise the intended
1.5.0 support lane.

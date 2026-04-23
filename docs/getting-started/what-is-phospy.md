# What Is PhosPy?

PhosPy is a Python package for phosphoproteomics workflows.

In practical terms, it helps you do three things:

1. turn phosphosite data into a validated analysis-ready dataset
2. score and predict kinase activity from that dataset
3. optionally group results into signalome outputs

PhosPy does **not** expose a web server or REST API. You use it through Python
or the `phospy` CLI.

## Who this is for

PhosPy is aimed at scientists who want a reproducible workflow, plus engineers
and maintainers who need clear public boundaries.

## The supported beginner lane

For the smoothest first run:

- use a `phospho` matrix and matching `site_metadata`
- set `organism=Organism.RAT`
- run kinase with `ReferencePreset.AUTO`
- add `protein_id` only if you plan to run signalome

Bundled runtime references are rat-only in this release. Human and mouse lanes
need an explicit `ReferenceBundle` in Python.

## The public shape

The supported public shape is intentionally small:

- `AnalysisReadyDatasetBuilder().run(DatasetBuildRequest(...))`
- `KinaseWorkflow().run(KinaseWorkflowRequest(...))`
- `SignalomeWorkflow().run(SignalomeWorkflowRequest(...))`

Continue with [Quickstart: first workflow](quickstart-first-workflow.md).

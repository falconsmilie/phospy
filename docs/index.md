# PhosPy Documentation

PhosPy docs are organised to get a new user to a successful first run quickly,
without hiding the stricter contract details that advanced users and
maintainers need later.

You do **not** need to read everything in order.

## Start Here

If you are brand new to PhosPy, use this sequence:

1. [What is PhosPy?](getting-started/what-is-phospy.md)
2. [Quickstart: first workflow](getting-started/quickstart-first-workflow.md)
3. [Troubleshooting first run](getting-started/troubleshooting-first-run.md)
4. [Choose your next path](learning-paths/choose-your-path.md)

## Common Tasks

| I want to... | Go here |
| --- | --- |
| Build and run my first analysis | [Quickstart](getting-started/quickstart-first-workflow.md) |
| Check whether my inputs match the supported lane | [Quickstart](getting-started/quickstart-first-workflow.md) |
| Run workflows from the command line | [CLI Guide](cli.md) |
| Understand dataset -> kinase -> signalome flow | [Workflow guides](workflow-guides/index.md) |
| Integrate PhosPy in Python code | [API Guide](api.md) |
| Diagnose validation or runtime errors | [Troubleshooting](getting-started/troubleshooting-first-run.md) |
| Understand what 1.5.0 ships | [Release Notes 1.5.0](release_notes/1.5.0.md) |

## Docs Areas

- [Getting started](getting-started/index.md) for onboarding and first success.
- [Workflows and usage](workflow-guides/index.md) for day-to-day execution docs.
- [API and validation](reference/index.md) for strict contract details.
- [Troubleshooting](getting-started/troubleshooting-first-run.md) for common errors and recovery paths.
- [Maintainer docs](contributor/index.md) for architecture, ADRs, fixtures, and process-heavy material.

## A Good First Outcome

After the getting-started pages, you should be able to:

- prepare a valid `phospho` table and matching `site_metadata`
- build `AnalysisReadyPhosphoDataset`
- run kinase in the supported bundled rat lane
- understand when signalome needs `protein_id`

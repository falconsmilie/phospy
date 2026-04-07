# Changelog

All notable changes to this project are documented here.

## [1.2.1] - 2026-04-04

Current package metadata version reflected in this repository snapshot.

### Documentation

- simplified the README and docs pages so the main workflows are easier to find
- aligned API docs with the current public classes, methods, CLI options, and output files
- clarified where `predMat` validation happens in `PhosRPipeline`
- tightened validation and parity docs to focus on user-facing behaviour

## [1.0.0] - 2026-03-26

First supported PhosPy release.

### Scope

PhosPy 1.0.0 covers:

- core preprocessing from total and phospho inputs to corrected phosphosite matrices
- downstream kinase analysis from `predMat`
- a native `KinaseWorkflow` with seam-level validation against committed references
- a small supported root-level public API
- CLI support for the core preprocessing plus `predMat` path

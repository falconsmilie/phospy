# Roadmap

This roadmap describes likely direction, not a release promise.

PhosPy is intentionally small. The goal is not to mirror all of PhosR as quickly as possible. The goal is to grow in
ways that fit the current codebase: narrow public APIs, fixture-backed behaviour, and clear scope.

## Guiding Principles

The most credible next work is the work that keeps PhosPy:

- small at the public API boundary
- well tested against fixtures and seam-specific checks
- honest about which workflows are implemented and validated
- useful for real phosphoproteomics analysis without over-claiming package-wide parity

In practice, that means PhosPy is more likely to deepen the workflows it already supports than to chase broad feature
parity in one jump.

## Release Themes

### v1.2: Deepen the Current Package

This release would focus on making the current package easier to use and more complete without changing what it
fundamentally is.

Most likely themes:

- wider preprocessing helpers that fit the current dataset model
- better QC and downstream reporting
- native `KinaseWorkflow` CLI support
- smoother workflow inputs, outputs, and reproducible configuration handling

This is the most natural next step because PhosPy already has a clear core path:

- preprocess total and phospho tables
- build corrected phosphosite matrices
- analyse kinase activity from `predMat`
- run a native kinase workflow from scoring through prediction

The biggest near-term win is to make those workflows more ergonomic, more diagnosable, and easier to run end to end.

Likely examples of this work:

- more explicit filtering helpers instead of pushing every decision through one core path
- optional scaling or standardisation helpers where they fit naturally
- clearer QC summaries and exportable diagnostics
- richer output bundles for downstream interpretation
- command-line support for the native workflow with validated request inputs

### v1.3: Add Biological Interpretation

Once the current workflows feel smoother and more complete, the next likely step is to improve interpretation.

Most likely themes:

- pathway enrichment
- phosphosite-to-gene collapsing helpers
- annotation utilities
- packaged reference-data loaders for common workflow inputs

This is where PhosPy would start to feel less like a narrow computational port and more like a practical analysis
package.

Likely examples of this work:

- over-representation analysis from filtered phosphosite or gene sets
- rank-based enrichment from differential or activity-derived rankings
- helper utilities that add gene, site, motif, or known-substrate context
- better handling of packaged or registry-based reference inputs instead of requiring users to assemble every resource
  manually

This kind of work adds interpretation value without forcing a package-wide redesign.

### v2.0: Tackle the Heavier Statistical Work

A later milestone may take on the larger PhosR-inspired areas that carry more statistical and implementation risk.

Most likely themes:

- imputation
- SPS-related groundwork
- batch-correction workflows
- signalome or network foundations, if the validation story is strong enough

These features are valuable, but they are not good candidates for rushed implementation. They depend on solid
preprocessing, strong fixtures, clear defaults, and careful documentation.

If PhosPy grows into these areas, it should do so conservatively and with narrow claims.

Likely examples of this work:

- imputation APIs with explicit, well-documented defaults
- stronger diagnostics around missingness and filtering before correction steps
- SPS-aware data handling or reference support
- early signalome or network export layers backed by reproducible examples

## Likely Next Ports From PhosR-Inspired Work

The most natural ports are the ones that extend workflows PhosPy already supports.

### Wider Preprocessing Helpers

Version `1.0.0` focuses on the core path from total and phospho tables to corrected phosphosite matrices. A sensible
next wave is to add carefully chosen helpers around filtering, scaling, standardisation, and QC where they fit the
current dataset model cleanly.

### Richer Downstream Kinase Summaries

The current downstream analysis is centred on `predMat`, weighted activity, KSEA-style summaries, and target counts.
A natural extension is a broader reporting layer for interpretation, export, and inspection.

### Smoother Native Workflow Support

The native workflow already covers profile construction, motif scoring, score combination, candidate selection, and
adaptive SVM prediction. A practical next step is to make it easier to move between raw tables, validated request
objects, saved configuration, and reproducible output bundles.

### Pathway and Annotation Layers

Pathway enrichment and annotation are strong candidates for the next major analytical expansion because they build on
the outputs PhosPy already produces and improve interpretation without requiring a full package redesign.

## Less Likely Near-Term Work

Some directions are possible, but less likely in the near term.

### Broad PhosR Parity Claims

PhosPy is deliberately narrow. A fast attempt to claim broad parity with PhosR is less likely than incremental,
well-tested growth.

### Large Public API Expansion

A much wider public surface would make the package harder to validate and harder to support. Growth is more likely to
happen behind a small, stable set of public classes.

### Heavy Object-Model Redesign

A larger experiment-container abstraction may make sense later, but only if the current DataFrame-first approach
clearly becomes a limitation. It is not the most likely near-term change.

### Signalome Visualisation as a First Priority

Signalome and network-style outputs are attractive, but they sit on top of multiple lower-level capabilities. They are
more likely to follow stronger preprocessing, annotation, and interpretation layers than to arrive first.

## Reading This Roadmap

Treat the items above as direction rather than commitment.

The best guide to what ships next is still the same:

- does it fit the current package shape?
- does it improve real workflows?
- can it be validated with strong fixtures and clear examples?
- can PhosPy ship it without pretending to be broader than it is?

If the answer is yes, it is the kind of roadmap item that is most likely to happen.
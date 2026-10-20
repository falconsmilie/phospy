# Welcome to PhosPy

PhosPy helps scientists and analysts prepare and analyse quantitative phosphoproteomics data in Python. It works with
familiar pandas tables and records the choices that shape each analysis.

Phosphoproteomics measures phosphorylation: the addition of phosphate groups at particular amino-acid positions on
proteins. Each measured position is a **phosphosite**. Kinases are enzymes that add phosphate groups to proteins.

Data commonly arrive as an intensity table with phosphosites in rows and samples in columns. These measurements are
affected by localisation uncertainty, missing values, sample preparation, and other experimental factors.

PhosPy helps organise and analyse this evidence, but it does not replace scientific judgement. Biological interpretation
still depends on the experiment, reference data, assumptions, and independent validation.

## Start Here

The [Quickstart](quickstart.md) introduces the main PhosPy concepts in context. You can begin without understanding
request objects, configuration classes, builders, workflows, or result models.

If this is your first visit:

1. [Install PhosPy](installation.md).
2. Follow the [Quickstart](quickstart.md) from two small pandas tables to a differential-analysis result.
3. Read [Build an Analysis-Ready Dataset](api/dataset-build-workflow.md) when you are ready to prepare your own data and
   choose preprocessing policies.

## Choose and Follow an Analysis Path

Most analyses begin by preparing phosphosite data and end by inspecting the results in the context of the original
experiment. The workflow you choose between those points depends on the scientific question.

### Prepare or Import Phosphosite Data

Start with an intensity table and aligned site metadata. The intensity table contains the measured signal for each
phosphosite and sample. The metadata identifies the sites and provides the biological information needed by later steps.

If you have an exported results table, see [Import Phosphosite Data](importers.md) for the supported import routes.

### Build an Analysis-Ready Dataset

PhosPy aligns the intensity table with its metadata, applies the preprocessing choices you make, and checks that the
result satisfies the requirements for downstream analysis.

Building a dataset is therefore more than loading a DataFrame. It establishes the quantitative and biological context
that later workflows rely on.

See [Build an Analysis-Ready Dataset](api/dataset-build-workflow.md) for the guided workflow.

### Compare Experimental Conditions

Use differential analysis to compare fitted phosphosite intensities between named experimental conditions.

The result describes statistical evidence under the fitted model. An intensity difference does not by itself establish a
change in phosphorylation occupancy or explain the biological mechanism responsible for it.

See [Differential Analysis](api/differential-analysis.md).

### Test Biological-Set Overlap

Use enrichment analysis when you have selected a list of gene, protein, or phosphosite identifiers and want to test
whether it overlaps supplied biological sets more than expected.

Enrichment also requires a measured background: the identifiers that could reasonably have been selected in the
experiment. This background is essential for interpreting the overlap.

See [Enrichment](api/enrichment.md).

### Explore Kinase Support

Use kinase analysis to compare the relative support for candidate kinases across phosphosites. Kinase scores are
evidence within the analysis rather than direct proof of kinase activity.

See [Kinase Analysis](api/kinase.md).

### Describe Coordinated Signalome Patterns

Signalome analysis builds on a completed kinase analysis. It groups kinase-supported phosphosites into modules and
summarises relationships in their kinase score patterns.

These modules and associations are descriptive results. Their biological relevance must be assessed using the experiment
and appropriate supporting evidence.

See [Signalome Analysis](api/signalome.md).

### Correct Unwanted Variation

Stable phosphosite discovery identifies sites whose relative signal remains consistently small across suitable
independent reference datasets. These sites can then serve as controls for RUV-style correction of unwanted variation in
a target experiment.

See [SPS Discovery and RUV-Style Correction](api/sps-ruv.md).

### Use Suitable Reference Data

Some analyses require a **reference bundle**, which is a collection of biological reference data used during scoring and
interpretation.

Read [Reference Data](reference_bundles.md) before running a reference-dependent analysis. Bundled runtime references
are rat-only. Analyses of other organisms require an appropriate explicit reference bundle.

### Inspect Results and Provenance

Inspect the main result tables together with the diagnostics, warnings, caveats, and preprocessing reports provided by
the workflow.

**Provenance** records the inputs, processing choices, and analysis settings used to produce a dataset or result. Keep
it with the result so that another reader can understand what was analysed and how.

## Keep the Analysis Traceable

PhosPy validates inputs at dataset and workflow boundaries. Inconsistent alignment, quantitative scales, identifiers, or
workflow requirements therefore fail with explicit errors instead of silently producing a result.

Built datasets record their preprocessing provenance. Workflow results expose applicable diagnostics, caveats, and
provenance.

Documented deterministic modes are repeatable for the same inputs, analysis settings, and package version. When a
supported policy uses randomness, its seed makes that choice explicit and reproducible.

## Interpret Results Carefully

Phosphosite-intensity differences do not by themselves prove a change in the proportion of a protein population
phosphorylated at that site, known as **phosphorylation occupancy**. They also do not prove kinase activity, pathway
activation, or a biological mechanism.

Keep the following distinctions in mind:

- Differential-analysis statistics describe evidence under the fitted model.
- Enrichment describes overlap relative to the supplied background and set collection.
- Kinase scores provide relative support within an analysis.
- Signalome modules and associations are descriptive summaries derived from kinase results.

Preprocessing and validation cannot repair a confounded experimental design or make an unsuitable reference biologically
appropriate. Choose thresholds, backgrounds, controls, and reference data for your experiment rather than treating
examples in the guides as universal recommendations.

Read [Scientific Interpretation and Limitations](scientific-interpretation.md) before drawing biological
conclusions. [Scientific Coverage](scientific-coverage.md) describes the supported methods and their boundaries in more
detail.
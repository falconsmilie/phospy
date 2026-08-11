# Scientific Interpretation and Limitations

PhosPy helps you organise and analyse phosphoproteomics evidence. Its outputs
support interpretation and follow-up; they do not replace experimental design,
biological context, or independent validation.

Each workflow guide keeps its most important caveats beside the request and
response contract. This page collects the limits shared across workflows.

## Differential Analysis

Differential analysis estimates an explicit condition contrast on an
established log2 phosphosite-intensity scale. `logFC` is the fitted contrast.
`P.Value` and `adj.P.Val` describe statistical evidence under the selected
model and multiple-testing method.

A small *p* value is not an effect size, and statistical significance does not
establish biological importance. Review the effect estimate, uncertainty,
replicate structure, preprocessing history, and the number of tested sites
together.

The current workflow is limited to tested design and contrast envelopes. It is
not full limma or PhosR parity. Batch covariates are fixed model terms, not a
general batch-correction system.

## Kinase Analysis

Kinase scores express relative support within a run. Higher values indicate
stronger support under the selected method, reference data, and thresholds.
They are not calibrated probabilities and do not prove that a kinase caused an
observed phosphosite change.

Optional activity outputs summarise substrate-supported patterns. Treat them as
prioritisation evidence, not direct proof of kinase activation or causal
pathway activity.

## Signalome Analysis

Signalome modules and network-style edges are summaries derived from upstream
kinase outputs. They describe score-supported groupings and score-profile
associations within the analyzed data.

They do not establish causality, physical interaction, direction of signaling,
or experimental validation of a pathway relationship. Module labels are local
to the run and should not be treated as universal biological classes.

## Enrichment

`EnrichmentWorkflow` performs offline over-representation analysis (ORA) with
caller-supplied identifiers, gene sets, and an explicit background universe.
Results depend strongly on identifier mapping, the selected background, and the
contents of the supplied collection.

Overlap statistics do not prove pathway activation, regulation, or causality.
ORA does not imply GSEA or PTM-SEA support.

## References and Organisms

Bundled runtime references are rat-only in this release. Human, mouse, and
custom analyses must use an explicit `ReferenceBundle`. Keep the source,
version, license, identifier namespace, and organism with the analysis record.

A technically compatible reference is not automatically the best biological
reference for every tissue, condition, or experimental platform.

## Continue Reading

- [Scientific Coverage](scientific-coverage.md)
- [Reference Data](reference_bundles.md)
- [Differential Analysis](api/differential-analysis.md)
- [Enrichment](api/enrichment.md)
- [Kinase Analysis](api/kinase.md)
- [Signalome Analysis](api/signalome.md)

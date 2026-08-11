# Scientific interpretation and limitations

This page collects the main interpretation limits that apply across PhosPy
workflows. Each workflow page keeps the workflow-specific caveats beside its
request and result contract.

## Differential analysis

Differential analysis estimates explicit condition contrasts on an established
log2 phosphosite intensity scale. `logFC` is the fitted contrast on that scale;
`P.Value` and `adj.P.Val` describe statistical evidence for the tested contrast
under the configured model and multiple-testing method.

The current differential workflow is limited to tested design and contrast
envelopes and is not full limma or PhosR parity. Batch covariates are ordinary
fixed-effect model terms, not broad batch correction.

## Kinase analysis

Kinase scoring and prediction outputs are relative support within a run. Higher
scores mean stronger support under the selected scoring and reference context;
they are not calibrated probabilities and are not causal proof.

Optional kinase activity outputs are substrate-supported summaries. They are
useful for prioritising follow-up, but they are not direct proof of kinase
activation or causal pathway activity.

## Signalome analysis

Signalome modules and network-style edges are derived summaries from kinase
outputs. They describe score-supported groupings and score-profile
associations. They do not establish causality, physical interaction, direction,
or experimental validation of signalling relationships.

## Enrichment

`EnrichmentWorkflow` runs offline over-representation analysis (ORA) using
caller-supplied selected identifiers, set collections, and an explicit
background universe. ORA overlap statistics do not prove pathway activation,
regulation, or biological causality. ORA does not imply GSEA or PTM-SEA
support.

## References and organisms

Bundled runtime references are rat-only in this release. For human, mouse, or
custom reference contexts, pass an explicit `ReferenceBundle` and keep the
reference source, version, and identifier namespace with your analysis record.

## More detail

- [Scientific Coverage](scientific-coverage.md)
- [Reference Bundles](reference_bundles.md)
- [Differential analysis](api/differential-analysis.md)
- [Enrichment](api/enrichment.md)
- [Kinase analysis](api/kinase.md)
- [Signalome analysis](api/signalome.md)

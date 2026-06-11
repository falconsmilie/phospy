# Public Examples

These scripts show the recommended beginner lane.

Read them in this order:

1. `dataset_builder_demo.py`
2. `kinase_workflow_demo.py`
3. `signalome_workflow_demo.py`

Use `reference_bundle_builder_demo.py` when you need a small local-file example
for building an explicit human/mouse-style `ReferenceBundle`.

They intentionally keep the story small:

- rat bundled-reference lane first
- explicit localisation policy (`localisation_confidence >= 0.75`) at dataset build
- `phospy.api` for requests and enums
- explicit protein context (`organism`, `protein_namespace`, and
  `protein_identifier`) so the builder can derive `site_key`
- explicit `protein_id` where the signalome lane needs it
- output tables that expose `site_key`, `display_id`, `gene_symbol`, `site`, and
  protein-context metadata
- explicit reporting of both numeric scale and quantitative meaning (for example, `log2` with `phosphosite_log_abundance`)

They are not meant to show every advanced option. For that, use the main docs.

# Differential Contract Negative-Case Fixture Provenance

Reference source: PhosPy differential workflow contract and ADR-0019
(`docs/adr/adr_0019_experimental_design_and_contrast_contract.md`).

These fixtures are intentionally not limma output fixtures. They encode
unsupported-input and boundary-rejection scenarios that must fail clearly in
PhosPy:

- rank-deficient design matrix at core executor boundary
- missing values rejected before differential execution at
  `AnalysisReadyPhosphoDataset` boundary

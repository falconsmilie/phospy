from __future__ import annotations


def test_duplicate_correlation_scientific_contracts_are_not_public_exports() -> None:
    import phospy
    import phospy.advanced as advanced_api
    import phospy.api as public_api
    import phospy.api.configs as api_configs
    import phospy.contracts.configs as contract_configs
    import phospy.contracts.results as contract_results
    import phospy.science.differential.models as differential_models

    duplicate_correlation_symbols = {
        "DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN",
        "DUPLICATE_CORRELATION_TRIM_FRACTION",
        "DuplicateCorrelationBlockStructureSummary",
        "DuplicateCorrelationConsensusResult",
        "DuplicateCorrelationConsensusSummary",
        "DuplicateCorrelationFailureReason",
        "DuplicateCorrelationFeatureEstimate",
        "DuplicateCorrelationFeatureStatus",
        "DuplicateCorrelationWorkflowProvenance",
        "PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION",
    }
    public_modules = (
        phospy,
        public_api,
        advanced_api,
        api_configs,
        contract_configs,
        contract_results,
        differential_models,
    )

    assert "duplicate_correlation" not in (
        contract_configs.SUPPORTED_PAIRED_DESIGN_POLICIES
    )
    for module in public_modules:
        exported_names = set(getattr(module, "__all__", ()))
        assert duplicate_correlation_symbols.isdisjoint(exported_names)
        assert all(
            not hasattr(module, symbol_name)
            for symbol_name in duplicate_correlation_symbols
        )

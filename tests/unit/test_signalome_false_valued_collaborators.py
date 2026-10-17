from __future__ import annotations

from dataclasses import dataclass

import pytest

from phospy.quantitative.validation import WorkflowQuantitativeInputValidator
from phospy.signalome.alignment_diagnostics import (
    SignalomeAlignmentDiagnosticsBuilder,
)
from phospy.signalome.clustering_runner import SignalomeClusteringRunner
from phospy.signalome.context_tables import SignalomeContextTableBuilder
from phospy.signalome.executor import SignalomeWorkflowExecutor
from phospy.signalome.interpreter import SignalomeWorkflowInterpreter
from phospy.signalome.matrix_alignment import SignalomeMatrixAligner
from phospy.signalome.module_tables import SignalomeModuleTableBuilder
from phospy.signalome.network_builder import SignalomeNetworkBuilder
from phospy.signalome.protein_resolution import SignalomeProteinResolver
from phospy.signalome.provenance import SignalomeProvenanceBuilder
from phospy.signalome.public import SignalomeWorkflow
from phospy.signalome.result_assembly import (
    SignalomeExpandedSignalomeBuilder,
    SignalomeResultAssembler,
)
from phospy.signalome.score_matrix_selection import SignalomeScoreMatrixSelector
from phospy.signalome.score_preconditioning import SignalomeScorePreconditioner
from phospy.signalome.validator import (
    SignalomeConfigValidator,
    SignalomeWorkflowValidator,
)


@dataclass(frozen=True)
class _FalsyCollaborator:
    def __bool__(self) -> bool:
        return False


@pytest.mark.parametrize(
    ("factory", "attribute", "expected_type"),
    (
        (
            lambda value: SignalomeWorkflow._with_components(validator=value),
            "_validator",
            SignalomeWorkflowValidator,
        ),
        (
            lambda value: SignalomeWorkflow._with_components(interpreter=value),
            "_interpreter",
            SignalomeWorkflowInterpreter,
        ),
        (
            lambda value: SignalomeWorkflow._with_components(executor=value),
            "_executor",
            SignalomeWorkflowExecutor,
        ),
        (
            lambda value: SignalomeWorkflowValidator(config_validator=value),
            "_config_validator",
            SignalomeConfigValidator,
        ),
        (
            lambda value: SignalomeWorkflowValidator(
                quantitative_input_validator=value
            ),
            "_quantitative_input_validator",
            WorkflowQuantitativeInputValidator,
        ),
        (
            lambda value: SignalomeWorkflowInterpreter(score_matrix_selector=value),
            "_score_matrix_selector",
            SignalomeScoreMatrixSelector,
        ),
        (
            lambda value: SignalomeWorkflowInterpreter(matrix_aligner=value),
            "_matrix_aligner",
            SignalomeMatrixAligner,
        ),
        (
            lambda value: SignalomeWorkflowInterpreter(score_preconditioner=value),
            "_score_preconditioner",
            SignalomeScorePreconditioner,
        ),
        (
            lambda value: SignalomeWorkflowInterpreter(protein_resolver=value),
            "_protein_resolver",
            SignalomeProteinResolver,
        ),
        (
            lambda value: SignalomeWorkflowInterpreter(
                alignment_diagnostics_builder=value
            ),
            "_alignment_diagnostics_builder",
            SignalomeAlignmentDiagnosticsBuilder,
        ),
        (
            lambda value: SignalomeWorkflowExecutor(clustering_runner=value),
            "_clustering_runner",
            SignalomeClusteringRunner,
        ),
        (
            lambda value: SignalomeWorkflowExecutor(module_table_builder=value),
            "_module_table_builder",
            SignalomeModuleTableBuilder,
        ),
        (
            lambda value: SignalomeWorkflowExecutor(network_builder=value),
            "_network_builder",
            SignalomeNetworkBuilder,
        ),
        (
            lambda value: SignalomeWorkflowExecutor(context_table_builder=value),
            "_context_table_builder",
            SignalomeContextTableBuilder,
        ),
        (
            lambda value: SignalomeWorkflowExecutor(expanded_signalome_builder=value),
            "_expanded_signalome_builder",
            SignalomeExpandedSignalomeBuilder,
        ),
        (
            lambda value: SignalomeWorkflowExecutor(provenance_builder=value),
            "_provenance_builder",
            SignalomeProvenanceBuilder,
        ),
        (
            lambda value: SignalomeWorkflowExecutor(result_assembler=value),
            "_result_assembler",
            SignalomeResultAssembler,
        ),
    ),
)
def test_false_valued_optional_collaborators_preserve_default_fallbacks(
    factory, attribute: str, expected_type: type[object]
) -> None:
    instance = factory(_FalsyCollaborator())

    assert isinstance(getattr(instance, attribute), expected_type)

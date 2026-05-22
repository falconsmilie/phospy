SHELL := /usr/bin/env bash
.SHELLFLAGS := -eu -o pipefail -c
.DEFAULT_GOAL := help

PYTHON ?= python
PIP ?= $(PYTHON) -m pip
PYTEST ?= $(PYTHON) -m pytest
RUFF ?= $(PYTHON) -m ruff
BUILD ?= $(PYTHON) -m build
PRE_COMMIT ?= $(PYTHON) -m pre_commit
RSCRIPT ?= Rscript
MKDIR_P ?= mkdir -p
RM ?= rm -rf

TRACE_KINASES ?= PRKAA1,MAPK1
TRACE_TOP_N ?= 10
FIXTURES_ROOT ?= tests/fixtures
REWRITE_PARITY_ROOT ?= $(FIXTURES_ROOT)/rewrite_parity
R_L6_OUTDIR ?= $(REWRITE_PARITY_ROOT)/r_reference_l6
PUBLIC_WORKFLOW_OUTDIR ?= $(FIXTURES_ROOT)/public_workflow_reference
ACTIVE_SCRIPTS_DIR ?= scripts/active

.PHONY: help \
	check-tools check-r-tools fixtures-dirs \
	install install-dev lint format type-check pre-commit test tests-all test-unit test-parity test-performance test-release-gate test-seams build clean \
	fixtures fixtures-r-l6 traces-r \
	fixtures-public-workflow-reference fixtures-provenance-goldens fixtures-all \
	dataset-builder-demo kinase-workflow-demo signalome-workflow-demo demo-all

help:
	@echo Available targets:
	@echo   make install                       Install the package in editable mode
	@echo   make install-dev                   Install editable package with dev and test extras
	@echo   make lint                          Run Ruff checks
	@echo   make format                        Run Ruff formatter
	@echo   make type-check                    Run the same Pyright entrypoint used by CI and pre-commit
	@echo   make pre-commit                    Run all pre-commit hooks
	@echo   make test-unit                     Run the non-parity pytest suite
	@echo   make test-parity                   Run the parity pytest suite
	@echo   make test-performance              Run the performance contract suite
	@echo   make test-release-gate             Run release validation (unit/integration, reproducibility goldens, parity, performance)
	@echo   make test                          Run unit and parity tests
	@echo   make tests-all                     Alias for all-tests
	@echo   make test-seams                    Run seam-focused rewrite parity tests
	@echo   make dataset-builder-demo          Run examples.dataset_builder_demo.main()
	@echo   make kinase-workflow-demo          Run examples.kinase_workflow_demo.main()
	@echo   make signalome-workflow-demo       Run examples.signalome_workflow_demo.main()
	@echo   make build                         Build source and wheel distributions
	@echo   make clean                         Remove common local build and test artefacts
	@echo   make fixtures-r-l6                 Generate the main L6 R-backed fixture family
	@echo   make traces-r                      Regenerate the committed R L6 prediction trace
	@echo   make fixtures-public-workflow-reference Regenerate public workflow signalome fixtures
	@echo   make fixtures-provenance-goldens   Regenerate provenance golden hash fixtures
	@echo   make fixtures-all                  Bootstrap active maintainer fixture families from scratch
	@echo   make fixtures                      Alias for fixtures-all

check-tools:
	@command -v "$(PYTHON)" >/dev/null 2>&1 || { printf 'Python executable not found: %s\n' "$(PYTHON)" >&2; exit 1; }

check-r-tools:
	@command -v "$(RSCRIPT)" >/dev/null 2>&1 || { printf 'Rscript executable not found: %s\n' "$(RSCRIPT)" >&2; exit 1; }

fixtures-dirs:
	@$(MKDIR_P) "$(FIXTURES_ROOT)" "$(REWRITE_PARITY_ROOT)"

install: check-tools
	$(PIP) install -e .

install-dev: check-tools
	$(PIP) install -e ".[dev,test]"

lint: check-tools
	$(RUFF) check .

format: check-tools
	$(RUFF) format .

pre-commit: check-tools
	$(PRE_COMMIT) run --all-files

type-check: check-tools
	$(PYTHON) scripts/run_pyright.py

test: test-unit test-parity

tests-all: check-tools
	$(PYTEST) -o addopts= tests/

test-unit: check-tools
	$(PYTEST) -m "not parity"

test-parity: check-tools
	$(PYTEST) tests/parity -m parity -s

test-performance: check-tools
	$(PYTEST) tests/performance -m "performance or release_gate"

test-release-gate: check-tools
	$(PYTEST) tests/unit tests/integration -m "not parity and not performance and not release_gate"
	$(PYTEST) tests/unit/test_provenance_regressions.py tests/integration/test_kinase_workflow_integration.py::test_kinase_public_predmat_provenance_matches_golden_contract tests/integration/test_signalome_workflow_integration.py::test_signalome_l6_provenance_matches_golden_contract -m "release_gate and (reproducibility or golden)"
	$(PYTEST) tests/parity -m "parity and not parity_diagnostic" -s
	$(PYTEST) tests/performance -m "performance or release_gate" -q

dataset-builder-demo: check-tools
	PYTHONPATH=src $(PYTHON) -c "from examples.dataset_builder_demo import main; main()"

kinase-workflow-demo: check-tools
	PYTHONPATH=src $(PYTHON) -c "from examples.kinase_workflow_demo import main; main()"

signalome-workflow-demo: check-tools
	PYTHONPATH=src $(PYTHON) -c "from examples.signalome_workflow_demo import main; main()"

demo-all: dataset-builder-demo kinase-workflow-demo signalome-workflow-demo

fixtures: fixtures-all

fixtures-r-l6: check-r-tools fixtures-dirs
	$(RSCRIPT) $(ACTIVE_SCRIPTS_DIR)/generate_r_l6_fixtures.R --outdir "$(R_L6_OUTDIR)"

traces-r: check-r-tools fixtures-dirs
	$(RSCRIPT) $(ACTIVE_SCRIPTS_DIR)/generate_r_l6_fixtures.R --outdir "$(R_L6_OUTDIR)" --trace_kinases "$(TRACE_KINASES)" --trace_top_n "$(TRACE_TOP_N)"

fixtures-public-workflow-reference: check-tools fixtures-dirs
	$(PYTHON) $(ACTIVE_SCRIPTS_DIR)/generate_signalome_public_workflow_reference.py --outdir "$(PUBLIC_WORKFLOW_OUTDIR)"

fixtures-provenance-goldens: check-tools fixtures-dirs
	$(PYTHON) $(ACTIVE_SCRIPTS_DIR)/generate_provenance_goldens.py

test-seams: check-tools
	$(PYTEST) -q \
		tests/parity/test_prediction_science_parity.py \
		tests/parity/test_adaptive_prediction_parity.py \
		tests/parity/test_adaptive_replay_parity.py

fixtures-all: fixtures-r-l6 fixtures-public-workflow-reference

build: check-tools
	$(BUILD)

clean:
	$(RM) .pytest_cache .ruff_cache build dist .eggs
	find . -type d \( -name '__pycache__' -o -name '*.egg-info' \) -prune -exec $(RM) {} +
	find . -type f \( -name '*.pyc' -o -name '*.pyo' \) -delete

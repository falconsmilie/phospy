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
PYTHON_TRACE_ROOT ?= $(FIXTURES_ROOT)/python_reference_l6
PYTHON_TRACE_OUTDIR ?= $(PYTHON_TRACE_ROOT)/prediction_trace
R_L6_OUTDIR ?= $(FIXTURES_ROOT)/r_reference_l6
R_SMALL_OUTDIR ?= $(FIXTURES_ROOT)/r_reference
FRAGILE_OUTDIR ?= $(FIXTURES_ROOT)/fragile_support_reference
L6_STRESS_OUTDIR ?= $(FIXTURES_ROOT)/r_reference_l6_seam_stress
SYNTHETIC_EDGE_OUTDIR ?= $(FIXTURES_ROOT)/synthetic_adaptive_sampling_edge
PUBLIC_WORKFLOW_OUTDIR ?= $(FIXTURES_ROOT)/public_workflow_reference

.PHONY: help \
	check-tools check-r-tools fixtures-dirs \
	install install-dev lint format pre-commit test test-unit test-parity test-seams build clean \
	fixtures fixtures-r-small fixtures-r-l6 traces-r fixtures-fragile fixtures-r-l6-seam-stress \
	traces-python traces-python-replay fixtures-synthetic-edge fixtures-public-workflow-reference fixtures-all \
	dataset-builder-demo kinase-workflow-demo signalome-workflow-demo demo-all

help:
	@echo Available targets:
	@echo   make install                       Install the package in editable mode
	@echo   make install-dev                   Install editable package with dev and test extras
	@echo   make lint                          Run Ruff checks
	@echo   make format                        Run Ruff formatter
	@echo   make pre-commit                    Run all pre-commit hooks
	@echo   make test-unit                     Run the non-parity pytest suite
	@echo   make test-parity                   Run the parity pytest suite
	@echo   make test                          Run unit and parity tests
	@echo   make test-seams                    Run seam-focused rewrite parity tests
	@echo   make dataset-builder-demo          Run examples.dataset_builder_demo.main()
	@echo   make kinase-workflow-demo          Run examples.kinase_workflow_demo.main()
	@echo   make signalome-workflow-demo       Run examples.signalome_workflow_demo.main()
	@echo   make build                         Build source and wheel distributions
	@echo   make clean                         Remove common local build and test artefacts
	@echo   make fixtures-r-small              Generate the small R-backed fixture family
	@echo   make fixtures-r-l6                 Generate the main L6 R-backed fixture family
	@echo   make traces-r                      Regenerate the committed R L6 prediction trace
	@echo   make fixtures-fragile              Generate the curated fragile-support seam fixture
	@echo   make fixtures-r-l6-seam-stress     Generate the smaller R-backed L6 seam-stress fixture
	@echo   make traces-python                 Export Python prediction traces
	@echo   make traces-python-replay          Export Python traces replaying R sampling rows
	@echo   make fixtures-synthetic-edge       Generate the synthetic adaptive-sampling edge fixture
	@echo   make fixtures-public-workflow-reference Regenerate public workflow signalome fixtures
	@echo   make fixtures-all                  Bootstrap every committed fixture family from scratch
	@echo   make fixtures                      Alias for fixtures-all

check-tools:
	@command -v "$(PYTHON)" >/dev/null 2>&1 || { printf 'Python executable not found: %s\n' "$(PYTHON)" >&2; exit 1; }

check-r-tools:
	@command -v "$(RSCRIPT)" >/dev/null 2>&1 || { printf 'Rscript executable not found: %s\n' "$(RSCRIPT)" >&2; exit 1; }

fixtures-dirs:
	@$(MKDIR_P) "$(FIXTURES_ROOT)" "$(PYTHON_TRACE_ROOT)"

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

test: test-unit test-parity

test-unit: check-tools
	$(PYTEST) -m "not parity"

test-parity: check-tools
	$(PYTEST) tests/parity -m parity -s

dataset-builder-demo: check-tools
	PYTHONPATH=src $(PYTHON) -c "from examples.dataset_builder_demo import main; main()"

kinase-workflow-demo: check-tools
	PYTHONPATH=src $(PYTHON) -c "from examples.kinase_workflow_demo import main; main()"

signalome-workflow-demo: check-tools
	PYTHONPATH=src $(PYTHON) -c "from examples.signalome_workflow_demo import main; main()"

demo-all: dataset-builder-demo kinase-workflow-demo signalome-workflow-demo

fixtures: fixtures-all

fixtures-r-small: check-r-tools fixtures-dirs
	$(RSCRIPT) scripts/generate_r_fixtures.R --outdir "$(R_SMALL_OUTDIR)"

fixtures-r-l6: check-r-tools fixtures-dirs
	$(RSCRIPT) scripts/generate_r_l6_fixtures.R --outdir "$(R_L6_OUTDIR)"

traces-r: check-r-tools fixtures-dirs
	$(RSCRIPT) scripts/generate_r_l6_fixtures.R --outdir "$(R_L6_OUTDIR)" --trace_kinases "$(TRACE_KINASES)" --trace_top_n "$(TRACE_TOP_N)"

fixtures-fragile: check-tools fixtures-r-l6
	$(PYTHON) scripts/generate_fragile_support_reference.py --source-dir "$(R_L6_OUTDIR)" --outdir "$(FRAGILE_OUTDIR)"

fixtures-r-l6-seam-stress: check-tools fixtures-r-l6
	$(PYTHON) scripts/generate_l6_seam_stress_reference.py --outdir "$(L6_STRESS_OUTDIR)"

traces-python: check-tools fixtures-r-l6 fixtures-dirs
	$(PYTHON) scripts/export_python_prediction_traces.py \
		--trace-kinases "$(TRACE_KINASES)" \
		--svm-mode r_parity \
		--debug-top-n "$(TRACE_TOP_N)" \
		--outdir "$(PYTHON_TRACE_OUTDIR)"

traces-python-replay: check-tools fixtures-r-l6 traces-python fixtures-dirs
	$(PYTHON) scripts/export_python_prediction_traces.py \
		--trace-kinases "$(TRACE_KINASES)" \
		--svm-mode r_parity \
		--debug-top-n "$(TRACE_TOP_N)" \
		--sampling-trace-dir "$(R_L6_OUTDIR)/prediction_trace" \
		--outdir "$(PYTHON_TRACE_OUTDIR)"

fixtures-synthetic-edge: check-tools fixtures-dirs
	$(PYTHON) scripts/generate_synthetic_adaptive_sampling_edge_fixtures.py --outdir "$(SYNTHETIC_EDGE_OUTDIR)"

fixtures-public-workflow-reference: check-tools fixtures-dirs
	$(PYTHON) scripts/generate_signalome_public_workflow_reference.py --outdir "$(PUBLIC_WORKFLOW_OUTDIR)"

test-seams: check-tools
	$(PYTEST) -q \
		tests/parity/test_prediction_science_parity.py \
		tests/parity/test_adaptive_prediction_parity.py \
		tests/parity/test_adaptive_replay_parity.py

fixtures-all: fixtures-r-small traces-python fixtures-fragile fixtures-r-l6-seam-stress fixtures-synthetic-edge fixtures-public-workflow-reference

build: check-tools
	$(BUILD)

clean:
	$(RM) .pytest_cache .ruff_cache build dist .eggs
	find . -type d \( -name '__pycache__' -o -name '*.egg-info' \) -prune -exec $(RM) {} +
	find . -type f \( -name '*.pyc' -o -name '*.pyo' \) -delete

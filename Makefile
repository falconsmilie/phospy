ifeq ($(OS),Windows_NT)
BASH ?= C:/Program Files/Git/bin/bash.exe
SHELL := $(BASH)
else
SHELL := /usr/bin/env bash
endif
.SHELLFLAGS := -eu -o pipefail -c
.DEFAULT_GOAL := help

PYTHON ?= python
PIP ?= $(PYTHON) -m pip
PYTEST ?= $(PYTHON) -m pytest
RUFF ?= $(PYTHON) -m ruff
BUILD ?= $(PYTHON) -m build --no-isolation
PRE_COMMIT ?= $(PYTHON) -m pre_commit
MKDOCS ?= $(PYTHON) -m mkdocs
RSCRIPT ?= Rscript
MKDIR_P ?= mkdir -p
RM_RF ?= rm -rf

TRACE_KINASES ?= PRKAA1,MAPK1
TRACE_TOP_N ?= 10
FIXTURES_ROOT ?= tests/fixtures
REWRITE_PARITY_ROOT ?= $(FIXTURES_ROOT)/rewrite_parity
R_L6_OUTDIR ?= $(REWRITE_PARITY_ROOT)/r_reference_l6
PUBLIC_WORKFLOW_OUTDIR ?= $(FIXTURES_ROOT)/public_workflow_reference
RELEASE_VALIDATION_OUTDIR ?= $(FIXTURES_ROOT)/release_validation_regression
RELEASE_VALIDATION_SEED ?= 20260724
RELEASE_VALIDATION_TIMESTAMP ?= 2026-07-24T00:00:00Z
LARGE_DIFFERENTIAL_LIMMA_TREND_OUTDIR ?= $(REWRITE_PARITY_ROOT)/differential_limma_trend_large
LARGE_DIFFERENTIAL_LIMMA_TREND_FEATURES ?= 1600
ACTIVE_SCRIPTS_DIR ?= scripts/active
PYTEST_DURATION_ARGS ?= --durations=25 --durations-min=0.01
PYTEST_REPORT_DIR ?= build/reports
TWINE ?= $(PYTHON) -m twine

.PHONY: help \
	check-tools check-r-tools fixtures-dirs \
	install install-dev lint format type-check pre-commit test tests-all test-unit test-contract test-parity test-performance test-release-gates docs-build validate-reference-bundles release-check benchmark-release-scale test-seams build clean \
	verify-installed-distributions \
	fixtures fixtures-r-l6 traces-r \
	fixtures-public-workflow-reference fixtures-provenance-goldens fixtures-release-validation-regression fixtures-large-differential-limma-trend fixtures-all \
	dataset-builder-demo kinase-workflow-demo signalome-workflow-demo demo-all

help:
	@printf '%s\n' 'Available targets:'
	@printf '%s\n' '  make install                       Install the package in editable mode'
	@printf '%s\n' '  make install-dev                   Install editable package with dev and test extras'
	@printf '%s\n' '  make lint                          Run Ruff checks'
	@printf '%s\n' '  make format                        Run Ruff formatter'
	@printf '%s\n' '  make type-check                    Run the same Pyright entrypoint used by CI and pre-commit'
	@printf '%s\n' '  make pre-commit                    Run all pre-commit hooks'
	@printf '%s\n' '  make test-unit                     Run the non-parity pytest suite'
	@printf '%s\n' '  make test-contract                 Run external-consumer public API contract tests'
	@printf '%s\n' '  make test-parity                   Run the parity pytest suite'
	@printf '%s\n' '  make test-performance              Run the performance contract suite'
	@printf '%s\n' '  make test-release-gates            Run release/golden/reproducibility gates'
	@printf '%s\n' '  make docs-build                    Build documentation with strict link/error checks'
	@printf '%s\n' '  make validate-reference-bundles    Validate checked-in reference bundle manifests and files'
	@printf '%s\n' '  make verify-installed-distributions Install and execute the built wheel and sdist outside the checkout'
	@printf '%s\n' '  make release-check                 Run maintainer release checks'
	@printf '%s\n' '  make test                          Run unit, contract, and parity tests'
	@printf '%s\n' '  make tests-all                     Alias for all-tests'
	@printf '%s\n' '  make test-seams                    Run seam-focused rewrite parity tests'
	@printf '%s\n' '  make benchmark-release-scale       Optional local 50,000x48 builder+differential benchmark'
	@printf '%s\n' '  make dataset-builder-demo          Run examples.dataset_builder_demo.main()'
	@printf '%s\n' '  make kinase-workflow-demo          Run examples.kinase_workflow_demo.main()'
	@printf '%s\n' '  make signalome-workflow-demo       Run examples.signalome_workflow_demo.main()'
	@printf '%s\n' '  make build                         Build and validate source/wheel distributions'
	@printf '%s\n' '  make clean                         Remove common local build and test artefacts'
	@printf '%s\n' '  make fixtures-r-l6                 Generate the main L6 R-backed fixture family'
	@printf '%s\n' '  make traces-r                      Regenerate the committed R L6 prediction trace'
	@printf '%s\n' '  make fixtures-public-workflow-reference Regenerate public workflow signalome fixtures'
	@printf '%s\n' '  make fixtures-provenance-goldens   Regenerate provenance golden hash fixtures'
	@printf '%s\n' '  make fixtures-release-validation-regression Regenerate compact release-validation regression fixtures'
	@printf '%s\n' '  make fixtures-large-differential-limma-trend Regenerate large R/limma trend parity fixture'
	@printf '%s\n' '  make fixtures-all                  Bootstrap active maintainer fixture families from scratch'
	@printf '%s\n' '  make fixtures                      Alias for fixtures-all'

check-tools:
	@command -v "$(PYTHON)" >/dev/null 2>&1 || { printf 'Python executable not found: %s\n' "$(PYTHON)" >&2; exit 1; }

check-r-tools:
	@command -v "$(RSCRIPT)" >/dev/null 2>&1 || { printf 'Rscript executable not found: %s\n' "$(RSCRIPT)" >&2; exit 1; }

fixtures-dirs:
	@$(MKDIR_P) "$(FIXTURES_ROOT)" "$(REWRITE_PARITY_ROOT)"

install: check-tools
	$(PIP) install -e .

install-dev: check-tools
	$(PIP) install -e ".[dev,test,docs]"

lint: check-tools
	$(RUFF) check .

format: check-tools
	$(RUFF) format .

pre-commit: check-tools
	$(PRE_COMMIT) run --all-files

type-check: check-tools
	$(PYTHON) scripts/run_pyright.py

test: test-unit test-contract test-parity

tests-all: check-tools
	$(PYTEST) -o addopts= tests/

test-unit: check-tools
	$(PYTEST) -m "not parity"

test-contract: check-tools
	$(MKDIR_P) "$(PYTEST_REPORT_DIR)"
	$(PYTEST) -o addopts= tests/contract --junitxml "$(PYTEST_REPORT_DIR)/contract.xml"

test-parity: check-tools
	$(PYTEST) tests/parity -m "parity and not parity_diagnostic" -s

test-performance: check-tools
	$(MKDIR_P) "$(PYTEST_REPORT_DIR)"
	$(PYTEST) $(PYTEST_DURATION_ARGS) tests/performance -m "performance or release_gate" --junitxml "$(PYTEST_REPORT_DIR)/performance.xml"

test-release-gates: check-tools
	$(MKDIR_P) "$(PYTEST_REPORT_DIR)"
	$(PYTEST) -o addopts= tests/release tests/golden \
		-m "release_gate or golden or reproducibility" \
		--junitxml "$(PYTEST_REPORT_DIR)/release-gates.xml"

validate-reference-bundles: check-tools
	$(PYTHON) scripts/validate_reference_bundle_index.py --repo-root .

docs-build: check-tools
	$(MKDOCS) build --strict

release-check: lint type-check test-unit test-contract test-parity test-performance docs-build validate-reference-bundles test-release-gates build verify-installed-distributions

benchmark-release-scale: check-tools
	$(PYTHON) benchmarks/measure_release_scale_builder_differential.py

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

fixtures-release-validation-regression: check-tools fixtures-dirs
	$(PYTHON) $(ACTIVE_SCRIPTS_DIR)/generate_release_validation_regression_fixtures.py --outdir "$(RELEASE_VALIDATION_OUTDIR)" --timestamp "$(RELEASE_VALIDATION_TIMESTAMP)" --seed "$(RELEASE_VALIDATION_SEED)"

fixtures-large-differential-limma-trend: check-r-tools fixtures-dirs
	$(RSCRIPT) $(ACTIVE_SCRIPTS_DIR)/generate_large_differential_limma_trend_fixture.R --outdir "$(LARGE_DIFFERENTIAL_LIMMA_TREND_OUTDIR)" --seed "$(RELEASE_VALIDATION_SEED)" --timestamp "$(RELEASE_VALIDATION_TIMESTAMP)" --n_features "$(LARGE_DIFFERENTIAL_LIMMA_TREND_FEATURES)"

test-seams: check-tools
	$(PYTEST) -q \
		tests/parity/test_prediction_science_parity.py \
		tests/parity/test_adaptive_prediction_parity.py \
		tests/parity/test_adaptive_replay_parity.py

fixtures-all: fixtures-r-l6 fixtures-public-workflow-reference fixtures-release-validation-regression fixtures-large-differential-limma-trend

build: check-tools
	$(RM_RF) dist
	$(BUILD)
	@shopt -s nullglob; wheels=(dist/*.whl); sdists=(dist/*.tar.gz); \
	if (( $${#wheels[@]} != 1 )); then printf 'Expected exactly one wheel in dist/, found %s\n' "$${#wheels[@]}" >&2; exit 1; fi; \
	if (( $${#sdists[@]} != 1 )); then printf 'Expected exactly one sdist in dist/, found %s\n' "$${#sdists[@]}" >&2; exit 1; fi
	$(TWINE) check dist/*
	$(PYTHON) scripts/validate_reference_bundle_distribution.py --no-git-index-compare dist/*

verify-installed-distributions: check-tools build
	$(PYTHON) scripts/verify_installed_distributions.py --dist-dir dist --repo-root . --constraint constraints/ci.txt

clean:
	$(RM_RF) .pytest_cache .ruff_cache build dist .eggs
	find . -type d \( -name '__pycache__' -o -name '*.egg-info' \) -prune -exec $(RM_RF) {} +
	find . -type f \( -name '*.pyc' -o -name '*.pyo' \) -delete

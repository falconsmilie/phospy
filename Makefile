PYTHON ?= python
RSCRIPT ?= Rscript

R_SMALL_OUTDIR ?= tests/fixtures/r_reference
R_L6_OUTDIR ?= tests/fixtures/r_reference_l6
FRAGILE_OUTDIR ?= tests/fixtures/fragile_support_reference
PYTHON_TRACE_OUTDIR ?= tests/fixtures/python_reference_l6/prediction_trace
SYNTHETIC_EDGE_OUTDIR ?= tests/fixtures/synthetic_adaptive_sampling_edge

TRACE_KINASES ?= PRKAA1,MAPK1
TRACE_TOP_N ?= 10
SVM_MODE ?= r_parity

.PHONY: help \
	fixtures-r-small \
	fixtures-r-l6 \
	traces-r \
	fixtures-fragile \
	traces-python \
	traces-python-replay \
	fixtures-synthetic-edge \
	fixtures-all

help:
	@printf "Available targets:\n\n"
	@printf "  %-24s %s\n" "fixtures-r-small" "Generate tests/fixtures/r_reference"
	@printf "  %-24s %s\n" "fixtures-r-l6" "Generate tests/fixtures/r_reference_l6"
	@printf "  %-24s %s\n" "traces-r" "Generate R L6 prediction traces"
	@printf "  %-24s %s\n" "fixtures-fragile" "Generate tests/fixtures/fragile_support_reference"
	@printf "  %-24s %s\n" "traces-python" "Generate Python prediction traces"
	@printf "  %-24s %s\n" "traces-python-replay" "Replay Python traces against committed R sampling traces"
	@printf "  %-24s %s\n" "fixtures-synthetic-edge" "Generate tests/fixtures/synthetic_adaptive_sampling_edge"
	@printf "  %-24s %s\n\n" "fixtures-all" "Run the main fixture generators and trace exporters"
	@printf "Common overrides:\n\n"
	@printf "  make traces-r TRACE_KINASES=PRKAA1,MAPK1 TRACE_TOP_N=10\n"
	@printf "  make traces-python TRACE_KINASES=PRKAA1,MAPK1 SVM_MODE=r_parity\n"
	@printf "  make fixtures-r-l6 RSCRIPT=/path/to/Rscript\n"
	@printf "  make fixtures-fragile PYTHON=python3\n"

fixtures-r-small:
	$(RSCRIPT) scripts/generate_r_fixtures.R --outdir $(R_SMALL_OUTDIR)

fixtures-r-l6:
	$(RSCRIPT) scripts/generate_r_l6_fixtures.R --outdir $(R_L6_OUTDIR)

traces-r:
	$(RSCRIPT) scripts/generate_r_l6_fixtures.R \
		--outdir $(R_L6_OUTDIR) \
		--trace_kinases $(TRACE_KINASES) \
		--trace_top_n $(TRACE_TOP_N)

fixtures-fragile:
	$(PYTHON) scripts/generate_fragile_support_reference.py

traces-python:
	$(PYTHON) scripts/export_python_prediction_traces.py \
		--trace-kinases $(TRACE_KINASES) \
		--svm-mode $(SVM_MODE) \
		--debug-top-n $(TRACE_TOP_N) \
		--outdir $(PYTHON_TRACE_OUTDIR)

traces-python-replay:
	$(PYTHON) scripts/export_python_prediction_traces.py \
		--trace-kinases $(TRACE_KINASES) \
		--svm-mode $(SVM_MODE) \
		--debug-top-n $(TRACE_TOP_N) \
		--sampling-trace-dir $(R_L6_OUTDIR)/prediction_trace \
		--outdir $(PYTHON_TRACE_OUTDIR)

fixtures-synthetic-edge:
	$(PYTHON) scripts/generate_synthetic_adaptive_sampling_edge_fixtures.py

fixtures-all: fixtures-r-small fixtures-r-l6 traces-r fixtures-fragile traces-python fixtures-synthetic-edge

PYTHON ?= python
RSCRIPT ?= Rscript
TRACE_KINASES ?= PRKAA1,MAPK1
TRACE_TOP_N ?= 10
PYTHON_TRACE_OUTDIR ?= tests/fixtures/python_reference_l6/prediction_trace
R_L6_OUTDIR ?= tests/fixtures/r_reference_l6
R_SMALL_OUTDIR ?= tests/fixtures/r_reference
FRAGILE_OUTDIR ?= tests/fixtures/fragile_support_reference
L6_STRESS_OUTDIR ?= tests/fixtures/r_reference_l6_seam_stress
SYNTHETIC_EDGE_OUTDIR ?= tests/fixtures/synthetic_adaptive_sampling_edge

.PHONY: help fixtures-r-small fixtures-r-l6 traces-r fixtures-fragile fixtures-r-l6-seam-stress traces-python traces-python-replay fixtures-synthetic-edge fixtures-all test-seams

help:
	@printf '%s\n' \
	  'Available targets:' \
	  '  make fixtures-r-small              Generate the small R-backed fixture family' \
	  '  make fixtures-r-l6                 Generate the main L6 R-backed fixture family' \
	  '  make traces-r                      Regenerate the committed R L6 prediction trace' \
	  '  make fixtures-fragile              Generate the curated fragile-support seam fixture' \
	  '  make fixtures-r-l6-seam-stress     Generate the smaller R-backed L6 seam-stress fixture' \
	  '  make traces-python                 Export Python prediction traces' \
	  '  make traces-python-replay          Export Python prediction traces replaying R sampling rows' \
	  '  make fixtures-synthetic-edge       Generate the synthetic adaptive-sampling edge fixture' \
	  '  make test-seams                    Run the seam-focused parity tests' \
	  '  make fixtures-all                  Generate the fixture families used by seam tests'

fixtures-r-small:
	$(RSCRIPT) scripts/generate_r_fixtures.R --outdir $(R_SMALL_OUTDIR)

fixtures-r-l6:
	$(RSCRIPT) scripts/generate_r_l6_fixtures.R --outdir $(R_L6_OUTDIR)

traces-r:
	$(RSCRIPT) scripts/generate_r_l6_fixtures.R --outdir $(R_L6_OUTDIR) --trace_kinases $(TRACE_KINASES) --trace_top_n $(TRACE_TOP_N)

fixtures-fragile:
	$(PYTHON) scripts/generate_fragile_support_reference.py

fixtures-r-l6-seam-stress:
	$(PYTHON) scripts/generate_l6_seam_stress_reference.py --outdir $(L6_STRESS_OUTDIR)

traces-python:
	$(PYTHON) scripts/export_python_prediction_traces.py \
		--trace-kinases $(TRACE_KINASES) \
		--svm-mode r_parity \
		--debug-top-n $(TRACE_TOP_N) \
		--outdir $(PYTHON_TRACE_OUTDIR)

traces-python-replay:
	$(PYTHON) scripts/export_python_prediction_traces.py \
		--trace-kinases $(TRACE_KINASES) \
		--svm-mode r_parity \
		--debug-top-n $(TRACE_TOP_N) \
		--sampling-trace-dir tests/fixtures/r_reference_l6/prediction_trace \
		--outdir $(PYTHON_TRACE_OUTDIR)

fixtures-synthetic-edge:
	$(PYTHON) scripts/generate_synthetic_adaptive_sampling_edge_fixtures.py

test-seams:
	pytest -q tests/test_reference_workflow_seams.py tests/test_fragile_support_reference.py

fixtures-all: fixtures-r-small fixtures-r-l6 fixtures-fragile fixtures-r-l6-seam-stress fixtures-synthetic-edge

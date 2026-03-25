# L6 seam-stress reference dataset

This directory is a smaller R-backed seam fixture derived by slicing the committed L6 reference outputs.

It is intentionally not a second independent biological source dataset. Instead, it keeps the reference provenance R-backed while narrowing the row set to stress different native seam decisions:
- thinner candidate pools for selected kinases under stricter candidate-selection settings
- retained exact R sampling replay for a smaller traced kinase subset
- full overlap-kinase weight behaviour preserved by keeping the committed L6 score columns intact

Row count: 431
Kinase count: 28
Trace kinases: MAPK1, IRAK1
Candidate-selection settings: top=50, score_threshold=0.8, inclusion=20
Trace replay settings: top=30, score_threshold=0.6, inclusion=5, ensemble_size=10, n_iterations=5, random_state=1, debug_top_n=10

Provenance:
- profile_scores.csv / motif_scores.csv / combined_scores.csv are direct row slices of the committed L6 R reference outputs
- motif_sizes.csv / profile_sizes.csv / combined_weights.csv are direct L6 R-backed seam metadata tables
- predMat.csv and prediction_top30.csv are R-backed prediction outputs sliced to the seam-stress row set
- prediction_trace/* is filtered from the committed L6 R trace for the traced kinases only

Files:
- profile_scores.csv / motif_scores.csv / combined_scores.csv: seam-score reference tables
- motif_sizes.csv / profile_sizes.csv / combined_weights.csv: score-combination weight inputs and outputs
- candidate_substrates.csv / screening_summary.csv: stricter candidate-selection seam references
- predMat.csv / prediction_top30.csv: R-backed prediction ranking references on the seam-stress row set
- prediction_trace/*: filtered R sampling and debug-trace reference tables for replay checks

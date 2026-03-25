# Fragile-support reference dataset

This dataset is a curated L6-derived reference family designed to stress decision fragility rather than broad coverage.

It is intentionally smaller and more uneven than the full L6 reference family:
- mixed kinase support counts
- uneven motif/profile agreement
- smaller and more fragile candidate pools
- at least one dropped kinase, one kinase below inclusion, one just above inclusion, and multiple robust kinases

Selected kinases: MAPK1, AKT1, IRAK1, PRKAA1, PRKAA2, LCK
Row count: 105
Candidate selection settings: top=50, score_threshold=0.8, inclusion=20

This directory is not a blanket parity claim. It is a committed curated dataset for widening evidence beyond the main L6 path and for future seam expansion.

Files:
- phospho_matrix.csv: curated phosphosite matrix
- site_sequences.csv: centralized sequences for the curated site index
- substrate_map.csv: selected kinase-to-site mapping
- motif_sequences.csv: flattened kinase motif sequences for the selected substrate sites
- profile_matrix.csv / profile_sizes.csv: Python-built substrate profile reference outputs
- profile_scores.csv / motif_scores.csv: deterministic scoring seam outputs
- motif_sizes.csv / combined_scores.csv / combined_weights.csv: score-combination seam outputs
- candidate_substrates.csv: candidate selection output under the configured threshold settings
- screening_summary.csv: selection-summary table for the curated dataset

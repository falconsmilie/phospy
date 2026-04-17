# Synthetic adaptive-sampling decision seam fixtures

These fixtures are intentionally small and deterministic.
They complement the R-backed L6 parity traces by pinning edge-case replay behaviour:

- tied candidate scores with stable mergesort ordering
- tiny candidate and negative-pool sizes
- explicit per-iteration sampling overrides
- exact final top-site decisions on a deterministic replay path

They are not a standalone claim of PhosR parity.

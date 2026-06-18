# Test Audit Rubric

Use this rubric to classify each test file during the PhosPy test audit.

Core audit question:

> What bug, scientific invariant, public contract, or reproducibility guarantee would become unsafe if this test file disappeared?

If the answer is "none" or "unclear," the file should usually be rewritten, consolidated, or deleted.

Increasing coverage alone is not a sufficient reason to add tests.

## Categories

### Scientific Invariant

- definition: Verifies a scientific rule that must always hold (for example, monotonicity, conservation, boundedness, or deterministic invariants under fixed inputs/seeds).
- default action: keep.
- examples: Invariant checks for scoring behavior, normalization constraints, and stable ranking semantics for known scientific conditions.
- deletion/rewrite guidance: Delete only if the invariant is no longer part of supported behavior; otherwise rewrite to assert the invariant directly instead of incidental implementation details.

### Public API Contract

- definition: Verifies behavior promised to users at public boundaries (function/class signatures, return schema, error contracts, documented flags, and documented workflow outputs).
- default action: keep.
- examples: Tests asserting documented parameters, output column presence, stable exception types/messages where documented.
- deletion/rewrite guidance: Rewrite when tests pin undocumented internals; delete only when the corresponding public contract is intentionally removed and docs are updated.

### Regression for Past Bug

- definition: Encodes a previously observed failure so it cannot silently return.
- default action: keep.
- examples: A test added after a null-handling crash, shape mismatch, seed nondeterminism, or edge-case overflow bug.
- deletion/rewrite guidance: Keep unless fully subsumed by a broader invariant/property test that demonstrably covers the failure mode; if rewritten, preserve the original failing scenario as data.

### Provenance/Golden Contract

- definition: Verifies reproducibility-critical output/provenance artifacts against committed golden references.
- default action: keep.
- examples: Provenance hash checks, workflow metadata contract checks, and stable fixture-backed output comparisons.
- deletion/rewrite guidance: Rewrite when golden fixtures include incidental/non-contract fields; delete only when reproducibility policy explicitly changes and replacement controls are introduced.

### Parity

- definition: Compares Python rewrite behavior to the accepted R/PhosR reference behavior for defined seams.
- default action: keep.
- examples: Reference fixture comparisons for prediction or workflow stages designated as parity-bearing.
- deletion/rewrite guidance: Rewrite when tolerance/fixture strategy is brittle or duplicates stronger contract coverage; delete only for intentionally retired parity seams with explicit sign-off.

### Performance/Memory Contract

- definition: Verifies supported runtime and memory envelopes for key workflows or kernels.
- default action: keep.
- examples: Tests guarding algorithmic regressions in large-input clustering, scoring, or workflow execution paths.
- deletion/rewrite guidance: Rewrite environment-fragile checks into robust contract thresholds; delete only when the feature/scale contract itself is removed from supported behavior.

### Property-Based Invariant

- definition: Uses generated inputs to verify broad invariants across many cases.
- default action: keep.
- examples: Hypothesis tests validating idempotence, shape invariants, domain constraints, and permutation-insensitive properties where expected.
- deletion/rewrite guidance: Rewrite if flakey, over-constrained, or asserting internals; delete only when the protected invariant is obsolete or fully covered elsewhere with equal breadth.

### Internal Implementation Detail

- definition: Asserts private helper structure, exact call paths, internal naming, or incidental ordering that is not part of contract behavior.
- default action: rewrite or delete.
- examples: Tests that fail on harmless refactors (mock call-count pinning without contract meaning, private function name coupling).
- deletion/rewrite guidance: Rewrite against public behavior or true scientific invariants; delete if no meaningful external or scientific guarantee is being protected.

### Duplicate Validation Matrix

- definition: Repeats the same assertions across many files/parameter combinations without adding unique risk coverage.
- default action: consolidate.
- examples: Near-identical test matrices across modules validating the same schema/validator branch.
- deletion/rewrite guidance: Consolidate into parameterized/shared tests with one authoritative matrix; delete exact duplicates after equivalence is confirmed.

### Obsolete Behavior

- definition: Protects behavior no longer supported, documented, or desired.
- default action: delete.
- examples: Tests for removed flags, deprecated data shapes, or legacy compatibility paths no longer promised.
- deletion/rewrite guidance: Delete together with any stale fixtures/helpers; if partial behavior remains supported, rewrite narrowly to the current contract.

## Action-Oriented Examples

- keep: A parity or provenance test that protects reproducibility-critical outputs and catches scientifically meaningful drift.
- rewrite: A test that currently asserts private call order but should instead assert documented output behavior.
- consolidate: Multiple files that validate the same input schema matrix can become one parameterized contract test.
- delete: A test that only checks removed/undocumented legacy behavior and has no current contract value.

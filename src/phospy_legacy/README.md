# phospy_legacy

This package is the pre-cutover implementation kept as migration reference
material only.

- It is not a supported public API target.
- It must not receive new architectural or feature work.
- All new design and implementation must land under `src/phospy/`.

When logic is needed from this tree, copy it deliberately into the new package
with explicit review rather than extending legacy modules in place.

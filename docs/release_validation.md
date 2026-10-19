# Release validation

PhosPy validates supported releases on Ubuntu, Windows, and macOS with Python
3.11 and 3.12. The permanent source lane runs repository hygiene, Ruff, strict
Pyright coverage, the generated architecture inventory check, the complete
non-network tests, public-consumer contracts, scientific parity, the non-wall-clock
performance selection, reference-bundle validation, release and golden tests, and
a strict MkDocs build. Absolute wall-clock contracts run in a separate fresh-process
matrix across the same six operating-system/Python lanes with the same constrained
dependency profile. Together the source/runtime and isolated selectors partition
the complete performance suite without dropping or duplicating a node within a
lane.

Run the source lane locally with:

```console
make architecture-release-source-check
```

The source lane leaves enough process headroom for the two fresh-process L6
evidence runs and their exact repeatability comparison after the functional
selection completes.

The complete release check also builds one wheel and one sdist, validates their
metadata with Twine, compares packaged production files byte-for-byte with the
committed source blobs, and installs each distribution outside the checkout for
consumer verification:

```console
make architecture-release-check
```

`make release-check` is the maintained alias for the complete check. Network
acquisition and all six operating-system/Python combinations remain required in
hosted CI even when a local environment cannot acquire dependencies.

The exact batch-correction provenance-byte regression has a separate
authoritative hosted lane. Its existing test predicate requires Python
`3.11.15`, NumPy `2.2.6`, pandas `2.2.3`, SciPy `1.15.2`, and scikit-learn
`1.6.1`, Windows AMD64, Anaconda Python, and the OpenBLAS architecture recorded
for the immutable capture.
The pinned bytes were produced on Windows, and hosted execution with the exact
runtime on Linux does not reproduce them. The passing capture environment uses
the Anaconda Windows x64 Python build `3.11.15-hb00fc5c_1` and NumPy's OpenBLAS
`SkylakeX` computation target. CI recreates that environment with Miniconda and
checks the declared target against OpenBLAS's loaded runtime architecture rather
than trusting the environment-variable label alone. The fixture retains its
pre-regeneration Git blob identity. CI uploads a JUnit report for the exact node
so a skip cannot be mistaken for execution in the authoritative lane.

Disposable local output belongs under `build/` or another ignored path.
Verification logs and downloadable workflow artifacts are retained by GitHub
Actions; they are not committed to the repository. The repository hygiene gate
rejects tracked migration-process material and tracked verification artifacts.

## Accepted Wave 10 gate

Wave 10 was accepted on 19 October 2026 at commit
`a4b9929fb41f6507cab3011dafc1012e44c5ce19`, tree
`ecaba2700d978f31e170f02efe38d2a83f11ee6d`, with sole parent
`f2c00c9290187dd6b73a97114240808eb0f3abc0`. Authoritative hosted run
[`39161306082`](https://github.com/falconsmilie/phospy/actions/runs/39161306082)
completed with 43 of 43 jobs successful. Its required operating-system and
Python matrix comprised six source/runtime jobs, six isolated wall-clock
performance jobs, and six installed-distribution jobs. All twelve wheel and
source-distribution consumers passed, together with the minimum-dependency,
release/golden, typing, Ruff, documentation, audit, fixture, reference-bundle,
hygiene, and packaging gates.

The accepted exact-byte provenance lane recorded the explicit disposition
`ineligible_missing_required_cpu_features`: the loaded runtime did not expose
`AVX512F`, `AVX512CD`, or `AVX512_SKX`. Portable semantic provenance, stable
serialization, and both fresh-process repeatability checks still executed and
passed. The supplemental exact-byte predicate did not execute and was not
reported as a pass.

The accepted macOS Python 3.12 SPS measurement was 18.955900 seconds under the
unchanged 20-second ceiling. Its narrow margin is an accepted qualification;
the workload, warm-up, iteration count, and ceiling were not relaxed.

This gate records architecture-migration acceptance only. It does not create
or authorise a tag, package publication, deployment, or release publication.

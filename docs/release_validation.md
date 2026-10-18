# Release validation

PhosPy validates supported releases on Ubuntu, Windows, and macOS with Python
3.11 and 3.12. The permanent source lane runs repository hygiene, Ruff, strict
Pyright coverage, the generated architecture inventory check, the complete
non-network tests, public-consumer contracts, scientific parity and performance
contracts, reference-bundle validation, release and golden tests, and a strict
MkDocs build.

Run the source lane locally with:

```console
make architecture-release-source-check
```

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
`1.6.1`; the predicate itself does not impose an operating-system condition.
The pinned bytes were produced on Windows, and hosted execution with the exact
runtime on Linux does not reproduce them. The passing capture environment uses
the Anaconda Windows x64 Python build `3.11.15-hb00fc5c_1` and NumPy's OpenBLAS
`SkylakeX` computation target. CI recreates that environment with Miniconda,
checks the platform, interpreter build, and dependency versions before pytest,
and uploads a JUnit report for the exact node so a skip cannot be mistaken for
execution in the authoritative lane.

Disposable local output belongs under `build/` or another ignored path.
Verification logs and downloadable workflow artifacts are retained by GitHub
Actions; they are not committed to the repository. The repository hygiene gate
rejects tracked migration-process material and tracked verification artifacts.

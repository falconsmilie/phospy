# Contributing

Thank you for your interest in PhosPy.

PhosPy 1.0.0 is intentionally narrow. The most helpful contributions keep that scope clear: small changes, good tests,
and straightforward review.

## Get Set Up

Install the package with test and development tooling:

```bash
pip install -e ".[test,dev]"
```

Install the Git hooks:

```bash
pre-commit install
```

Run the usual local checks:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

## Code Style

PhosPy uses Ruff for linting and formatting.

- `ruff check --fix` handles linting, import sorting, and safe auto-fixes
- `ruff format` handles formatting
- `pre-commit` runs the repository checks before commit

The current repository rules live in:

- `.pre-commit-config.yaml`
- `pyproject.toml`

## Tests

PhosPy has two main test layers:

- **non-parity tests** for package behaviour that should stay stable without R
- **parity tests** for fixture-backed comparisons against committed R/PhosR reference outputs

Run the parity layer with:

```bash
pytest -m parity
```

If you want the parity suite to print its optional comparison summaries while you are debugging, these environment
variables are available:

- `PHOSPY_SHOW_PARITY`: master switch for parity metrics output
- `PHOSPY_SHOW_PROFILE_CONSTRUCTION`: adds the optional profile-construction summary
- `PHOSPY_SHOW_PREDICTION_MODE_COMPARISON`: adds default-versus-`r_parity` prediction comparison metrics
- `PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON`: adds replayed prediction comparison metrics

The more specific flags only take effect when `PHOSPY_SHOW_PARITY` is also enabled. Run pytest with `-s` if you want
to see the printed summaries in the terminal.

Quick Linux or macOS example:

```bash
PHOSPY_SHOW_PARITY=1 PHOSPY_SHOW_PROFILE_CONSTRUCTION=1 PHOSPY_SHOW_PREDICTION_MODE_COMPARISON=1 PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON=1 pytest -m parity -s
```

For the full cross-platform command set, sample output, and notes on how to read the printed metrics, see
[`docs/parity.md`](docs/parity.md).

## Release Gate

The practical 1.0.0 release gate is:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

The non-parity suite includes the documented example smoke workflow. The parity suite covers the committed R-backed
seams. The `make test-parity` shortcut enables `PHOSPY_SHOW_PARITY=1` and
`PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON=1`, so it prints the core parity summaries plus the replayed
mode-comparison block rather than every optional metrics block.

## CI Expectations

GitHub Actions currently runs three checks:

- `pre-commit`
- non-parity tests
- parity tests against the committed fixture snapshots

If you intentionally change fixture-producing behaviour, regenerate the affected fixtures, commit them, and explain the
change clearly in the pull request.

## Good First Directions

The public roadmap lives in [`docs/roadmap.md`](docs/roadmap.md). The most useful contributions are the ones that
extend the current supported surface without over-claiming: native-workflow CLI coverage, better seam-level validation,
clearer diagnostics, and carefully scoped PhosR-inspired additions.

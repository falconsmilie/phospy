# Contributing to PhosPy

Thank you for helping improve PhosPy. Small, focused changes are easier to
review and safer to release than broad rewrites that mix unrelated concerns.

## Before You Open a Pull Request

Please make sure the change is:

- accurate for the current source and tests;
- clear to scientists and first-time users;
- limited to one understandable purpose;
- honest about supported, advanced, and internal behavior; and
- covered by the relevant tests or documentation checks.

## Documentation Style

PhosPy uses an APA-informed style adapted for technical documentation. The goal
is readable scientific writing, not manuscript formatting.

### Voice and Structure

- Write in a warm, direct, professional voice.
- Lead with the user’s task or scientific question.
- Prefer active voice and short paragraphs.
- Use **you** for instructions and **PhosPy** for software behavior.
- Put the working example before advanced implementation detail.
- Keep one authoritative user page for each workflow.
- Define an abbreviation the first time it appears on a page.
- Use inclusive, specific language; singular **they** is acceptable.

### Headings and Terms

- Use title case for page titles and headings.
- Keep heading levels in order; do not skip from level 2 to level 4.
- Preserve public class names, parameters, enum values, table columns, and file
  paths exactly as code.
- Use established scientific terminology, then explain how PhosPy applies it.
- Avoid architecture language in beginner guides unless users must act on it.

### Numbers and Statistics

- Use numerals for measurements, versions, thresholds, sample counts, and code
  values.
- Use a leading zero for ordinary decimal values below 1. For statistical
  quantities that cannot exceed 1, follow APA convention in prose, such as
  *p* = .03; preserve exact field names and machine-readable values in code.
- Italicise statistical symbols in prose, including *p*, *t*, *F*, and *r*.
  Do not italicise code identifiers such as `P.Value` or `adj.P.Val`.
- State effect sizes, uncertainty, assumptions, and limitations where they
  matter; do not equate statistical significance with biological importance.

### Citations and References

Use author–date citations for external scientific claims, for example,
`(Author & Author, 2024)`. Add a **References** section to the page when it cites
external literature. Keep entries alphabetical, use sentence case for article
and book titles, and include a DOI or stable source link when available.

Do not add citations merely to describe PhosPy’s API. Current public exports,
source definitions, tests, and release contracts are the source of truth for
software behavior.

### Documentation Locations

- Keep welcoming and task-based guides in `docs/`.
- Keep each complete workflow contract in its page under `docs/api/`.
- Keep architecture decisions in `docs/adr/`.
- Keep testing-audit material in `docs/testing/`.
- Use `phospy.api` for public request, configuration, result, enum, reference,
  and error imports unless a documented route says otherwise.
- Use `phospy.advanced` only for supported advanced configuration.
- Never present private validators or internal execution modules as public API.

## Local Setup

```bash
pip install -e ".[dev]"
pip install -e ".[dev,parquet]"  # Optional Parquet support
```

## Tests to Run

For most changes:

```bash
pytest -m "not parity"
```

For public documentation and examples:

```bash
pytest tests/unit/test_public_contract_import_routes.py
pytest tests/unit/test_public_examples_contract.py
pytest tests/integration/test_public_examples_smoke.py
python -m mkdocs build --strict
```

Run parity tests when scientific logic or fixture-backed behavior changes:

```bash
pytest tests/parity -m "parity and not parity_diagnostic" -s
```

## Release Checks

Normal contributor checks are the focused commands above: the default
non-parity suite for most changes, documentation checks for public docs, and
parity checks when scientific logic or fixture-backed behavior changes.

Run the aggregate release command before publishing, and also after changes to
scientific, provenance, performance, distribution, reference-bundle, or
public-contract behavior:

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet,docs]"
make release-check
```

The detailed procedure lives in [Maintenance](maintenance.md). Ordinary
CI/build success and source-tree tests provide normal development confidence;
they are not sufficient for publishing. Final release verification is stricter
because it requires a Git-backed checkout for staged-byte validation and
verifies the freshly built wheel and sdist. A successful source-tree test run
alone does not prove the built distributions are valid.

The optional 50,000 x 48 scale benchmark is informational and
machine-dependent. Run it locally with `make benchmark-release-scale`; it is not
part of the release gate or continuous integration.

## Code Style and Type Checking

Use Ruff for linting and formatting:

```bash
ruff check --fix
ruff format
```

Run the configured strict-type coverage and Pyright checks:

```bash
python tools/testing/pyright_strict_coverage.py --check
python scripts/run_pyright.py
```

Prefer precise boundary types over `Any`. Keep suppressions narrow and local.
When Pyright cannot model correct runtime behavior, use the required inline
format with explicit rule names and a concrete technical reason:

```python
# pyright: ignore[reportRuleName] - Concrete reason the runtime operation is safe
```

Blanket ignores, placeholder reasons, and file-wide strictness downgrades are
rejected by the strict-coverage check.

For the architecture behind these boundaries, see
[ADR 0002: Internal Workflow Architecture for PhosPy](adr/adr_0002_internal_workflow_architecture.md).

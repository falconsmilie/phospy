# Documentation Audit and Restructure Plan

## Executive Summary

PhosPy documentation contained strong advanced material but lacked a clear beginner entry flow. The main issue was information architecture, not technical accuracy. This restructure introduces progressive disclosure:

- Layer 1 (beginner): guided onboarding path from "What is PhosPy?" to first workflow and core concepts.
- Layer 2 (advanced): preserved detailed API/validation/science/architecture/maintainer content grouped into explicit categories.

Advanced pages were preserved and re-linked rather than flattened.

## Documentation Audit

### Beginner-facing (new)

- `docs/index.md`
- `docs/getting-started/index.md`
- `docs/getting-started/what-is-phospy.md`
- `docs/getting-started/quickstart-first-workflow.md`
- `docs/concepts/core-concepts.md`
- `docs/learning-paths/choose-your-path.md`

### Advanced/reference (existing, preserved)

- `docs/api.md`
- `docs/cli.md`
- `docs/output_bundles.md`
- `docs/validation.md`
- `docs/frame_ownership.md`

### Scientific/parity/governance (existing, preserved)

- `docs/parity.md`
- `docs/roadmap.md`
- `docs/architecture/legacy_science_gap_audit.md`
- `docs/architecture/activity_science_port_review.md`

### Architecture/ADRs (existing, preserved)

- `docs/architecture/*.md`
- `docs/adr/*.md`

### Contributor/maintainer-oriented (existing + grouped)

- `docs/fixtures.md`
- `docs/architecture/rewrite_cutover_boundary.md`
- `docs/architecture/phospy_architecture_reset_notes.md`
- governance docs in parity/architecture

### Entry-point information previously buried

- High-value onboarding information was primarily in `README.md`, `docs/api.md`, and dense contract pages.
- There was no docs homepage, no beginner path, and no category index pages.

### Navigation pain points

- Flat docs root with many advanced pages.
- No explicit separation between "how to start" and "deep contract/governance."
- ADR and architecture materials required prior context to navigate efficiently.

### Mixed-audience pages identified

- `README.md` (beginner onboarding mixed with dense contract/governance details)
- `docs/api.md` (includes quick-usage onboarding plus deep contract inventory)
- `docs/validation.md` (mixes practical troubleshooting with deep boundary internals)
- `docs/parity.md` (important context for advanced users; too dense for first read)

## Proposed Documentation Map / Hierarchy

- Getting started
  - `docs/getting-started/index.md`
  - `docs/getting-started/what-is-phospy.md`
  - `docs/getting-started/quickstart-first-workflow.md`
- Tutorials/user guides
  - `docs/user-guides/index.md`
  - `docs/cli.md`
  - `docs/output_bundles.md`
- Concepts
  - `docs/concepts/core-concepts.md`
- Workflow guides
  - `docs/workflow-guides/index.md`
- API/reference
  - `docs/reference/index.md`
  - `docs/api.md`
- Validation/contracts
  - `docs/contracts/index.md`
  - `docs/validation.md`
  - `docs/frame_ownership.md`
- Architecture/ADRs
  - `docs/architecture/index.md`
  - `docs/adr/index.md`
- Scientific/parity/governance
  - `docs/science/index.md`
  - `docs/parity.md`
  - `docs/roadmap.md`
  - `docs/architecture/legacy_science_gap_audit.md`
- Contributor/maintainer
  - `docs/contributor/index.md`
  - `docs/fixtures.md`
  - architecture/governance links above

## Beginner Journey

Primary path:

1. `docs/index.md`
2. `docs/getting-started/what-is-phospy.md`
3. `docs/getting-started/quickstart-first-workflow.md`
4. `docs/concepts/core-concepts.md`
5. `docs/learning-paths/choose-your-path.md`

Outcome: a first-time user can understand product shape and run the supported happy path before diving into architecture/parity governance details.

## Advanced/Reference Journey

Primary path:

1. `docs/reference/index.md`
2. `docs/api.md` and `docs/validation.md`
3. `docs/science/index.md` and `docs/architecture/index.md`
4. `docs/adr/index.md` for design rationale and historical decisions

Outcome: advanced users can reach deep contracts, scientific policy, and architecture intentionally, without being routed through beginner material.

## Linking and Navigation Recommendations

- Beginner pages should link outward to contract/reference pages at decision points.
- Advanced pages should include "if you are new, start here" links back to beginner pages.
- Category index pages should include explicit audience and purpose.
- Dense docs should include onward links so no page ends as a dead end.

Implemented via:

- `mkdocs.yml` navigation hierarchy
- new index hubs under `docs/` categories
- cross-link additions in dense advanced pages

## Page-by-Page Rewrite / Move / Split Recommendations

- `README.md`: rewritten as beginner-oriented project entry with links to layered docs.
- `docs/api.md`: keep advanced depth; add audience/context and beginner back-links.
- `docs/validation.md`: keep deep contract detail; add context links and route from/to concepts.
- `docs/parity.md`: keep governance depth; position as advanced science/governance page.
- `docs/cli.md`, `docs/output_bundles.md`, `docs/fixtures.md`, `docs/frame_ownership.md`, `docs/roadmap.md`: keep content, add audience and onward links.
- ADR and architecture pages: preserved; linked through new `docs/adr/index.md` and `docs/architecture/index.md`.

## Sample Rewritten Docs Produced

- Homepage: `docs/index.md`
- Beginner introduction: `docs/getting-started/what-is-phospy.md`
- Quickstart: `docs/getting-started/quickstart-first-workflow.md`
- Navigation/index pages:
  - `docs/getting-started/index.md`
  - `docs/user-guides/index.md`
  - `docs/workflow-guides/index.md`
  - `docs/reference/index.md`
  - `docs/contracts/index.md`
  - `docs/science/index.md`
  - `docs/architecture/index.md`
  - `docs/adr/index.md`
  - `docs/contributor/index.md`

## Risks, Tradeoffs, and Open Questions

- If docs are rendered outside MkDocs, `mkdocs.yml` nav will not apply; Markdown index pages still provide navigability on GitHub.
- Existing direct links to advanced pages remain valid, but users may bypass onboarding unless docs landing is highlighted in README and site root.
- As new docs are added, category placement and cross-link rules should be maintained to avoid drift back to mixed-audience pages.

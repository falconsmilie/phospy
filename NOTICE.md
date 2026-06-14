This repository is an unofficial Python implementation of selected functionality
from the PhosR project.

Upstream project:

- Name: PhosR
- Upstream repository: https://github.com/PYangLab/PhosR
- Upstream license: GPL-3.0

Attribution:

- Scientific credit for the original methods, package design, and workflow belongs to the PhosR authors and maintainers.
- This repository should cite and acknowledge the original PhosR publications and repository.

Licensing position:

- This repository is distributed under GPL-3.0.
- Copyright in the original PhosR code remains with its original authors.
- Copyright in changes made in this repository belongs to the contributors to this repository.

Status:

- This is an independent, unofficial implementation.
- It is not affiliated with or endorsed by the PhosR authors unless they explicitly say so.

Packaged reference data:

- Runtime reference bundles are package data only when committed under
  `src/phospy/data/reference_bundles`.
- The current packaged rat lane is `rat/l6_native`. Its manifest records a
  PhosR-derived lineage from `phospho.L6.ratio.pe`, `PhosphoSite.mouse`, and
  `motif.mouse.list`, generated with
  `scripts/active/generate_r_l6_fixtures.R` and packaged on 2026-04-16.
- PhosR package metadata declares GPL-3 + file `LICENSE`, but upstream PhosR
  documentation identifies PhosphoSitePlus-derived annotations and
  PRIDE-sourced phosphoproteomics data in the objects used for the rat lane.
  This repository does not record independent written redistribution approval
  for the exact derived rat CSV snapshot.
- This release does not include packaged human or mouse reference lanes.
- Human or mouse packaged lanes require manifest metadata documenting source
  provenance, license, retrieval method, redistribution status, redistribution
  basis, limitations, and supported uses before they can be bundled.
- Restricted third-party scientific reference datasets, including
  PhosphoSitePlus or Kinase Library data, must not be copied into this
  repository unless redistribution is explicitly permitted by license or written
  permission.

# Notice

PhosPy is an independent Python implementation of selected phosphoproteomics
workflow ideas inspired by the PhosR project. It is not affiliated with or
endorsed by the PhosR authors unless they explicitly say so.

## Upstream Project

- **Name:** PhosR
- **Repository:** https://github.com/PYangLab/PhosR
- **Licence:** GPL-3 + file LICENSE

## Attribution

Scientific credit for the original methods, package design, and workflow ideas
belongs to the PhosR authors and maintainers. Users should cite the original
PhosR publications, the PhosR repository, and the PhosPy software release when
PhosPy is used in scientific work.

## Licensing Position

- PhosPy is distributed under GPL-3.0.
- Copyright in the original PhosR code remains with its original authors.
- Copyright in changes made in this repository belongs to this repository's
  contributors.

## Project Status

PhosPy supports selected, documented workflow lanes. It does not claim full
PhosR package parity, full PhosR API compatibility, or endorsement by upstream
PhosR maintainers.

## Packaged Reference Data

Runtime reference bundles are package data only when they are committed under
`src/phospy/data/reference_bundles`.

The packaged rat reference lane `rat/l6_native` contains derived CSV snapshots
generated from documented data objects in upstream PhosR 1.20.0:

- `phospho.L6.ratio.pe`
- `PhosphoSite.mouse`
- `motif.mouse.list`

The upstream PhosR package is distributed under GPL-3 + file `LICENSE`. PhosPy
is distributed under GPL-3.0 and includes this exact derived snapshot with
attribution to PhosR and the PhosR publications.

This attribution and redistribution statement applies only to the exact files
committed under:

`src/phospy/data/reference_bundles/rat/l6_native`

It does not grant or imply permission to redistribute arbitrary PhosR,
PhosphoSitePlus, PRIDE, Kinase Library, or other third-party datasets outside
this bundled snapshot.

This release does not include packaged human or mouse reference lanes. Human or
mouse packaged lanes require manifest metadata documenting source provenance,
licence, retrieval method, redistribution status, redistribution basis,
limitations, and supported uses before they can be bundled.

Restricted third-party scientific reference datasets, including PhosphoSitePlus
or Kinase Library data, must not be copied into this repository unless
redistribution is explicitly permitted by licence or written permission.

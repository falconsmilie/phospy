# Script Layout

Repository scripts are split by maintenance status:

- `active/`: supported maintainer generators used by default workflows
- `support/`: helper modules used by active generators
- `archive/`: historical parity/debug tooling kept for forensic reference

Default maintainer fixture bootstrap paths are wired through `Makefile`
targets that execute scripts in `active/`.

Archived scripts are intentionally out of the default maintainer lane. See
`archive/README.md` for retained historical tooling and purpose.

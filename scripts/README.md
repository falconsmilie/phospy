# Script Layout

Repository scripts are split by maintenance status:

- `active/`: supported maintainer generators used by default workflows
- `support/`: helper modules used by active generators

Default maintainer fixture bootstrap paths are wired through `Makefile`
targets that execute scripts in `active/`.

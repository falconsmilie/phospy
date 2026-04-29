# Security Policy

## Supported Versions

Security reports should target the current public release line shown in
`pyproject.toml` and `CHANGELOG.md`.

## Reporting a Vulnerability

Please do not open a public issue for a private security concern. Use GitHub's
private vulnerability reporting flow when available, or contact the maintainer
through the repository owner profile.

Helpful reports include:

- affected version or commit
- smallest reproduction you can share safely
- expected and actual behaviour
- impact and any known workaround

## Good Fits for This Channel

- unsafe file handling
- unexpected path traversal
- dependency vulnerability with a concrete impact on PhosPy
- a way for crafted input to execute code or expose private data

## Not a Security Report

Please use normal issues for scientific accuracy questions, documentation gaps,
performance concerns, validation behaviour, or feature requests.

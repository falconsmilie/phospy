# Security Policy

## Supported Versions

PhosPy is a small project, and security fixes are only expected for the latest released version.

Older releases may not receive patches. Users should upgrade to the most recent published version where possible.

| Version | Supported |
|---------|-----------|
| Latest release | Yes |
| Older releases | No |

## Reporting a Vulnerability

If you believe you have found a security vulnerability in PhosPy, please report it privately.

Please do not open a public GitHub issue for suspected security problems.

Report security concerns to:

shane @ rededitor dot net

Please include as much relevant detail as you can, such as:

- a clear description of the issue
- steps to reproduce it
- the affected version or versions
- proof of concept, sample input, or screenshots where appropriate
- the likely impact

## What to Expect

Reports will be reviewed as promptly and carefully as possible.

The maintainers will aim to:

- acknowledge receipt of the report
- investigate whether the issue is reproducible
- assess whether it is a genuine security problem
- prepare a fix or mitigation where appropriate
- publish the fix in a normal project release

Response times may vary depending on maintainer availability and the complexity of the report.

## Coordinated Disclosure

Please allow time for the issue to be investigated and addressed before disclosing it publicly.

Avoid publishing exploit details, proof-of-concept code, or public discussion of the issue until a fix or mitigation
has been released.

## Scope

This policy applies to security issues in the PhosPy codebase and its published packages.

The following are generally out of scope unless they directly create a vulnerability in PhosPy itself:

- problems in third-party platforms or services
- vulnerabilities in dependencies that are not caused by PhosPy's own code or packaging
- installation support requests, usage questions, or feature requests
- theoretical concerns without a clear and reproducible security impact

## Dependencies and Environment

PhosPy depends on third-party Python packages and developer tooling. Reasonable effort will be made to keep those
dependencies current, but no project can guarantee the absence of security issues.

Users should:

- keep PhosPy updated
- review dependency advisories in their own environment
- use isolated environments where appropriate
- apply normal good practice when handling scientific and research data

## Final Note

PhosPy is an open-source scientific software project. It is provided on an "as is" basis, and users are responsible for
deciding whether it is suitable for their own environment and risk profile.

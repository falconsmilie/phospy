# Security Policy

## Supported Versions

PhosPy is a small project. Security fixes are only expected for the latest released version.

| Version | Supported |
| --- | --- |
| Latest release | Yes |
| Older releases | No |

## Reporting a Vulnerability

If you believe you found a security issue in PhosPy, please report it privately.

Do not open a public GitHub issue for suspected security problems.

Report security concerns to:

shane @ rededitor dot net

Please include:

- a clear description of the issue
- steps to reproduce it
- the affected version or versions
- the likely impact
- proof of concept, sample input, or screenshots where useful

## What Counts as a Security Report

Good fits for this channel include issues such as:

- unsafe file handling that could expose data
- code execution paths that should not be reachable
- packaging or distribution problems with security impact
- dependency or configuration behaviour that creates a real vulnerability in PhosPy

The following should use the normal issue tracker instead:

- scientific correctness or parity regressions without a security impact
- documentation bugs
- feature requests
- installation help or general usage questions

## What to Expect

The maintainers will aim to:

- acknowledge receipt of the report
- confirm whether the issue is reproducible
- assess whether it is a genuine security problem
- prepare a fix or mitigation where appropriate
- ship the fix in a normal project release

Response times may vary with maintainer availability and the complexity of the report.

## Scope

This policy applies to security issues in the PhosPy codebase and its published packages.

The following are generally out of scope unless they directly create a vulnerability in PhosPy itself:

- problems in third-party platforms or services
- vulnerabilities in dependencies not caused by PhosPy code or packaging
- theoretical concerns without a clear and reproducible security impact

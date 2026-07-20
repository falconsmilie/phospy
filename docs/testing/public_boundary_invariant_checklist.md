# Public-Boundary Invariant Checklist

PhosPy release candidates must pass the public-boundary integrity gate in both
source tests and isolated installed-artifact verification.

Required invariants:

- Root and `phospy.api` exports must not expose validators, interpreters,
  executors, internal views, ownership-transfer helpers or tokens, private
  constructor switches, warning-suppression switches, validation bypasses, or
  fingerprint bypasses.
- Supported public signatures are inventoried automatically from
  `phospy.__all__` and `phospy.api.__all__`; failures must name the exported
  symbol and parameter.
- Public dataset construction must reject stale provenance when a represented
  table is mutated after provenance creation.
- Public DataFrame-bearing dataset, importer, differential, kinase, signalome,
  and reference result boundaries must isolate caller inputs, nested
  object-dtype cells, and public exports.
- Public JSON-like result and evidence payloads must be recursively immutable
  from both constructor inputs and returned payloads.
- Exported JSON-like dataclass/model fields require either adversarial coverage
  or a written exemption explaining why the field is not nested mutable JSON
  state.
- The installed-artifact verifier must run outside the checkout under
  `python -I`, import only the installed `phospy` distribution and runtime
  dependencies, and reject any loaded `phospy` module from the source checkout.

The stable installed-artifact check is `public-boundary-integrity`. Its required
detail outcomes are:

- `public-signature-boundary`
- `dataset-provenance-binding`
- `public-dataframe-ownership`
- `public-json-immutability`

Final release attestation must fail if the stable check or any required detail
outcome is absent, renamed, skipped, or failed.

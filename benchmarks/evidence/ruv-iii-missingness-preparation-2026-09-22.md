# RUV-III missingness-preparation evidence (2026-09-22)

The [baseline](ruv-iii-missingness-preparation-baseline-2026-09-22.json)
and [final](ruv-iii-missingness-preparation-final-2026-09-22.json) reports were
captured on the same Windows/Python 3.12.10, NumPy 2.5.1, and pandas 3.0.3
environment. Each report records the command, Git revision and dirty state,
tracked-source diff digest, and SHA-256 of every affected production source.
The baseline used a temporary behavior-preserving scalar production version of
the private materialisation and completion stages; its RUV-III executor source
hash is `61e2f885...e6f6e15`. The final vectorised source hash is
`f69cf0e1...fb97a`.

Times are medians in seconds. Speedup is baseline divided by final. Kernel
timing is present in the JSON reports and excluded from preparation timing.

| Tier | Case | Mask baseline -> final | Completion baseline -> final | Complete preparation baseline -> final |
| --- | --- | ---: | ---: | ---: |
| Small | actual missing | 0.012055 -> 0.000157 (76.6x) | 0.136299 -> 0.000439 (310.5x) | 0.149407 -> 0.000657 (227.5x) |
| Small | upstream imputed | 0.012477 -> 0.000158 (79.2x) | n/a | 0.012681 -> 0.000194 (65.2x) |
| Small | mixed | 0.012458 -> 0.000158 (78.6x) | 0.093243 -> 0.000406 (229.9x) | 0.105587 -> 0.000653 (161.7x) |
| Representative | actual missing | 0.924341 -> 0.005310 (174.1x) | 3.958745 -> 0.010690 (370.3x) | 4.879850 -> 0.019571 (249.3x) |
| Representative | upstream imputed | 0.931631 -> 0.005447 (171.0x) | n/a | 0.934997 -> 0.006190 (151.0x) |
| Representative | mixed | 0.935384 -> 0.005336 (175.3x) | 2.387174 -> 0.009830 (242.8x) | 3.305919 -> 0.019442 (170.0x) |
| Stress | actual missing | 11.654747 -> 0.069928 (166.7x) | 43.855400 -> 0.083695 (524.0x) | 54.269054 -> 0.226680 (239.4x) |
| Stress | upstream imputed | 11.770748 -> 0.137317 (85.7x) | n/a | 11.790750 -> 0.107192 (110.0x) |
| Stress | mixed | 11.806960 -> 0.070270 (168.0x) | 23.316754 -> 0.069040 (337.7x) | 34.556357 -> 0.158248 (218.4x) |

Direct validation-stage speedups were 1.8x at small scale, 1.9-2.7x at
representative scale, and 0.9-1.3x at stress scale for governed cases. The
decision gate is met independently by both rules: scalar mask construction and
completion exceed 20% of preparation at representative and stress scale, and
each affected stage improves by much more than 2x.

The one-repeat all-case stress observation showed complete-input preparation at
0.009404 seconds baseline and 0.009046 seconds final. The retained nine-repeat
[baseline](ruv-iii-missingness-preparation-complete-baseline-2026-09-22.json)
and [final](ruv-iii-missingness-preparation-complete-final-2026-09-22.json)
complete-input reports confirm the result: 0.010683 seconds baseline and
0.008914 seconds final, a 1.20x improvement.

Peak process RSS for representative actual/upstream/mixed cases changed from
179.5/177.0/181.4 MiB to 186.9/172.7/180.7 MiB (+4.1%/-2.4%/-0.4%). Stress
actual/upstream/mixed changed from 375.5/353.7/390.1 MiB to
405.3/329.7/407.3 MiB (+7.9%/-6.8%/+4.4%). The actual-missing and mixed
per-call RSS deltas rose because vectorised `nanmedian` uses a dense,
cell-linear scratch workspace instead of repeated row-sized pandas objects. The
largest absolute process-peak increase was 29.8 MiB (7.9%) for the 3,000,000
cell stress matrix. This bounded memory tradeoff is accepted for the measured
218-239x complete-preparation speedups; no new retained data structure or
super-linear memory growth was introduced.

Exact scalar-reference preparation tests and full executor equivalence tests
cover prepared and corrected matrices, governed and actual masks, ordered
coordinates, restored positions, output masks, statuses, diagnostics, warnings,
provenance payloads, and fingerprints. The final implementation retains only
the vectorised production path.

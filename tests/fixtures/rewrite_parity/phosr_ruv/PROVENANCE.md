# PhosR SPS and ruv RUV-III external reference evidence

Generated: 2026-09-17T00:00:00Z
R: R version 4.5.2 (2025-10-31 ucrt)
PhosR: 1.13.1 commit 1be74902b775833c64f5833e70538eaf843cf6a5
ruv: 0.9.7.1

The inputs are deterministic synthetic values authored for PhosPy. No external
biological dataset or bundled PhosR reference dataset is redistributed.

SPS source files contain absolute processed log2 abundance. For every site, the
mean of baseline_1 and baseline_2 was subtracted from all four samples. The
resulting condition-relative log2 files, not the absolute files, were supplied
unchanged to both PhosR getSPS and the PhosPy parity test.

RUV-III consumes the checked-in, already-selected negative controls. It never
reruns SPS. Corrected matrices are compared, while latent factors are not
compared directly because SVD coordinates admit sign/rotation invariances.
External numerical parity applies only to inputs supported by both implementations.
PhosPy rejects singleton replicate sets and a requested k beyond the estimable
latent rank rather than silently reducing it. These are stricter input-validation
contracts, not corrected-matrix parity failures for the supported domain.

Regeneration is separate from CI and requires the exact environment recorded in
REFERENCE_ENVIRONMENT.json. The generator refuses version drift by default.
Normal test execution neither starts R, installs packages, nor accesses a network.
PhosPy has no runtime dependency on R.

Scientific citations:
- Kim HJ et al. (2021), Cell Reports 34(8):108771, doi:10.1016/j.celrep.2021.108771.
- Gagnon-Bartsch J (2019), ruv 0.9.7.1, doi:10.32614/CRAN.package.ruv.

# Data for Violation of S-duality in classical Q-cohomology

Chi-Ming Chang and Ying-Hsuan Lin, [arXiv:2510.24008](https://arxiv.org/abs/2510.24008).
Version 1.0.0; paper release tag `2510.24008-v1.0.0`.

This directory contains the complete symmetry-reduced cohomology count scan
through L=18 for SO(7) and USp(6), four final cluster logs supporting Table I,
two stored Cartan-image bases and portable exports, provenance and file hashes.
The code and validation tools are in the [matching code directory](https://github.com/yinhslin/bps-counting-code/tree/2510.24008-v1.0.0/papers/2510.24008).

- [Table I](data/derived/table_I.csv)
- [Cohomology scan](data/derived/cohomology_scan.csv)
- [Validation results](data/derived/validation.json)
- [Data formats and provenance](DATA_FORMATS.md)
- [Raw-file SHA-256 manifest](data/raw/manifest.json)
- [Complete directory manifest](release_manifest.json)

To validate from `papers/2510.24008/` in the code repository:

```sh
python3 tools/validate_results.py --data-root /path/to/bps-counting/2510.24008/data
```

For a small self-contained download, use `bps_cohomology_code_data.zip` from
the [code release](https://github.com/yinhslin/bps-counting-code/releases/tag/2510.24008-v1.0.0).
It contains both this dataset and the matching code, without older unrelated
repository contents. The combined archive records the exact commits of both
repositories.

All 3,232 included raw result/log files have SHA-256 checksums. The 50,027-file
external inventory lists names and sizes only; the full 3.13 TB of large bases,
intermediates and backups is not included. Large L=18 Q-complex calculations
have not been independently rerun during preparation. The exact rank checks
of the saved Cartan-image matrices do not reconstruct the full restriction map.

The MIT license in this directory applies to this released companion material;
it does not change licensing for files elsewhere in this repository. Please
cite the associated paper and the archived version of this computational companion.

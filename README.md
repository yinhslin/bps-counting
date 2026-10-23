# bps-counting

Counting of 1/16 BPS states and one-loop anomalous dimensions in 3+1d maximally supersymmetric Yang-Mills with gauge group SU(N)

A Mathematica notebook read.nb is provided for reading the data

## BPS counting

Publication: https://doi.org/10.1007/JHEP02(2023)109

Preprint: https://arxiv.org/abs/2209.06728

Progress:
N=2 n=25
N=3 n=19
N=4 n=15

Data is stored in cohomology_N.csv, with entries

n, #derivatives (5 entries), Y, N, #all, #multi-gravitons

## One-loop anomalous dimensions

Preprint: https://arxiv.org/abs/2306.04673

Progress:
N=0 (planar singletrace) n=22
N=2 n=24
N=3 n=17
N=4 n=15

Data is stored in ad_N.csv, with entries (no N for planar)

n, #derivatives (5 entries), Y, N, anomalous dimensions

## Violation of S-duality in classical Q-cohomology (arXiv:2510.24008)

The [2510.24008/](2510.24008/README.md) directory contains the SO(7)/USp(6)
cohomology scan through L=18, Table I job logs, saved Cartan-image bases and
verification records. The [matching code](https://github.com/yinhslin/bps-counting-code/tree/2510.24008-v1.0.0/papers/2510.24008) is in `bps-counting-code`.
Both repositories use the paper release tag `2510.24008-v1.0.0`. Earlier SU(N) datasets
remain in their existing locations.

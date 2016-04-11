# LDpop
### version 0.1.0
LDpop is a program for computing 2-locus likelihoods under the coalescent with recombination. Unlike previous methods, LDpop can correctly account for variable (but piecewise-constant) population size history.

LDpop produces lookup tables that can be used by other programs to estimate recombination maps. Other programs that can use the output of LDpop include:
* [LDhat](https://github.com/auton1/LDhat)
* [LDhot](https://github.com/auton1/LDhot)
* A modified [LDhelmet](https://sourceforge.net/projects/ldhelmet/), to be released soon (the original LDhelmet computes its lookup table internally).

LDpop also provides utilities for efficient posterior sampling of 2-locus ARGs.

## Installation and Dependencies

Prerequisites:
* Scientific distribution of Python 2.7, e.g. [Anaconda](http://continuum.io/downloads), [Enthought Canopy](https://www.enthought.com/products/canopy/)
  * Alternatively, custom installation of pip, the SciPy stack
* Optional: Java 8
  * Not required for computing lookup tables for LDhat/LDhelmet.
  * Required for posterior sampling of 2-locus ARGs.

To install, in the top-level directory of LDpop (where "setup.py" lives), type
```
pip install .
```

## Getting started
Use `run/ldtable.py` to create a lookup table. See
```
run/ldtable.py --help
```
for usage.

By default `run/ldtable.py` uses an exact algorithm to compute the likelihoods.
To use a reasonable approximation that is much faster and scales to larger sample sizes,
use the flag `--approx`.

`run/ldproposal.py` and `run/ImportanceSampler.jar` are for importance sampling from the posterior distribution of 2-locus ARGs.
`run/ldproposal.py` creates a proposal distribution, that `run/ImportanceSampler.jar` uses to sample the ARGs. See their `--help` for instructions.

Also, see the [examples](example/).

## Authors

[John Kamm](mailto:jkamm@stat.berkeley.edu), Jeffrey Spence, Jeffrey Chan, Yun S. Song

## Reference

Kamm, J.A., Spence, J.P., Chan, J., and Song, Y.S. Two-Locus Likelihoods under Variable Population Size and Fine-Scale Recombination Rate Estimation. http://arxiv.org/abs/1510.06017


## License

LDpop is free software under conditions of GNU GPL v3.
# LDpop
### version 1.0.0
LDpop is a program for computing 2-locus likelihoods under the coalescent with recombination. Unlike previous methods, LDpop can correctly account for variable (but piecewise-constant) population size history.

LDpop produces lookup tables that can be used by other programs to estimate recombination maps. Other programs that can use the output of LDpop include:
* [LDhat](https://github.com/auton1/LDhat)
* [LDhot](https://github.com/auton1/LDhot)
* [LDhelmet](https://sourceforge.net/projects/ldhelmet/)

LDpop also provides utilities for efficient posterior sampling of 2-locus ARGs.

## Installation and Dependencies

Prerequisites:
* Python 2.7, 3.5, or 3.6
* Optional: Java 8
  * Not required for computing lookup tables.
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

[Jack Kamm](mailto:jkamm@stat.berkeley.edu), Jeffrey Spence, Jeffrey Chan, Yun S. Song

## Reference

Kamm, J.A.\*, Spence, J.P.\*, Chan, J., and Song, Y.S. Two-Locus Likelihoods under Variable Population Size and Fine-Scale Recombination Rate Estimation. *Genetics*, 2016, in press. http://www.genetics.org/content/early/2016/05/09/genetics.115.184820

# ldpop
Two locus likelihoods and ARGs under changing population size

## Installation and Dependencies

Prerequisites:
* Scientific distribution of Python 2.7, e.g. [Anaconda](http://continuum.io/downloads), [Enthought Canopy](https://www.enthought.com/products/canopy/)
  * Alternatively, custom installation of pip, the SciPy stack
* For importance sampling: java

To install, in the top-level directory of ldpop (where "setup.py" lives), type
```
pip install .
```

To uninstall, do
```
pip uninstall ldpop
```

## Getting started

```
bin/ldtable.py --help
```

## Authors

[John Kamm](mailto:jkamm@stat.berkeley.edu), Jeffrey Spence, Jeffrey Chan, Yun S. Song

## License

ldpop is not yet publicly released; please do not share with others.

When ldpop is publicly released, it will be free software under conditions of GNU GPL v3.

# TODO

In roughly the order of importance
* Create bin/ImportanceSampler.jar
  * Merge Jeff Chan's changes to the main branch of svn repository, remove unneeded files, and create a single jar
  * Should take output from [bin/ldproposal.py](bin/ldproposal.py), and do importance sampling for a single config or all configs
  * Let the config be chosen through the command line, not through the moranGrid file (so that the same moranGrid could be reused for different configs)
* Finish [example/example_importance_sampling.sh](example/example_importance_sampling.sh)
* Clean up [bin/ldproposal.py](bin/ldproposal.py)
  * Write help/documentation
  * Make command line syntax the same as in [bin/ldtable.py](bin/ldtable.py)
  * Have config chosen through bin/ImportanceSampler.jar, instead of [bin/ldproposal.py](bin/ldproposal.py)
  * Use a single rho, instead of a grid of rhos
* Test ldpop installation for clean install
  * Can use virtualenv, or a virtual machine  
* Make ISProposal a class (a la LookupTable)
  * currently, it is a function returning a str
  * define `ISProposal.__str__` so that `print ISProposal(...)` still works
* Write docstrings for everything imported into [ldpop/__init__.py](ldpop/__init__.py)
  * LookupTable docstring already written
* Write better tests in test/
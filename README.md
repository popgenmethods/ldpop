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
Also, see the [examples](example/)

## Authors

[John Kamm](mailto:jkamm@stat.berkeley.edu), Jeffrey Spence, Jeffrey Chan, Yun S. Song

## License

ldpop is not yet publicly released; please do not share with others.

When ldpop is publicly released, it will be free software under conditions of GNU GPL v3.

# TODO

Tasks, in roughly the order of importance (edited by Jeff S. 3/9/16):
* Test with ldhat, ldhelmet  
* Write docstrings for everything imported into [ldpop/\_\_init\_\_.py](ldpop/__init__.py)
  * LookupTable docstring already written
  * ISProposal docstring already written
* Write better tests in test/

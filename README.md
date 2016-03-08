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

For Jeff & Jeff: `bin/ldtable.py` corresponds to `build_ldhat_table.py` from before. `bin/ldproposal.py` corresponds to `build_timed_grid_table.py` from before. The major remaining task is to add a callable to do the importance sampling -- I am thinking a jar file like `bin/ImportanceSampler.jar`

Tasks, in roughly the order of importance (edited by Jeff S. 3/8/16):
  * Merge Jeff Spence's changes to the main branch of svn repository
* Clean up [bin/ldproposal.py](bin/ldproposal.py)
* Finish [example/example_importance_sampling.sh](example/example_importance_sampling.sh)  
* Test with ldhat, ldhelmet  
* Test ldpop installation for clean install
  * Can use virtualenv, or a virtual machine  
* Make ISProposal a class (a la LookupTable)
  * currently, it is a function returning a str
  * define `ISProposal.__str__` so that `print ISProposal(...)` still works
* Write docstrings for everything imported into [ldpop/\_\_init\_\_.py](ldpop/__init__.py)
  * LookupTable docstring already written
* Write better tests in test/

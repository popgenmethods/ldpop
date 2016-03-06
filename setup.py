#!/usr/bin/env python

from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy

extensions = [Extension("moran_model",
                        sources=["moran_model.pyx"],
                        include_dirs=[numpy.get_include()])]

setup(name='moran2locus',
      version='1.0',
      description='Computing the two-locus sampling probability for variable population size',
      author='Jack Kamm, Jeffrey Spence, Jeffrey Chan, Yun S. Song',
      author_email='jkamm@stat.berkeley.edu, spence.jeffrey@berkeley.edu, chanjed@cs.berkeley.edu, yss@eecs.berkeley.edu',
      packages=['momi'],      
      install_requires=['numpy','scipy'],
      keywords=['population genetics','statistics','moran model','coalescent'],
      ext_modules=cythonize(extensions),      
      )
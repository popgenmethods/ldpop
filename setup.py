#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
import numpy


setup(name='ldpop',
      version='0.0.1',
      description='Computing the two-locus sampling probability for variable population size',
      author='Jack Kamm, Jeffrey Spence, Jeffrey Chan, Yun S. Song',
      author_email='jkamm@stat.berkeley.edu, spence.jeffrey@berkeley.edu, chanjed@cs.berkeley.edu, yss@eecs.berkeley.edu',
      packages=['ldpop'],      
      install_requires=['numpy','scipy','pandas'],
      keywords=['population genetics','statistics','moran model','coalescent'],
      )

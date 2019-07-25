'''
Created on Jan 20, 2015

@author: jkamm
'''
from __future__ import division
from __future__ import absolute_import
from builtins import str
from builtins import zip
from builtins import range

from .compute_stationary import stationary

import numpy
import logging
import time
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import norm


class NumericalError(Exception):
    def __init__(self, message):
        super(NumericalError, self).__init__(message)


# call this before dividing out pi_c
def assert_valid_likelihoods(likelihoods, pi_c, moranRates):
    EPSILON = 1e-7
    if numpy.abs(numpy.sum(likelihoods) - 1.) >= EPSILON:
        raise NumericalError('%g' % (numpy.sum(likelihoods) - 1.0))
    for currC in range(len(pi_c)):
        if pi_c[currC] != 0:
            val = numpy.abs(numpy.log(numpy.sum(
                likelihoods[moranRates.numC == currC]) / pi_c[currC]))
            if val >= EPSILON:
                raise NumericalError(str(
                    [numpy.log(sum(likelihoods[moranRates.numC == currC])
                               / pi_c[currC])
                     for currC in range(len(pi_c)) if pi_c[currC] != 0]))


def folded_likelihoods(moranRates, rho, theta, popSizes, timeLens,
                       gridPointsPerEpoch=0, lastEpochInit=None):
    assert len(popSizes) == len(timeLens) + 1
    timeStart = time.time()

    pi_c = moranRates.get_pi_c(popSize=popSizes[-1], theta=theta, rho=rho)
    renormalize = pi_c[moranRates.numC]
    not_zero = renormalize != 0.

    rates = moranRates.getRates(popSize=popSizes[-1], rho=rho, theta=theta)
    init_stationary = stationary(Q=rates, init=lastEpochInit)
    likelihoods = numpy.copy(init_stationary)

    assert_valid_likelihoods(likelihoods, pi_c, moranRates)
    likelihoods[not_zero] /= renormalize[not_zero]

    ret = {}
    currTime = sum(timeLens)

    for t, popSize in reversed(list(zip(timeLens, popSizes))):
        rates = moranRates.getRates(popSize=popSize, theta=theta, rho=rho)

        pi_c = moranRates.get_pi_c(popSize=popSize, theta=theta, rho=rho)
        renormalize = pi_c[moranRates.numC]
        likelihoods *= renormalize
        if gridPointsPerEpoch == 0:
            start = time.time()
            likelihoods = expm_multiply(rates.transpose() * t, likelihoods)
            end = time.time()
            logging.info('Computed action in %f seconds for rho=%f,N=%f,t=%f'
                         % (end-start, rho, popSize, t))

            not_zero = renormalize != 0.

            assert_valid_likelihoods(likelihoods, pi_c, moranRates)
            likelihoods[not_zero] /= renormalize[not_zero]

            assert norm(likelihoods[numpy.logical_not(not_zero)], 1) < 1e-300
            currTime -= t
        else:
            start = time.time()
            likelihoods = expm_multiply(rates.transpose(),
                                        likelihoods,
                                        endpoint=True,
                                        stop=t,
                                        start=0.0,
                                        num=gridPointsPerEpoch+1)
            end = time.time()
            logging.info('Computed action in %f seconds for rho=%f,N=%f,t=%f'
                         % (end-start, rho, popSize, t))

            for i in range(gridPointsPerEpoch+1):
                currLik = likelihoods[i]
                not_zero = renormalize != 0.

                assert_valid_likelihoods(currLik, pi_c, moranRates)
                currLik[not_zero] /= renormalize[not_zero]

                assert norm(currLik[numpy.logical_not(not_zero)], 1) < 1e-300

                if currTime in ret:
                    assert ret[currTime] == currLik
                else:
                    ret[currTime] = currLik
                currTime -= t / float(gridPointsPerEpoch+1)

            likelihoods = currLik

    assert abs(currTime) < 1e-15, str(currTime)

    ret[0.0] = likelihoods

    timeEnd = time.time()
    logging.info('Finished likelihoods in %f seconds for rho=%f'
                 % (timeEnd - timeStart, rho))

    if gridPointsPerEpoch == 0:
        ret = ret[0.0]
    return (ret, init_stationary)

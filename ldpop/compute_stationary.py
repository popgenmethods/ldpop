from __future__ import division
from builtins import range

import numpy
import logging
import math
import time
from numpy import array
import scipy.sparse as sp
from scipy.linalg import norm


def assertValidProbs(x):
    assert numpy.all(x >= 0) and abs(1.0 - numpy.sum(x)) < 1e-10, numpy.sum(x)


# Q must be tridiagonal and have no absorbing states
def stationary1d_tridiagonal(Q):
    assert Q.shape[0] == Q.shape[1]
    n = Q.shape[0] - 1
    # check Q is tridiagonal
    non_tri_nonzero_entries = Q != 0.0
    non_tri_nonzero_entries[numpy.arange(n+1), numpy.arange(n+1)] = False
    non_tri_nonzero_entries[numpy.arange(n), numpy.arange(1, n+1)] = False
    non_tri_nonzero_entries[numpy.arange(1, n+1), numpy.arange(n)] = False

    if numpy.any(non_tri_nonzero_entries.dot(numpy.ones(n+1)) > 0):
        raise Exception('Q must be strictly tridiagonal')

    ret = [0.0]
    for i in range(n):
        # pi[i] * Q[i,i+1] = pi[i+1] * Q[i+1,i]
        ret.append(ret[i] + math.log(Q[i, i+1]) - math.log(Q[i+1, i]))
    ret = numpy.array(ret)
    ret = ret - numpy.max(ret)
    ret = numpy.exp(ret)
    ret = ret / numpy.sum(ret)
    return ret


def stationary(Q, init=None, norm_order=1, epsilon=1e-8):
    start = time.time()

    size = Q.shape[0]
    assert size == Q.shape[1]
    # Compute a suitable stochastic matrix by means of uniformization
    ell = Q.min() * 1.001
    P = sp.eye(size, size) - Q/ell
    # compute Pi
    P = P.tocsr()
#     pi = zeros(size);  pi1 = zeros(size)
#     pi[0] = 1;
#     n = norm(pi - pi1,1)
    if init is not None:
        assertValidProbs(init)
        pi = init
    else:
        pi = array([1 / float(size) for i in range(size)])
    n = float('inf')

    i = 0

    while n > epsilon:
        pi1 = pi*P
        pi = pi1*P   # avoid copying pi1 to pi

        not_zero = numpy.logical_and(pi != 0., pi1 != 0.)
        if numpy.all((pi == 0.) == (pi1 == 0.)):
            n = norm(numpy.log(pi[not_zero] / pi1[not_zero]), norm_order)
            i += 1
        else:
            n = float('inf')
    end = time.time()
    logging.info('Computed stationary with L-%f stopping, epsilon=%g,'
                 'in %d iterations, %f seconds '
                 % (norm_order, epsilon, i, end - start))
#     logging.info( 'Done with Stationary' )
    return pi

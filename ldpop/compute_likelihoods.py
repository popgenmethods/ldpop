'''
Created on Jan 20, 2015

@author: jkamm
'''

from compute_stationary import stationary

import numpy, scipy, math, logging, time
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import norm

# def ordered_log_likelihoods(states, rho, theta, popSizes, timeLens, gridPointsPerEpoch=0,  lastEpochInit=None):
#     return states.ordered_log_likelihoods(folded_likelihoods(states, rho, theta, popSizes, timeLens, gridPointsPerEpoch, lastEpochInit))

class NumericalError(StandardError):
    def __init__(self, message):
        super(NumericalError, self).__init__(message)
    

''' call this before dividing out pi_c'''
def assert_valid_likelihoods(likelihoods, pi_c, states):
    EPSILON = 1e-7
    if abs(sum(likelihoods) - 1.) >= EPSILON:
        raise NumericalError( "%g" % (sum(likelihoods) - 1.0))
    if not all([abs( math.log(sum(likelihoods[states.numC == currC]) / pi_c[currC]) ) < EPSILON for currC in range(len(pi_c)) if pi_c[currC] != 0]):
        raise NumericalError(str([math.log(sum(likelihoods[states.numC == currC]) / pi_c[currC]) for currC in range(len(pi_c)) if pi_c[currC] != 0]))


def folded_likelihoods(states, rho, theta, popSizes, timeLens, gridPointsPerEpoch=0, lastEpochInit=None):
    assert len(popSizes) == len(timeLens) + 1
    timeStart = time.time()
        
    pi_c = states.get_pi_c(popSize=popSizes[-1], theta=theta, rho=rho)
    renormalize = pi_c[states.numC]
    not_zero = renormalize != 0.

    rates = states.getRates(popSize=popSizes[-1],rho=rho,theta=theta)
    likelihoods = stationary(Q=rates, init=lastEpochInit)
    
    assert_valid_likelihoods(likelihoods, pi_c, states)
    likelihoods[not_zero] /= renormalize[not_zero]

    ret = {}
    currTime = sum(timeLens)

    for t,popSize in reversed(zip(timeLens, popSizes)):
        rates = states.getRates(popSize=popSize,theta=theta,rho=rho)
                
        pi_c = states.get_pi_c(popSize=popSize, theta=theta, rho=rho)
        renormalize = pi_c[states.numC]
        likelihoods *= renormalize
        if gridPointsPerEpoch == 0:
            start = time.time()
            likelihoods = expm_multiply(rates.transpose() * t, likelihoods)
            end = time.time()
            logging.info("Computed action in %f seconds for rho=%f,N=%f,t=%f" % (end-start, rho, popSize, t))
            
            not_zero = renormalize != 0.
            
            assert_valid_likelihoods(likelihoods, pi_c, states)
            likelihoods[not_zero] /= renormalize[not_zero]
            
            assert norm(likelihoods[numpy.logical_not(not_zero)], 1) < 1e-300
            currTime -= t
        else:
            start = time.time()            
            likelihoods = expm_multiply(rates.transpose(), likelihoods, endpoint=True, stop=t, start=0.0, num=gridPointsPerEpoch+1)
            end = time.time()
            logging.info("Computed action in %f seconds for rho=%f,N=%f,t=%f" % (end-start, rho, popSize, t))
                
            for i in range(gridPointsPerEpoch+1):
                currLik = likelihoods[i]
                not_zero = renormalize != 0.
                
                assert_valid_likelihoods(currLik, pi_c, states)
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
    logging.info("Finished likelihoods in %f seconds for rho=%f" % (timeEnd - timeStart, rho))

    if gridPointsPerEpoch == 0:
        ret = ret[0.0]
    
    return ret

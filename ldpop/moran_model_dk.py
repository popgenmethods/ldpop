'''
Created on Jan 23, 2015

@author: jkamm
'''
from compute_stationary import stationary
from moran_model import AbstractMoranStates, makeFrozen, c_haps, makeAllConfigs

import logging, time, numpy, scipy
from scipy import sparse
from collections import Counter
from functools import partial
                
class MoranStatesDK(AbstractMoranStates):
    '''
    maintains a representation of the states of the 2 locus Moran model
    '''

    def __init__(self, n):
        '''
        Constructor
        '''
        start = time.time()
        super(MoranStatesDK, self).__init__(n)
        self.exact = False        
        
        # make all haplotypes
        self.hapList = []
        for allele1 in xrange(2):
            for allele2 in xrange(2):
                self.hapList.append((allele1,allele2))


        self.build_all_configs(n, exact=False)
        
        end = time.time()
        logging.info("Constructed approx states in %f seconds" % (end - start))

        #self.folded_configs = self.all_configs
        #self.allIdx_to_foldedIdx = range(len(self.folded_configs))
        self.build_symmetries()        
        
        start = time.time()        
        self.unscaled_recom_rates = self.build_rates(partial(config_recom_rates, self.n))
        logging.info("Constructed recombination rate matrix in %f seconds" % (time.time() - start))

        start = time.time()        
        self.unscaled_mut_rates = self.build_rates(config_mut_rates)
        logging.info("Constructed mut rate matrix in %f seconds" % (time.time() - start))

        start = time.time()
        self.unscaled_coal_rates = self.build_rates(config_copy_rates)
        logging.info("Constructed copying rate matrix in %f seconds" % (time.time() - start))
        
    # def getRates(self, popSize, theta, rho):
    #     return MoranRatesDK(states=self, rho=rho, popSize=popSize, theta=theta).rates
    
def config_copy_rates(state):
    ret = Counter()
    for hap,numHap in state:
        # hap copying event
        if numHap > 0:
            for otherHap in c_haps:
                otherState = dict(state)                
                numOtherHap = otherState[otherHap]

                if otherHap == hap or numOtherHap == 0:
                    continue

                rate = numHap * numOtherHap / float(2)

                otherState[otherHap] -= 1
                otherState[hap] += 1
                # make it hashable
                otherState = makeFrozen(otherState)

                ret[state, otherState] += rate
    return ret
    
def config_mut_rates(state):
    ret = Counter()
    for hap,numHap in state:
        # mutation
        if numHap > 0:
            rate = numHap
            for loc in xrange(2):
                otherHap = [hap[0], hap[1]]
                otherAllele = 1 - hap[loc]
                otherHap[loc] = otherAllele
                otherHap = tuple(otherHap)

                otherState = dict(state)
                otherState[hap] -=1
                otherState[otherHap] += 1
                # make it hashable
                otherState = makeFrozen(otherState)

                ret[state, otherState] += rate
    return ret
        
# def config_recom_rates(n, state):
#     ret = Counter()
#     # recombination
#     for leftHap in c_haps:
#         for rightHap in c_haps:
#             for replacedHap in c_haps:
#                 config = dict(state)
#                 numWays = config[replacedHap]
#                 config[replacedHap] -= 1
#                 numWays *= config[leftHap]
#                 config[leftHap] -= 1
#                 numWays *= config[rightHap]
#                 config[rightHap] -= 1

#                 if numWays == 0:
#                     continue

#                 config = dict(state)
#                 config[replacedHap] -= 1
#                 config[(leftHap[0], rightHap[1])] += 1
#                 otherState = makeFrozen(config)

#                 assert n >= 3
#                 rate = numWays / float(n-1) / float(n-2)
#                 ret[state, otherState] += rate
#     return ret

def config_recom_rates(n, state):
    ret = Counter()
    # recombination
    for mom in c_haps:
        for dad in c_haps:
            config = dict(state)
            numWays = config[mom]
            config[mom] -= 1
            numWays *= config[dad]
            config[dad] -= 1

            if numWays == 0:
                continue

            config[(mom[0], dad[1])] += 1
            config[(dad[0], mom[1])] += 1
            otherState = makeFrozen(config)

            rate = numWays / float(n-1) / float(2)
            ret[state, otherState] += rate
    return ret

'''
Created on Jan 23, 2015

@author: jkamm
'''
from __future__ import absolute_import
from builtins import range
from .compute_stationary import stationary
from .moran_augmented import AbstractMoranStates, c_haps, makeAllConfigs, build_mut_rates, build_copy_rates, get_rates, subtract_rowsum_on_diag

import logging, time, numpy, scipy
from scipy import sparse
from collections import Counter
from functools import partial
                
class MoranStatesFinite(AbstractMoranStates):
    '''
    maintains a representation of the states of the 2 locus Moran model
    '''

    def __init__(self, n):
        '''
        Constructor
        '''
        start = time.time()
        super(MoranStatesFinite, self).__init__(n)
        self.exact = False        
        
        # make all haplotypes
        self.hapList = []
        for allele1 in range(2):
            for allele2 in range(2):
                self.hapList.append((allele1,allele2))


        self.build_all_configs(n, exact=False)
        
        end = time.time()
        logging.info("Constructed approx states in %f seconds" % (end - start))

        #self.folded_configs = self.all_configs
        #self.allIdx_to_foldedIdx = range(len(self.folded_configs))
        self.build_symmetries()        
        
        start = time.time()        
        #self.unscaled_recom_rates = self.build_rates(partial(config_recom_rates, self.n))
        self.unscaled_recom_rates = build_recom_rates(self)
        logging.info("Constructed recombination rate matrix in %f seconds" % (time.time() - start))

        start = time.time()        
        self.unscaled_mut_rates = build_mut_rates(self)
        logging.info("Constructed mut rate matrix in %f seconds" % (time.time() - start))

        start = time.time()
        #self.unscaled_coal_rates = self.build_rates(config_copy_rates)
        self.unscaled_coal_rates = build_copy_rates(self)
        logging.info("Constructed copying rate matrix in %f seconds" % (time.time() - start))
        

def build_recom_rates(states):
    ret = sparse.csr_matrix(tuple([states.folded_config_array.shape[0]]*2))    
    confs = states.folded_config_array
    for mom in c_haps:
        for dad in c_haps:
            otherConfs = numpy.array(confs)
            
            rates = numpy.array(otherConfs[:,mom[0],mom[1]], dtype=float)
            otherConfs[:,mom[0],mom[1]] -= 1
            
            rates *= otherConfs[:,dad[0],dad[1]]
            otherConfs[:,dad[0],dad[1]] -= 1

            rates *= 1. / float(states.n-1) / float(2)
            
            otherConfs[:,mom[0], dad[1]] += 1
            otherConfs[:,dad[0], mom[1]] += 1

            ret = ret + get_rates(states, otherConfs, rates)

    return subtract_rowsum_on_diag(ret)

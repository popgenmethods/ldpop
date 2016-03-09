'''
Created on Jan 19, 2015

@author: jkamm
'''
from compute_stationary import stationary, stationary1d_tridiagonal, assertValidProbs

import logging, time, numpy, math, scipy, scipy.misc, scipy.special
from scipy import sparse
from collections import Counter
import multiprocessing as mp


## Have a separate class for the Rates
## so we don't have to pickle all of MoranStatesAugmented
## when doing multiprocessing
## (saves memory and communication time with worker processes)
class MoranRates(object):
    def __init__(self, states):
        self.exact = states.exact
        self.numC = states.numC
        self.n = states.n
        self.unscaled_recom_rates = states.unscaled_recom_rates
        self.unscaled_mut_rates = states.unscaled_mut_rates
        self.unscaled_coal_rates = states.unscaled_coal_rates

    def get_pi_c(self, popSize, theta, rho):
        if not self.exact:
            return numpy.array([0.0] * self.n + [1.0])
        n=self.n
        coalRate = 1. / popSize
        recomRate = float(rho) / 2.

        if rho == 0.0:
            return numpy.array([0.0] * self.n + [1.0])
        else:        
            numCoupledLinsRates = sparse.dok_matrix((n+1, n+1))
            for i in xrange(n+1):
                if i < n:
                    numCoupledLinsRates[i,i+1] = ((n-i)**2) * coalRate
                    numCoupledLinsRates[i,i] -= numCoupledLinsRates[i,i+1]
                if i > 0:
                    numCoupledLinsRates[i,i-1] = recomRate * i
                    numCoupledLinsRates[i,i] -= numCoupledLinsRates[i,i-1]
            return stationary1d_tridiagonal(numCoupledLinsRates)            

    def getRates(self, popSize, theta, rho):
        #return MoranRates(states=self, rho=rho, popSize=popSize, theta=theta).rates
        start = time.time()
        recomRate = float(rho) / 2.
        mutRate = float(theta) / 2.
        coalRate = 1. / float(popSize)
        #self.rates = nonRecomRates.rates + recomRate * get_recom_rates(self.n)
        ret = recomRate * self.unscaled_recom_rates + mutRate * self.unscaled_mut_rates + coalRate * self.unscaled_coal_rates
        end = time.time()
        logging.info("%f seconds to construct rates for rho=%f,theta=%f,N=%f" % (end-start, rho,theta,popSize))
        return ret

# make all haplotypes
a_haps = []; b_haps = []; c_haps = []
for allele1 in xrange(2):
    a_haps.append((allele1,-1))
    b_haps.append((-1, allele1))
    for allele2 in xrange(2):
        c_haps.append((allele1,allele2))

all_haps = a_haps + b_haps + c_haps

def makeAllConfigs(hapList, n):
    # make all configs
    # represent a config as a dict
    tmpConfigList = [{}]
    for hap in hapList:
        newConfigList = []
        for config in tmpConfigList:
            numHaps = sum([v for k,v in config.items()])
            assert numHaps <= n
            for i in xrange(n - numHaps + 1):
                newConfig = dict(config)
                newConfig[hap] = i
                newConfigList.append(newConfig)
        tmpConfigList = newConfigList
        
    configList = [[] for i in xrange(n+1)]
    for config in tmpConfigList:
        numHaps = sum([v for k,v in config.items()])
        configList[numHaps].append(config)
    return configList

        
def makeFrozen(dict):
    #return frozenset(dict.items())
    return tuple(sorted(dict.items()))

def one_locus_probs(popSize, theta, n):
    coalRate = 1. / popSize
    mutRate = float(theta) / 2.
    
    numOnesRates = sparse.dok_matrix((n+1,n+1))
    for i in xrange(n+1):
        if i < n:
            numOnesRates[i,i+1] = (n-i) * mutRate + i * (n-i) / 2.0 * coalRate
            numOnesRates[i,i] -= numOnesRates[i,i+1]
        if i > 0:
            numOnesRates[i,i-1] = i * mutRate + i * (n-i) / 2.0 * coalRate
            numOnesRates[i,i] -= numOnesRates[i,i-1]
    
    return stationary1d_tridiagonal(numOnesRates)
           
class AbstractMoranStates(object):
    def __init__(self, n):
        self.n = n
        self._stationary = {}

    def build_all_configs(self, n, exact):
        # now make a config list. make all the configs frozensets so they are hashable
        aConfigs = makeAllConfigs(a_haps, n)
        bConfigs = makeAllConfigs(b_haps, n)
        cConfigs = makeAllConfigs(c_haps, n)
        
        self.all_configs = []
        self.all_configs_idx = {}
        self.full_configs_all = []
        #self.numC = []

        if exact:
            cList = range(n+1)
        else:
            cList = [n]
        
        for numC in cList:
            for aConf in aConfigs[n - numC]:
                for bConf in bConfigs[n - numC]:
                    for cConf in cConfigs[numC]:
                        conf = {}
                        conf.update(aConf); conf.update(bConf); conf.update(cConf)
                        
                        assert len(conf) == len(all_haps)
                        conf = makeFrozen(conf)
                        
                        self.all_configs_idx[conf] = len(self.all_configs)
                        self.all_configs.append(conf)

                        if numC == n:
                            self.full_configs_all.append(conf)      
                    
        assert len(self.all_configs) == len(self.all_configs_idx)
        assert len(self.all_configs_idx) > 0

        self.build_config_array()

    def build_config_array(self):
        """
        Create self.config_array, defined by:
        self.config_array[i, a, b] = the count of haplotype (a,b) in the i-th config
        """
        idxs, vals = zip(*[zip(*config) for config in self.all_configs])
        idxs = numpy.array(idxs) # axis0=config, axis1=hap, axis2=locus
        vals = numpy.array(vals) # axis0=config, axis1=hap

        # now label each entry by the index of its config
        
        idxs = idxs.T # move the config axis at the end, the locus axis to the front
        # broadcast up the indices [0,1,2,3,...,max_config] to have same dimension as idxs
        # idxs becomes the pair (config_idxs, idxs)
        idxs = numpy.broadcast_arrays(numpy.arange(len(self.all_configs)),
                                      idxs)
        # now stack the config_idxs onto the idxs, along the locus axis
        idxs = numpy.vstack(idxs)[1:,:,:] # take [1:,...] because we copied the config_idx twice (for each locus)
        # now move config index to the front, locus index to the back
        idxs = idxs.T
        assert idxs.shape == (len(self.all_configs), 8, 3) # axis0=config, axis1=hap, axis2=configIdx/locus

        # flatten the counts
        vals = numpy.reshape(vals,-1,order='C')
        # flatten the configs (but keep the index/locus axis
        idxs = numpy.reshape(idxs,(-1,3),order='C')

        # now create the config array, and assign all the counts to it
        self.config_array = numpy.zeros((len(self.all_configs), 3, 3), dtype=int)
        self.config_array[idxs[:,0],idxs[:,1],idxs[:,2]] = vals

        # create dictionary mapping their hash values back to their index
        hash_vals = self.hash_config_array(self.config_array)
        assert len(set(hash_vals)) == len(hash_vals) # should be all unique
        self.hash_to_allIdx = {k:v for v,k in enumerate(hash_vals)}

    def hash_config_array(self, conf_arr):
        base = self.n+1
        hash_vals = conf_arr[:,0,0] + base * conf_arr[:,0,1] + (base**2) * (conf_arr[:,1,0])
        if self.exact:
            hash_vals += (base**3)*(conf_arr[:,1,1]) + (base**4)*(conf_arr[:,0,-1]) + (base**5)*(conf_arr[:,-1,0])
        return hash_vals
        
    def numOnes(self, loc):
        return self.folded_config_array.sum(axis=1+(1-loc))[:,1]
        
    def hapCount(self, hap):
        return numpy.array(self.folded_config_array[:,hap[0],hap[1]])
    
    def getUnlinkedStationary(self, popSize, theta):
        one_loc_probs = one_locus_probs(popSize=popSize, theta=theta, n=self.n)
        assertValidProbs(one_loc_probs)
        
        n = self.n
        leftOnes, rightOnes, bothOnes = self.numOnes(0), self.numOnes(1), self.hapCount((1,1))
#             joint = one_loc_probs[leftOnes] * one_loc_probs[rightOnes]
        joint = one_loc_probs[leftOnes] * one_loc_probs[rightOnes]
        if self.exact:
            joint[self.numC > 0] = 0
        else:
            joint = joint * scipy.misc.comb(rightOnes, bothOnes) * scipy.misc.comb(n-rightOnes, leftOnes-bothOnes)  / scipy.misc.comb(n, leftOnes)
           
        joint = joint * self.n_unfolded_versions

        assertValidProbs(joint)  
        return joint
   
    def build_symmetries(self):
        start = time.time()

        # the index of the folded version in all_configs
        folded_list = get_folded_config_idxs(self)
        
        # foldedIdx = the index in folded_configs, allIdx = the index in all_configs
        foldedIdx_to_allIdx = numpy.array(list(set(folded_list)))
        
        allIdx_to_foldedIdx = {v:k for k,v in enumerate(foldedIdx_to_allIdx)}       
        allIdx_to_foldedIdx = [allIdx_to_foldedIdx[x] for x in folded_list]

        self.hash_to_foldedIdx = {k: allIdx_to_foldedIdx[v] for k,v in self.hash_to_allIdx.iteritems()}
        self.folded_config_array = self.config_array[foldedIdx_to_allIdx,:,:]        
              
        self.numC = self.folded_config_array[:,0,0] + self.folded_config_array[:,0,1] + self.folded_config_array[:,1,0] + self.folded_config_array[:,1,1]
        
        symm_mat = sparse.dok_matrix((len(allIdx_to_foldedIdx), self.folded_config_array.shape[0]))
        symm_mat.update(dict(zip(enumerate(allIdx_to_foldedIdx), [1]*len(folded_list))))
        symm_mat = symm_mat.tocsc()
        
        antisymm_mat = symm_mat.transpose().tocsr(copy=True)
        # normalize rows
        self.n_unfolded_versions = numpy.array(antisymm_mat.sum(axis=1))[:,0]
        row_indices, col_indices = antisymm_mat.nonzero()
        antisymm_mat.data /= self.n_unfolded_versions[row_indices]
       
        self.symmetries = symm_mat.tocsr()
        self.antisymmetries = antisymm_mat.tocsr()
        
        logging.info("%f seconds to build symmetry matrices" % (time.time() - start))

    def unfold_likelihoods(self, liks):
        try:
            return {t: self.unfold_likelihoods(v) for t,v in liks.items()}
        except AttributeError:
            liks = liks * self.antisymmetries
            return {config : liks[self.all_configs_idx[config]] for config in self.full_configs_all}

    def ordered_log_likelihoods(self, liks):
        unordered_likelihoods = self.unfold_likelihoods(liks)

        try:
            return {time : single_ordered_log_likelihood(self.n, unordered) for time,unordered in unordered_likelihoods.items()}
        except AttributeError:
            return single_ordered_log_likelihood(self.n, unordered_likelihoods)    

        
def log_n_factorial(n):
    return scipy.special.gammaln(n+1)

def single_ordered_log_likelihood(n, unordered_likelihoods):
    ordered_ll = {}
    for config, likelihood in unordered_likelihoods.items():
        comb = log_n_factorial(n)
        for hap,count in config:
            comb -= log_n_factorial(count)
            
        ordered_ll[config] = math.log(likelihood) - comb
    return ordered_ll
   
class MoranStatesAugmented(AbstractMoranStates):
    '''
    maintains a representation of the states(possible configs) of the 2 locus Moran model
    '''
    def __init__(self, n):
        '''
        Constructor
        '''
        #os.system('taskset -p %s >&2' %os.getpid())
       
        start = time.time()
        super(MoranStatesAugmented, self).__init__(n)
        self.exact = True
                       
        self.build_all_configs(n, exact=True)
        
        end = time.time()
        logging.info("Constructed exact states in %f seconds" % (end - start))

        self.build_symmetries()
        
        start = time.time()        
        self.unscaled_recom_rates = build_recom_rates(self)
        logging.info("Constructed recombination rate matrix in %f seconds" % (time.time() - start))

        start = time.time()        
        #self.unscaled_mut_rates = self.build_rates(config_mut_rates)
        self.unscaled_mut_rates = build_mut_rates(self)
        logging.info("Constructed mut rate matrix in %f seconds" % (time.time() - start))

        start = time.time()
        #self.unscaled_coal_rates = self.build_rates(config_copy_rates) + self.build_rates(config_coal_rates)
        #self.unscaled_coal_rates = self.build_rates(config_copy_rates) + build_cross_coal_rates(self)
        self.unscaled_coal_rates = build_copy_rates(self) + build_cross_coal_rates(self)                
        logging.info("Constructed coalescent/copying rate matrix in %f seconds" % (time.time() - start))
               

def get_folded_config_idxs(states):
    arr = states.config_array
    # move the missing allele in between alleles 0,1
    arr = arr[:,(0,-1,1),:][:,:,(0,-1,1)]

    # relabel alleles 0,1 (4 ways to do this)
    symm_arrs = [arr, arr[:,::-1,:], arr[:,:,::-1], arr[:,::-1,::-1]]
    # swap the 2 loci
    symm_arrs += [numpy.transpose(a, axes=(0,2,1)) for a in symm_arrs]

    # swap back allele 1 with missing allele
    symm_arrs = [a[:,(0,-1,1),:][:,:,(0,-1,1)] for a in symm_arrs]
    
    # get hash val for each (folded) config
    hash_vals = numpy.vstack(map(states.hash_config_array, symm_arrs))
    # get the smallest hash val among all the folds
    hash_vals = numpy.amin(hash_vals, axis=0)
    assert len(hash_vals) == arr.shape[0]

    # return the corresponding indices
    return [states.hash_to_allIdx[h] for h in hash_vals]
        

def build_recom_rates(states):
    assert states.exact
    
    ret = sparse.csr_matrix(tuple([states.folded_config_array.shape[0]]*2))    
    confs = states.folded_config_array

    for hap in c_haps:
        rates = confs[:,hap[0],hap[1]]
        
        otherConfs = numpy.array(confs)
        otherConfs[:,hap[0],hap[1]] -= 1
        otherConfs[:,hap[0],-1] += 1
        otherConfs[:,-1,hap[1]] += 1

        ret = ret + get_rates(states, otherConfs, rates)

    return subtract_rowsum_on_diag(ret)

def build_mut_rates(states):
    ret = sparse.csr_matrix(tuple([states.folded_config_array.shape[0]]*2))    
    confs = states.folded_config_array

    if states.exact:
        hapList = all_haps
    else:
        hapList = c_haps
    
    for hap in hapList:
        rates = confs[:,hap[0],hap[1]]
        for loc in xrange(2):
            if hap[loc] == -1:
                continue
            otherHap = [hap[0], hap[1]]
            otherAllele = 1 - hap[loc]
            otherHap[loc] = otherAllele

            otherConfs = numpy.array(confs)
            otherConfs[:,hap[0],hap[1]] -= 1
            otherConfs[:,otherHap[0], otherHap[1]] += 1

            ret = ret + get_rates(states, otherConfs, rates)

    return subtract_rowsum_on_diag(ret)        

def build_copy_rates(states):
    ret = sparse.csr_matrix(tuple([states.folded_config_array.shape[0]]*2))    
    confs = states.folded_config_array

    if states.exact:
        hapList = all_haps
    else:
        hapList = c_haps
    
    for hap in hapList:
        for otherHap in hapList:
            # check if we can copy
            canCopy = True
            for loc in xrange(2):
                if hap[loc] == -1 and otherHap[loc] != -1:
                    canCopy = False
            if not canCopy:
                continue

            copiedHap = [hap[0], hap[1]]
            for loc in xrange(2):
                if otherHap[loc] == -1:
                    copiedHap[loc] = -1
            copiedHap = tuple(copiedHap)

            hapMissing = (hap[0] == -1) + (hap[1] == -1)
            otherMissing = (otherHap[0] == -1) + (otherHap[1] == -1)
            assert otherMissing >= hapMissing

            rates = confs[:,hap[0],hap[1]] * confs[:,otherHap[0],otherHap[1]] / 2.
            if otherMissing > hapMissing:
                rates *= 2

            otherConfs = numpy.array(confs)
            otherConfs[:,otherHap[0],otherHap[1]] -=1
            otherConfs[:,copiedHap[0],copiedHap[1]] += 1
            
            ret = ret + get_rates(states, otherConfs, rates)

    return subtract_rowsum_on_diag(ret)


def subtract_rowsum_on_diag(spmat):
    spmat = spmat.tocsr() - sparse.diags(numpy.array(spmat.sum(axis=1)).T, offsets=[0], format="csr")
    return spmat.tocsr()

def build_cross_coal_rates(states):
    assert states.exact
    ret = sparse.csr_matrix(tuple([states.folded_config_array.shape[0]]*2))

    confs = states.folded_config_array
    for hap in c_haps:
        # two lineages coalesce into hap        
        #leftHap = (hap[0],-1)
        #rightHap = (-1,hap[1])

        otherConfs = numpy.array(confs)
        
        rates = otherConfs[:,hap[0],-1] * otherConfs[:,-1,hap[1]]

        otherConfs[:,hap[0],hap[1]] += 1
        otherConfs[:,hap[0],-1] -= 1
        otherConfs[:,-1,hap[1]] -= 1

        ret = ret + get_rates(states, otherConfs, rates)

    return subtract_rowsum_on_diag(ret)
        
def get_rates(states, otherConfs, rates):
    otherConfs = otherConfs[rates != 0, :,:]
    otherConfs = states.hash_config_array(otherConfs)
    otherConfs = numpy.array([states.hash_to_foldedIdx[x] for x in otherConfs], dtype=int)
    
    confs = numpy.arange(states.folded_config_array.shape[0], dtype=int)
    confs = confs[rates != 0]

    rates = rates[rates!=0]

    ret = sparse.coo_matrix((rates, (confs, otherConfs)), shape=[states.folded_config_array.shape[0]]*2)
    return ret.tocsr()

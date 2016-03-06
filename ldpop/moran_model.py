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
## so we don't have to pickle all of MoranStates
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
        self.numC = []

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
                        self.numC.append(numC)
                        if numC == n:
                            self.full_configs_all.append(conf)      
        
        self.numC = numpy.array(self.numC)
            
        assert len(self.all_configs) == len(self.all_configs_idx)
        assert len(self.all_configs_idx) > 0

    def numOnes(self, loc):
        ret = []
        for config in self.folded_configs:
            curr = 0
            for hap,count in config:
                if hap[loc] == 1:
                    curr += count
            ret.append(curr)
        return numpy.array(ret)       

        
    # def numOnes(self, loc):
    #     ret = []
    #     for config in self.all_configs:
    #         curr = 0
    #         for hap,count in config:
    #             if hap[loc] == 1:
    #                 curr += count
    #         ret.append(curr)
    #     return numpy.array(ret)
    
    def hapCount(self, hap):
        ret = []
        #for config in self.all_configs:
        for config in self.folded_configs:
            config = dict(config)
            ret.append(config[hap])
        return numpy.array(ret)
            
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

    def get_folded_idx(self, state):
        return self.allIdx_to_foldedIdx[self.all_configs_idx[state]]
    
    def build_symmetries(self):
        start = time.time()

        # the index of the folded version in all_configs
        folded_list = map(get_folded_config, self.all_configs)
        folded_list = [self.all_configs_idx[folded_conf] for folded_conf in folded_list]

        # foldedIdx = the index in folded_configs, allIdx = the index in all_configs
        foldedIdx_to_allIdx = list(set(folded_list))
        allIdx_to_foldedIdx = {v:k for k,v in enumerate(foldedIdx_to_allIdx)}
        
        self.allIdx_to_foldedIdx = [allIdx_to_foldedIdx[x] for x in folded_list]
        self.folded_configs = [self.all_configs[foldedIdx_to_allIdx[x]] for x in range(len(foldedIdx_to_allIdx))]
        self.numC = numpy.array([self.numC[foldedIdx_to_allIdx[i]] for i in range(len(self.folded_configs))],
                                dtype=int)
        
        symm_mat = sparse.dok_matrix((len(self.allIdx_to_foldedIdx), len(self.folded_configs)))
        symm_mat.update(dict(zip(enumerate(self.allIdx_to_foldedIdx), [1]*len(folded_list))))
        symm_mat = symm_mat.tocsc()
        
        antisymm_mat = symm_mat.transpose().tocsr(copy=True)
        # normalize rows
        self.n_unfolded_versions = numpy.array(antisymm_mat.sum(axis=1))[:,0]
        row_indices, col_indices = antisymm_mat.nonzero()
        antisymm_mat.data /= self.n_unfolded_versions[row_indices]
       
        self.symmetries = symm_mat.tocsr()
        self.antisymmetries = antisymm_mat.tocsr()
        
        logging.info("%f seconds to build symmetry matrices" % (time.time() - start))

    def build_rates(self, rate_fun):
        rates = {}
        for rate in map(rate_fun, self.folded_configs):
            rates.update(rate)

        folded_rates = Counter()
        for (i,j),v in rates.iteritems():
            folded_rates[self.get_folded_idx(i), self.get_folded_idx(j)] += v

        ret = sparse.dok_matrix((len(self.folded_configs), len(self.folded_configs)))
        ret.update(folded_rates)

        ret = ret.tocsr() - sparse.diags(numpy.array(ret.sum(axis=1)).T, offsets=[0], format="csr")
        return ret.tocsr()

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
   
class MoranStates(AbstractMoranStates):
    '''
    maintains a representation of the states(possible configs) of the 2 locus Moran model
    '''
    def __init__(self, n):
        '''
        Constructor
        '''
        #os.system('taskset -p %s >&2' %os.getpid())
       
        start = time.time()
        super(MoranStates, self).__init__(n)
        self.exact = True
                       
        self.build_all_configs(n, exact=True)
        
        end = time.time()
        logging.info("Constructed exact states in %f seconds" % (end - start))

        self.build_symmetries()
        
        #self.unscaled_recom_rates = get_recom_rates(self)
        start = time.time()        
        self.unscaled_recom_rates = self.build_rates(config_recom_rates)
        logging.info("Constructed recombination rate matrix in %f seconds" % (time.time() - start))

        start = time.time()        
        self.unscaled_mut_rates = self.build_rates(config_mut_rates)
        logging.info("Constructed mut rate matrix in %f seconds" % (time.time() - start))

        start = time.time()
        self.unscaled_coal_rates = self.build_rates(config_copy_rates) + self.build_rates(config_coal_rates)
        logging.info("Constructed coalescent/copying rate matrix in %f seconds" % (time.time() - start))
               

def get_folded_config(state):
    arr = config2array(state)
    symm_arrs = [arr, arr[::-1,:], arr[:,::-1], arr[::-1,::-1]]
    symm_arrs += [a.T for a in symm_arrs]
    symm_dicts = map(array2dict, symm_arrs)
    return makeFrozen(min(symm_dicts))


def config_recom_rates(state):
    ret = Counter()
    for hap,numHap in state:
        if numHap == 0 or hap[0] == -1 or hap[1] == -1:
            continue

        leftHap = (hap[0],-1)
        rightHap = (-1,hap[1])

        # hap recombine
        otherState = dict(state)
        rate = numHap

        otherState[hap] -= 1
        otherState[leftHap] += 1
        otherState[rightHap] += 1                        

        # make it hashable
        otherState = makeFrozen(otherState)

        ret[state, otherState] += rate
        #ret[state, state] -= rate
    return dict(ret)

def config_mut_rates(state):
    fullRates = Counter()
    for hap, numHap in state:
        # hap mutates to something else
        if numHap == 0:
            continue
        rate = numHap
        for loc in xrange(2):
            if hap[loc] == -1:
                continue
            otherHap = [hap[0], hap[1]]
            otherAllele = 1 - hap[loc]
            otherHap[loc] = otherAllele
            otherHap = tuple(otherHap)

            otherState = dict(state)
            otherState[hap] -=1
            otherState[otherHap] += 1
            # make it hashable
            otherState = makeFrozen(otherState)

            fullRates[state, otherState] += rate
            #fullRates[state, state] -= rate
    return dict(fullRates)

def config_copy_rates(state):
    fullRates = Counter()

    for hap,numHap in state:
        # hap copying event
        if numHap == 0:
            continue

        for otherHap,numOtherHap in state:
            if numOtherHap == 0:
                continue
            
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

            rate = numHap * numOtherHap / 2.
            if otherMissing > hapMissing:
                rate *= 2

            otherState = dict(state)
            otherState[otherHap] -= 1
            otherState[copiedHap] += 1
            # make it hashable
            otherState = makeFrozen(otherState)

            fullRates[state, otherState] += rate
            #fullRates[state, state] -= rate
    return dict(fullRates)

def config_coal_rates(state):
    fullRates = Counter()

    for hap in all_haps:
        if hap[0] == -1 or hap[1] == -1:
            continue
        leftHap = (hap[0],-1)
        rightHap = (-1,hap[1])

        # two lineages coalesce into hap
        otherState = dict(state)
        if otherState[leftHap] == 0 or otherState[rightHap] == 0:
            continue
        
        rate = otherState[leftHap] * otherState[rightHap]

        otherState[hap] += 1
        otherState[leftHap] -= 1
        otherState[rightHap] -= 1
        # make it hashable
        otherState = makeFrozen(otherState)

        fullRates[state, otherState] += rate
        #fullRates[state, state] -= rate

    return dict(fullRates)


        
def allele2idx(allele):
    if allele == 0:
        return 0
    elif allele == 1:
        return 2    
    elif allele == -1:
        return 1
    else:
        assert False

def idx2allele(idx):
    if idx == 0:
        return 0
    elif idx == 2:
        return 1
    elif idx == 1:
        return -1
    else:
        assert False
        
def hap2idx(hap):
    return tuple(map(allele2idx,hap))

def idx2hap(idx):
    return tuple(map(idx2allele, idx))

def config2array(config):
    ret = numpy.zeros((3,3),dtype=int)
    for hap,count in config:
        ret[hap2idx(hap)] = count
    return ret

def array2dict(array):
    ret = {}
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            if i == 1 and j == 1:
                continue
            ret[idx2hap((i,j))] = array[i,j]
    return ret
    
    

    

    
    
# class MoranDemoRates(object):
#     ## TODO: this class is useless now
#     def __init__(self, states, popSize, theta, mymap=None):
#         if mymap is None:
#             mymap = map
        
#         start = time.time()
        
#         self.theta = theta
#         self.popSize = popSize
#         self.states = states               
    
#     def getRates(self, rho):
#         return MoranRates(nonRecomRates=self, rho=rho)

'''
Created on Jan 19, 2015

@author: jkamm
'''
from __future__ import division
from __future__ import absolute_import
from builtins import map
from builtins import zip
from builtins import range
from builtins import object
from .compute_stationary import stationary1d_tridiagonal
from .compute_stationary import assertValidProbs
import logging
import time
import numpy
import scipy
import scipy.special
from scipy import sparse


# Have a separate class for the Rates
# so we don't have to pickle all of MoranStatesAugmented
# when doing multiprocessing
# (saves memory and communication time with worker processes)
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
        n = self.n
        coalRate = 1. / popSize
        recomRate = float(rho) / 2.

        if rho == 0.0:
            return numpy.array([0.0] * self.n + [1.0])
        else:
            numCoupledLinsRates = sparse.dok_matrix((n+1, n+1))
            for i in range(n+1):
                if i < n:
                    numCoupledLinsRates[i, i+1] = ((n-i)**2) * coalRate
                    numCoupledLinsRates[i, i] -= numCoupledLinsRates[i, i+1]
                if i > 0:
                    numCoupledLinsRates[i, i-1] = recomRate * i
                    numCoupledLinsRates[i, i] -= numCoupledLinsRates[i, i-1]
            return stationary1d_tridiagonal(numCoupledLinsRates)

    def getRates(self, popSize, theta, rho):
        start = time.time()
        recomRate = float(rho) / 2.
        mutRate = float(theta) / 2.
        coalRate = 1. / float(popSize)
        ret = (recomRate * self.unscaled_recom_rates
               + mutRate * self.unscaled_mut_rates
               + coalRate * self.unscaled_coal_rates)
        end = time.time()
        logging.info('%f seconds to construct rates for '
                     'rho=%f,theta=%f,N=%f' % (end-start, rho, theta, popSize))
        return ret


# make all haplotypes
a_haps = []
b_haps = []
c_haps = []
for allele1 in range(2):
    a_haps.append((allele1, -1))
    b_haps.append((-1, allele1))
    for allele2 in range(2):
        c_haps.append((allele1, allele2))

all_haps = a_haps + b_haps + c_haps


def makeAllConfigs(hapList, n):
    # make all configs
    # represent a config as a dict
    tmpConfigList = [{}]
    for hapIdx, hap in enumerate(hapList):
        newConfigList = []
        for config in tmpConfigList:
            numHaps = sum([v for k, v in config.items()])
            assert numHaps <= n
            if hapIdx == len(hapList)-1:
                next_count = [n-numHaps]
            else:
                next_count = range(n - numHaps + 1)
            for i in next_count:
                newConfig = dict(config)
                newConfig[hap] = i
                newConfigList.append(newConfig)
        tmpConfigList = newConfigList

    return tmpConfigList


def one_locus_probs(popSize, theta, n):
    coalRate = 1. / popSize
    mutRate = float(theta) / 2.

    numOnesRates = sparse.dok_matrix((n+1, n+1))
    for i in range(n+1):
        if i < n:
            numOnesRates[i, i+1] = (n-i) * mutRate + i * (n-i) / 2.0 * coalRate
            numOnesRates[i, i] -= numOnesRates[i, i+1]
        if i > 0:
            numOnesRates[i, i-1] = i * mutRate + i * (n-i) / 2.0 * coalRate
            numOnesRates[i, i] -= numOnesRates[i, i-1]

    return stationary1d_tridiagonal(numOnesRates)


class AbstractMoranStates(object):
    def __init__(self, n):
        self.n = n
        self._stationary = {}

    def build_all_configs(self, n, exact):
        '''
        Create self.config_array, defined by:
        self.config_array[i, a, b] = the count of haplotype
                                     (a,b) in the i-th config
        '''
        if exact:
            cList = list(range(n+1))
        else:
            cList = [n]

        aConfigs = {n-c: makeAllConfigs(a_haps, n-c) for c in cList}
        bConfigs = {n-c: makeAllConfigs(b_haps, n-c) for c in cList}
        cConfigs = {c: makeAllConfigs(c_haps, c) for c in cList}

        all_configs = []

        for numC in cList:
            for aConf in aConfigs[n - numC]:
                for bConf in bConfigs[n - numC]:
                    for cConf in cConfigs[numC]:
                        conf = {}
                        conf.update(aConf)
                        conf.update(bConf)
                        conf.update(cConf)

                        all_configs.append(conf)

        self.config_array = numpy.zeros((len(all_configs), 3, 3), dtype=int)
        for idx, conf in enumerate(all_configs):
            for (i, j), count in conf.items():
                self.config_array[idx, i, j] = count

        # create dictionary mapping their hash values back to their index
        hash_vals = self.hash_config_array(self.config_array)
        assert len(set(hash_vals)) == len(hash_vals)  # should be all unique
        self.hash_to_allIdx = {k: v for v, k in enumerate(hash_vals)}

    def hash_config_array(self, conf_arr):
        base = self.n+1
        hash_vals = (conf_arr[:, 0, 0]
                     + base * conf_arr[:, 0, 1]
                     + (base**2) * (conf_arr[:, 1, 0]))
        if self.exact:
            hash_vals += ((base**3)*(conf_arr[:, 1, 1])
                          + (base**4)*(conf_arr[:, 0, -1])
                          + (base**5)*(conf_arr[:, -1, 0]))
        return hash_vals

    def numOnes(self, loc):
        return self.folded_config_array.sum(axis=1+(1-loc))[:, 1]

    def hapCount(self, hap):
        return numpy.array(self.folded_config_array[:, hap[0], hap[1]])

    def getUnlinkedStationary(self, popSize, theta):
        one_loc_probs = one_locus_probs(popSize=popSize, theta=theta, n=self.n)
        assertValidProbs(one_loc_probs)

        n = self.n
        leftOnes = self.numOnes(0)
        rightOnes = self.numOnes(1)
        bothOnes = self.hapCount((1, 1))
        joint = one_loc_probs[leftOnes] * one_loc_probs[rightOnes]
        if self.exact:
            joint[self.numC > 0] = 0
        else:
            joint *= (scipy.special.comb(rightOnes, bothOnes)
                      * scipy.special.comb(n-rightOnes, leftOnes-bothOnes)
                      / scipy.special.comb(n, leftOnes))

        joint *= self.n_unfolded_versions

        assertValidProbs(joint)
        return joint

    def build_symmetries(self):
        start = time.time()

        # the index of the folded version in all_configs
        folded_list = get_folded_config_idxs(self)

        # foldedIdx = the index in folded_configs
        # allIdx = the index in all_configs
        foldedIdx_to_allIdx = numpy.array(list(set(folded_list)))

        allIdx_to_foldedIdx = {v: k for k, v in enumerate(foldedIdx_to_allIdx)}
        allIdx_to_foldedIdx = [allIdx_to_foldedIdx[x] for x in folded_list]

        self.hash_to_foldedIdx = {k: allIdx_to_foldedIdx[v]
                                  for k, v in self.hash_to_allIdx.items()}
        self.folded_config_array = self.config_array[foldedIdx_to_allIdx, :, :]

        self.numC = (self.folded_config_array[:, 0, 0]
                     + self.folded_config_array[:, 0, 1]
                     + self.folded_config_array[:, 1, 0]
                     + self.folded_config_array[:, 1, 1])

        symm_mat = sparse.dok_matrix((len(allIdx_to_foldedIdx),
                                      self.folded_config_array.shape[0]))
        for i, j in enumerate(allIdx_to_foldedIdx):
            symm_mat[i, j] = 1
        symm_mat = symm_mat.tocsc()

        antisymm_mat = symm_mat.transpose().tocsr(copy=True)
        # normalize rows
        self.n_unfolded_versions = numpy.array(antisymm_mat.sum(axis=1))[:, 0]
        row_indices, col_indices = antisymm_mat.nonzero()
        antisymm_mat.data /= self.n_unfolded_versions[row_indices]

        self.symmetries = symm_mat.tocsr()
        self.antisymmetries = antisymm_mat.tocsr()

        logging.info('%f seconds to build symmetry matrices'
                     % (time.time() - start))

    def ordered_log_likelihoods(self, liks):
        try:
            return {time: self.ordered_log_likelihoods(l)
                    for time, l in liks.items()}
        except AttributeError:
            liks = liks * self.antisymmetries

            all_nC = self.config_array[:, :-1, :-1].sum(axis=(1, 2))
            liks = liks[all_nC == self.n]

            full_confs = self.config_array[:, :-1, :-1][all_nC == self.n, :, :]

            liks = numpy.log(liks)
            liks -= scipy.special.gammaln(self.n+1)
            for i in (0, 1):
                for j in (0, 1):
                    liks += scipy.special.gammaln(full_confs[:, i, j]+1)

            full_confs = [tuple(sorted(((i, j), cnf[i, j])
                                       for i in (0, 1) for j in (0, 1)))
                          for cnf in full_confs]
            return dict(zip(full_confs, liks))


class MoranStatesAugmented(AbstractMoranStates):
    '''
    maintains a representation of the states(possible configs)
    of the 2 locus Moran model
    '''
    def __init__(self, n):
        '''
        Constructor
        '''
        start = time.time()
        super(MoranStatesAugmented, self).__init__(n)
        self.exact = True

        self.build_all_configs(n, exact=True)

        end = time.time()
        logging.info('Constructed exact states in %f seconds' % (end - start))

        self.build_symmetries()

        start = time.time()
        self.unscaled_recom_rates = build_recom_rates(self)
        logging.info('Constructed recombination rate matrix in %f seconds'
                     % (time.time() - start))

        start = time.time()
        self.unscaled_mut_rates = build_mut_rates(self)
        logging.info('Constructed mut rate matrix in %f seconds'
                     % (time.time() - start))

        start = time.time()
        self.unscaled_coal_rates = (build_copy_rates(self)
                                    + build_cross_coal_rates(self))
        logging.info('Constructed coalescent/copying rate matrix in %f seconds'
                     % (time.time() - start))


def get_folded_config_idxs(states):
    arr = states.config_array
    # move the missing allele in between alleles 0,1
    arr = arr[:, (0, -1, 1), :][:, :, (0, -1, 1)]

    # relabel alleles 0,1 (4 ways to do this)
    symm_arrs = [arr, arr[:, ::-1, :], arr[:, :, ::-1], arr[:, ::-1, ::-1]]
    # swap the 2 loci
    symm_arrs += [numpy.transpose(a, axes=(0, 2, 1)) for a in symm_arrs]

    # swap back allele 1 with missing allele
    symm_arrs = [a[:, (0, -1, 1), :][:, :, (0, -1, 1)] for a in symm_arrs]

    # get hash val for each (folded) config
    hash_vals = numpy.vstack(list(map(states.hash_config_array, symm_arrs)))
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
        rates = confs[:, hap[0], hap[1]]

        otherConfs = numpy.array(confs)
        otherConfs[:, hap[0], hap[1]] -= 1
        otherConfs[:, hap[0], -1] += 1
        otherConfs[:, -1, hap[1]] += 1

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
        rates = confs[:, hap[0], hap[1]]
        for loc in range(2):
            if hap[loc] == -1:
                continue
            otherHap = [hap[0], hap[1]]
            otherAllele = 1 - hap[loc]
            otherHap[loc] = otherAllele

            otherConfs = numpy.array(confs)
            otherConfs[:, hap[0], hap[1]] -= 1
            otherConfs[:, otherHap[0], otherHap[1]] += 1

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
            for loc in range(2):
                if hap[loc] == -1 and otherHap[loc] != -1:
                    canCopy = False
            if not canCopy:
                continue

            copiedHap = [hap[0], hap[1]]
            for loc in range(2):
                if otherHap[loc] == -1:
                    copiedHap[loc] = -1
            copiedHap = tuple(copiedHap)

            hapMissing = (hap[0] == -1) + (hap[1] == -1)
            otherMissing = (otherHap[0] == -1) + (otherHap[1] == -1)
            assert otherMissing >= hapMissing

            rates = (confs[:, hap[0], hap[1]]
                     * confs[:, otherHap[0], otherHap[1]] / 2.)
            if otherMissing > hapMissing:
                rates *= 2

            otherConfs = numpy.array(confs)
            otherConfs[:, otherHap[0], otherHap[1]] -= 1
            otherConfs[:, copiedHap[0], copiedHap[1]] += 1

            ret = ret + get_rates(states, otherConfs, rates)

    return subtract_rowsum_on_diag(ret)


def subtract_rowsum_on_diag(spmat):
    spmat = spmat.tocsr() - sparse.diags(numpy.array(spmat.sum(axis=1)).T,
                                         offsets=[0],
                                         format='csr')
    return spmat.tocsr()


def build_cross_coal_rates(states):
    assert states.exact
    ret = sparse.csr_matrix(tuple([states.folded_config_array.shape[0]]*2))

    confs = states.folded_config_array
    for hap in c_haps:
        otherConfs = numpy.array(confs)

        rates = otherConfs[:, hap[0], -1] * otherConfs[:, -1, hap[1]]

        otherConfs[:, hap[0], hap[1]] += 1
        otherConfs[:, hap[0], -1] -= 1
        otherConfs[:, -1, hap[1]] -= 1

        ret = ret + get_rates(states, otherConfs, rates)

    return subtract_rowsum_on_diag(ret)


def get_rates(states, otherConfs, rates):
    otherConfs = otherConfs[rates != 0, :, :]
    otherConfs = states.hash_config_array(otherConfs)
    otherConfs = numpy.array([states.hash_to_foldedIdx[x] for x in otherConfs],
                             dtype=int)

    confs = numpy.arange(states.folded_config_array.shape[0], dtype=int)
    confs = confs[rates != 0]

    rates = rates[rates != 0]

    ret = sparse.coo_matrix((rates, (confs, otherConfs)),
                            shape=[states.folded_config_array.shape[0]]*2)
    return ret.tocsr()

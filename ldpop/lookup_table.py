from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import map
from builtins import range
from builtins import object
from .compute_likelihoods import folded_likelihoods, NumericalError
from .moran_augmented import MoranStatesAugmented, MoranRates
from .moran_finite import MoranStatesFinite
from .compute_stationary import stationary

from multiprocessing import Pool
import logging
import time
import pandas
import numpy


def getKey(num00, num01, num10, num11):
    key = {}
    key[(0, 0)] = num00
    key[(0, 1)] = num01
    key[(1, 0)] = num10
    key[(1, 1)] = num11
    # key[(0,-1)] = 0
    # key[(-1,0)] = 0
    # key[(1,-1)] = 0
    # key[(-1,1)] = 0
    # return frozenset(key.items())
    return tuple(sorted(key.items()))


# columns has the columns of the table
# plus possibly extraneous things, but we just pull out what we want
def getRow(num00, num01, num10, num11, columns, rhos):
    toReturn = []
    for rho in rhos:
        key = getKey(num00, num01, num10, num11)
        toReturn.append(columns[rho][key])
    return toReturn


def get_states(n, exact):
    if exact:
        return MoranStatesAugmented(n)
    else:
        return MoranStatesFinite(n)


def getColumnHelper(args):
    return getColumn(*args)


def getColumn(moranRates, rho, theta, popSizes, timeLens, init):
    try:
        return folded_likelihoods(moranRates,
                                  rho,
                                  theta,
                                  popSizes,
                                  timeLens,
                                  lastEpochInit=init)
    except NumericalError as err:
        print(rho)
        print(err)


def computeLikelihoods(n, exact, popSizes, theta, timeLens, rhoGrid, cores,
                       store_stationary=None, load_stationary=None):
    rhoGrid = list(rhoGrid)
    assert rhoGrid == sorted(rhoGrid)

    # make the pool first to avoid copying large objects.
    # maxtasksperchild=1 to avoid memory issues
    executor = Pool(cores, maxtasksperchild=1)

    # make the states and the rates
    states = get_states(n, exact)
    moranRates = MoranRates(states)

    # compute initial distributions and likelihoods
    prevInit = states.getUnlinkedStationary(popSize=popSizes[-1], theta=theta)
    inits = []

    if load_stationary:
        stationary_dists = numpy.load(load_stationary)
        for stationary_dist in stationary_dists:
            inits.append(stationary_dist)
    else:
        for rho in reversed(rhoGrid):
            rates = moranRates.getRates(rho=rho,
                                        popSize=popSizes[-1],
                                        theta=theta)
            prevInit = stationary(Q=rates,
                                  init=prevInit,
                                  norm_order=float('inf'),
                                  epsilon=1e-2)
            inits.append(prevInit)
    ret = executor.map(getColumnHelper,
                       [(moranRates, rho, theta, popSizes, timeLens, prevInit)
                        for rho, prevInit in zip(reversed(rhoGrid), inits)])
    logging.info('Cleaning up results...')
    if store_stationary:
        full_inits = numpy.array([result[1] for result in ret])
        numpy.save(store_stationary, full_inits)
    ret = [states.ordered_log_likelihoods(result[0]) for result in ret]
    executor.close()

    return [(rho, lik) for rho, lik in zip(rhoGrid, reversed(ret))]


class LookupTable(object):
    '''
    Lookup table of 2-locus likelihoods. Construct with
        lookup_table = LookupTable(n, theta, rhos, [pop_sizes],
                                   [times], [exact], [processes])
    (optional arguments in square bracket [])

    Printing
    --------
    str(lookup_table) returns a string in the format expected by ldhat
    (or ldhelmet, if the rhos are not a uniformly spaced grid starting at 0)
    to print the table to STDOUT in this format, do
        print lookup_table

    Attributes
    ----------
    lookup_table.table = pandas.DataFrame,
        with rows corresponding to configs and columns corresponding to rho. So
             lookup_table.table.ix['13 0 0 1', 1.0]
        returns the likelihood of sample config '13 0 0 1' at rho=1.0.
    lookup_table.n = number of haplotypes
    lookup_table.theta = 2*mutation rate
    lookup_table.column = [rho0,rho1,...] = the grid of rhos (2*recomb rate)
    lookup_table.index = [config0,config1,...] = the sample configs
    lookup_table.pop_sizes = [size0,size1,...,sizeD]
                           = coalescent scaled population sizes
         size0 = present size, sizeD = ancient size
    lookup_table.times = [t1,...,tD] = the times of size changes,
                                       going backwards from present
         must be increasing positive reals
    lookup_table.exact = boolean
         if False, use finite moran model, a reasonable approximation
         that is much faster than the exact formula. Accuracy of the
         approximation can be improved by taking n larger than needed,
         and subsampling (e.g. using ldhat/lkgen.c). As n->infty,
         the finite moran model -> 'exact' diffusion.

    Parallelization: use
        LookupTable(...,processes)
    to specify the number of parallel processes to use.
    '''
    def __init__(self, n, theta, rhos, pop_sizes=[1], times=[], exact=True,
                 processes=1, store_stationary=None, load_stationary=None):
        assert (list(rhos) == list(sorted(rhos))
                and len(rhos) == len(set(rhos))), 'rhos must be sorted, unique'
        assert len(pop_sizes) == len(times)+1

        timeLens = epochTimesToIntervalLengths(times)
        assert len(pop_sizes) == len(timeLens)+1

        start = time.time()
        minRho = rhos[0]

        # only use exact to compute rho > 0
        if exact and minRho == 0.0:
            rhos = rhos[1:]
        results = computeLikelihoods(n, exact, pop_sizes, theta, timeLens,
                                     rhos, processes, store_stationary,
                                     load_stationary)

        # use approx to compute rho == 0.0
        # because exact==approx and approx is faster
        if exact and minRho == 0.0:
            logging.info('Computing column for rho=0')
            results = computeLikelihoods(n,
                                         False,
                                         pop_sizes,
                                         theta,
                                         timeLens,
                                         [0.0],
                                         processes) + results
            rhos = [0.0] + rhos

        halfn = int(n) // 2
        numConfigs = (1
                      + halfn + halfn * (halfn - 1) * (halfn + 4) // 6
                      + (halfn - 1) * (halfn + 2) // 2)

        columns = dict(results)
        index, rows = [], []
        # make all these configs then print them out
        for i in range(1, halfn + 1):
                for j in range(1, i + 1):
                    for k in range(j, -1, -1):
                        hapMult11 = k
                        hapMult10 = j - k
                        hapMult01 = i - k
                        hapMult00 = n - i - j + k

                        index += ['%d %d %d %d' % (hapMult00,
                                                   hapMult01,
                                                   hapMult10,
                                                   hapMult11)]
                        rows += [getRow(hapMult00,
                                        hapMult01,
                                        hapMult10,
                                        hapMult11,
                                        columns,
                                        rhos)]
        self.table = pandas.DataFrame(rows, index=index, columns=rhos)

        assert self.table.shape[0] == numConfigs
        end = time.time()
        logging.info('Computed lookup table in %f seconds ' % (end-start))

        self.n = n
        self.theta = theta
        self.pop_sizes = list(pop_sizes)
        self.times = list(times)
        self.exact = exact

    def __str__(self):
        ret = []
        ret += [[self.n, self.table.shape[0]]]
        ret += [[1, self.theta]]

        ret += [rhos_to_string(self.table.columns).split()]

        ret += [[], []]
        for i, (config, row) in enumerate(self.table.iterrows(), start=1):
            ret += [[i, '#', config, ':'] + list(row)]

        return '\n'.join([' '.join(map(str, x)) for x in ret])


def epochTimesToIntervalLengths(epochTimes):
    if len(epochTimes) == 0:
        return []
    if epochTimes[0] == 0:
        raise IOError('Your first epoch time point should not be zero!')
    epochLengths = list(epochTimes)
    totalTime = 0.
    for i in range(0, len(epochLengths)):
        epochLengths[i] = epochLengths[i] - totalTime
        totalTime += epochLengths[i]
    return epochLengths


def rhos_to_string(rhos):
    rhos = numpy.array(rhos)
    if rhos[0] == 0 and numpy.allclose(rhos[1:] - rhos[:-1],
                                       rhos[-1] / float(len(rhos)-1),
                                       atol=0):
        rho_line = [len(rhos), rhos[-1]]
    else:
        rho_line = []
        prev_rho, prev_diff = rhos[0], float('inf')
        for rho in rhos[1:]:
            if not numpy.isclose(rho - prev_rho, prev_diff, atol=0):
                prev_diff = rho - prev_rho
                rho_line += [prev_rho, prev_diff]
            prev_rho = rho
        rho_line += [prev_rho]
    return ' '.join(map(str, rho_line))


def rhos_from_string(rho_string):
    '''
    Return a list of rhos obtained from a comma
    separated string of rhos in one of two formats:
    if the rho_string has two elements, <num_rho>,<max_rho>
    return a list of size num_rho [0, ...., max_rho]
    otherwise, the rho_string should be <rho_1>,<step_size_1>,<rho_2>...
    and return [rho_1, rho_1 + step_size_1, ..., rho_2,...]
    note that this implies that if rho_string is just <rho>, return [rho].
    '''

    rho_args = rho_string.split(',')

    if len(rho_args) == 2:
        n, R = int(rho_args[0]), float(rho_args[1])
        return list(numpy.arange(n, dtype=float) / float(n-1) * R)

    rhos = [float(rho_args[0])]
    arg_idx = 1
    while(arg_idx < len(rho_args)):
        assert arg_idx + 1 < len(rho_args)
        step_size = float(rho_args[arg_idx])
        endingRho = float(rho_args[arg_idx+1])
        arg_idx += 2
        cur_rho = rhos[-1]

        # 1e-13 deals with the numeric issues
        while(cur_rho < endingRho-1e-13):
            cur_rho += step_size
            rhos.append(cur_rho)
        if abs(cur_rho - endingRho) > 1e-13:
            print(cur_rho)
            print(endingRho)
            raise IOError('the Rhos you input are not so nice'
                          '(stepsize should divide difference in rhos)')
    return rhos

from __future__ import division
from __future__ import absolute_import
from builtins import zip
from builtins import map
from builtins import str
from builtins import range
from builtins import object
from .lookup_table import epochTimesToIntervalLengths, rhos_to_string
from .compute_likelihoods import folded_likelihoods
from .moran_finite import MoranStatesFinite
from .moran_augmented import MoranRates

from multiprocessing import Pool
import math
import pandas
import logging
import time

haps = []
for a in range(2):
    for b in range(2):
        haps.append((a, b))


class ISProposal(object):
    '''
    Tensor of 2-locus likelihoods, with axes
    corresponding to recombination rates, times, and configs.
    This is to be used as a proposal distribution
    with ImportanceSampler.jar Construct with
    ldproposal = ISProposal(n, theta, rhos, pop_sizes,
                            times, grid_points, [processes])
    (optional arguments in square bracket [])

    ISProposal is a subclass of pandas.Panel,

    Printing
    --------
    str(ldproposal) returns a string in the format
    expected by ImportanceSampler.jar
    to print the table to STDOUT in this format, do
        print ldproposal

    Attributes
    ----------
    ldproposal.panel = pandas.Panel,
         with the index axis corresponding to rho,
         rows corresponding to configs
         and columns corresponding to time points. So
                ldproposal.ix[5.0,'13 0 0 1', 1.0]
         returns the approximate likelihood of sample config
         '13 0 0 1' with rho=5.0 at time 1.0.
    ldproposal.num_haps = number of haplotypes
        in each config in the table (2*n in the constructor)
    ldproposal.theta = 2*mutation rate
    ldproposal.pop_sizes = [size0,size1,...,sizeD]
                         = coalescent scaled population sizes
         size0 = present size, sizeD = ancient size
    ldproposal.times = [t1,...,tD] = the times of size changes
        going backwards from present. must be increasing positive reals
    The data can be accessed using the methods of pandas.Panel,
        e.g. ldproposal.ix[5.0,'13 0 0 1', 1.0] as above
    Parallelization: use
        ISProposal(...,processes)
    to specify the number of parallel processes to use.
    '''
    def __init__(self, n, theta, rhos, pop_sizes,
                 times, numTimePointsPerEpoch, processes=1):
        start = time.time()
        epochLengths = epochTimesToIntervalLengths(times)
        # All possible configs
        states = MoranStatesFinite(2*n)
        moranRates = MoranRates(states)

        executor = Pool(processes)
        likelihoodDictList = list(map(
            states.ordered_log_likelihoods,
            executor.map(ordered_wrapper, [(moranRates, rho, theta,
                                            pop_sizes, epochLengths,
                                            numTimePointsPerEpoch)
                                           for rho in rhos])))
        executor.close()
        executor.join()

        data = {}
        for rho, likelihoodDict in zip(rhos, likelihoodDictList):
            timeList = list(likelihoodDict.keys())
            timeList.sort()
            index, rows = [], []
            for config in sorted(likelihoodDict[0.0].keys()):
                config_dict = dict(config)
                index += [' '.join([str(config_dict[hap]) for hap in haps])]
                rows += [[math.exp(likelihoodDict[disc_time][config])
                          for disc_time in timeList]]
            data[rho] = pandas.DataFrame(rows, index=index, columns=timeList)

        self.panel = pandas.Panel(data)
        end = time.time()
        logging.info('Computed lookup table in %f seconds ' % (end-start))
        self.num_haps = 2*n
        self.theta = theta
        self.pop_sizes = pop_sizes
        self.times = times

    # TODO: fix w.r.t. dealing with pandas stuff above
    def __str__(self):
        ret = []
        # header
        ret += [['numHaps', int(self.num_haps/2)]]
        ret += [['theta', self.theta]]
        ret += [['popSizes', ','.join(map(str, self.pop_sizes))]]
        ret += [['epochTimes', ','.join(map(str, self.times))]]
        ret += [['rhos', rhos_to_string(list(self.panel.keys()))]]
        ret += [['%']]   # delimiter
        # Iterate over each config x timepoints table for a given rho
        # then print.
        for rho in self.panel.keys():
            data_frame = self.panel[rho]
            ret += [['rho', rho]]
            ret += [['config', ' '.join(map(str, list(data_frame.keys())))]]
            for (config, row) in data_frame.iterrows():
                ret += [[config, ':'] + list(row)]
            ret += [['$']]  # delimiter

        return '\n'.join([' '.join(map(str, x)) for x in ret])


def ordered_wrapper(args_list):
    (moranRates,
     rho,
     theta,
     popSizes,
     epochLengths,
     numTimePointsPerEpoch) = args_list
    return folded_likelihoods(moranRates, rho, theta, popSizes, epochLengths,
                              gridPointsPerEpoch=numTimePointsPerEpoch)[0]
# Prints grid of the form
# config:    L@t1    L@t_2    L@time3...

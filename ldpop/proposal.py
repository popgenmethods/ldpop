from lookup_table import epochTimesToIntervalLengths, rhos_from_string, rhos_to_string
from compute_likelihoods import folded_likelihoods
from moran_finite import MoranStatesFinite
from moran_augmented import MoranRates

from multiprocessing import Pool
import logging, math, sys
from cStringIO import StringIO

haps = []
for a in range(2):
    for b in range(2):
        haps.append((a,b))

## TODO: change ISProposal to a class usable from python, similar to LookupTable. subclass of pandas.Panel?
def ISProposal(n, theta, rhos, pop_sizes, times, numTimePointsPerEpoch, processes):
    f = StringIO()
    epochLengths = epochTimesToIntervalLengths(times)    
    # Allows us to create a table for only one config instead of all of the configs
    num_haps = 2 * n
    f.write("\t".join(["numHaps"] + [str(n)]) + "\n")
    # All possible configs
    
    states = MoranStatesFinite(num_haps)
    
    # Rate matrices for each epoch. Not sure if the overhead is actually worth this but whatever
    #demoRatesList = [states.getDemoRates(theta, popSize) for popSize in popSizes]
    
    # Tell the importance sampler what parameters you used to compute stuff. Parameter header.
    f.write("\t".join(["theta"] + [str(theta)]) + "\n")
    f.write("\t".join(["popSizes"] + [",".join(map(str,pop_sizes))]) + "\n")
    f.write("\t".join(["epochTimes"] + [",".join(map(str,times))]) + "\n")
    f.write("\t".join(["rhos"] + [rhos_to_string(rhos)]) + "\n")
    f.write("%\n") # Ending delimiter of parameters

    moranRates = MoranRates(states)
    
    executor = Pool(processes)
    
    likelihoodDictList = map(states.ordered_log_likelihoods,
                             executor.map(ordered_wrapper, [(moranRates, rho, theta, pop_sizes, epochLengths, numTimePointsPerEpoch) for rho in rhos]))
    
    executor.close()
    executor.join()
    
    # Iterate over each config x timepoints table for a given rho, then print.
    for rho, likelihoodDict in zip(rhos,likelihoodDictList):   
        timeList = likelihoodDict.keys()
        timeList.sort()
        
        f.write("\t".join(["rho"] + [str(rho)]) + "\n")
        f.write("\t".join(["config"] + [str(t) for t in timeList]) + "\n")
        
        for config in sorted(likelihoodDict[0.0].keys()):
            configDict = dict(config)
            toPrint = []
            configString = " ".join([str(configDict[hap]) for hap in haps])
            toPrint.append(configString)
            toPrint.append(":")
            toPrint = toPrint + [str(math.exp(likelihoodDict[time][config])) for time in timeList]
            f.write("\t".join(toPrint) + "\n")
        f.write("$" + "\n")
    
    #f.close()
    return f.getvalue()
    

        
def ordered_wrapper(args_list):
    moranRates, rho, theta, popSizes, epochLengths, numTimePointsPerEpoch = args_list
    return folded_likelihoods(moranRates, rho, theta, popSizes, epochLengths, gridPointsPerEpoch=numTimePointsPerEpoch)
# Prints grid of the form
# config:    L@t1    L@t_2    L@time3...
# 

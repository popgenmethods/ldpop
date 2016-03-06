'''
Created on Feb 9, 2015

@author: jeffreyspence
'''

from build_ld_hat_table import epochTimesToIntervalLengths, ldhelmet_to_rho_array
from compute_likelihoods import folded_likelihoods
from moran_model_dk import MoranStatesDK
from moran_model import MoranRates

from multiprocessing import Pool
import argparse, logging, math, sys


haps = []
for a in range(2):
    for b in range(2):
        haps.append((a,b))

def ordered_wrapper(args_list):
    moranRates, rho, theta, popSizes, epochLengths, numTimePointsPerEpoch = args_list
    return folded_likelihoods(moranRates, rho, theta, popSizes, epochLengths, gridPointsPerEpoch=numTimePointsPerEpoch)
# Prints grid of the form
# config:    L@t1    L@t_2    L@time3...
# 
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--numHaps", type=int)
    parser.add_argument("--config", type=str) # specified in the order (0,0), (0,1), (1,0), (1,1), (-1,0), (-1,1), (0,-1), (1,-1)
    parser.add_argument("--theta", type=float)
    parser.add_argument("--popSizes", type=str)
    parser.add_argument("--epochTimes", type=str)
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--ldHelmetRhos", type=str)
    parser.add_argument("--timePointsPerEpoch", type = int)
    parser.add_argument("--outputFile",type=str)
    parser.add_argument("--log", type=str)    
    args = parser.parse_args()

    if args.log == ".":
        logging.basicConfig(level=logging.INFO)
    elif args.log is not None:
        logging.basicConfig(filename=args.log, level=logging.INFO)
    
    # Allow user to save moranGrid if they don't want to just pipe it in
    if args.outputFile == None:
        f = sys.stdout
    else:
        f = open(args.outputFile, 'w+')
    
    rhos, rho_string, true_rhos = ldhelmet_to_rho_array(args.ldHelmetRhos)
    
    theta = args.theta
    
    assert (args.popSizes is None) == (args.epochTimes is None)
    if args.popSizes is None:
        popSizes = [1]
        epochLengths = []
    else:
        popSizes = [float(i) for i in args.popSizes.split(",")]
        epochLengths = epochTimesToIntervalLengths([float(i) for i in args.epochTimes.split(",")])
    assert len(popSizes) == len(epochLengths)+1
    
    
    numTimePointsPerEpoch = args.timePointsPerEpoch
    assert len(popSizes) == len(epochLengths)+1
    
    # Allows us to create a table for only one config instead of all of the configs
    if args.config != None and args.numHaps == None:
        config = [int(i) for i in args.config.split(",")]
        assert len(config) == 8
        numHaps = 2 * sum(config)
        f.write("\t".join(["only_config"] + [str(args.config)]) + "\n")
    elif args.config == None and args.numHaps != None:
        numHaps = 2 * args.numHaps  #need twice as many guys to deal with missing data in the importance sampler
        f.write("\t".join(["numHaps"] + [str(args.numHaps)]) + "\n")
    else:
        raise IOError("You can only specify either a single config or the number of haplotypes!")
    # All possible configs
    
    states = MoranStatesDK(numHaps)
    
    # Rate matrices for each epoch. Not sure if the overhead is actually worth this but whatever
    #demoRatesList = [states.getDemoRates(theta, popSize) for popSize in popSizes]
    
    # Tell the importance sampler what parameters you used to compute stuff. Parameter header.
    f.write("\t".join(["theta"] + [str(theta)]) + "\n")
    f.write("\t".join(["popSizes"] + [str(args.popSizes)]) + "\n")
    f.write("\t".join(["epochTimes"] + [str(args.epochTimes)]) + "\n")
    f.write("\t".join(["cores"] + [str(args.cores)]) + "\n")
    f.write("\t".join(["ldHelmetRhos"] + [true_rhos]) + "\n")
    f.write("%\n") # Ending delimiter of parameters

    moranRates = MoranRates(states)
    
    #TODO: fix the False thing in case people want exact likelihoods thru time
    executor = Pool(args.cores)
    likelihoodDictList = map(states.ordered_log_likelihoods,
                             executor.map(ordered_wrapper, [(moranRates, rho, theta, popSizes, epochLengths, numTimePointsPerEpoch) for rho in rhos]))
    
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
    
    f.close()
    

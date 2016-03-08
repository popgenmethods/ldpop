#! /usr/bin/env python

from ldpop import ISProposal, rhos_from_string
import logging, argparse

## TODO: this needs to be cleaned up!
## make the command line args similar to ldtable.py

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--numHaps", type=int)
    parser.add_argument("--theta", type=float)
    parser.add_argument("--popSizes", type=str)
    parser.add_argument("--epochTimes", type=str)
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--ldHelmetRhos", type=str)
    parser.add_argument("--timePointsPerEpoch", type = int)
    parser.add_argument("--log", type=str)    
    args = parser.parse_args()

    if args.log == ".":
        logging.basicConfig(level=logging.INFO)
    elif args.log is not None:
        logging.basicConfig(filename=args.log, level=logging.INFO)
    
    rhos = rhos_from_string(args.ldHelmetRhos)
    
    theta = args.theta
    
    assert (args.popSizes is None) == (args.epochTimes is None)
    if args.popSizes is None:
        popSizes = [1]
        times = []
    else:
        popSizes = [float(i) for i in args.popSizes.split(",")]
        times = map(float, args.epochTimes.split(","))

    assert len(popSizes) == len(times)+1
    times = [float(i) for i in args.epochTimes.split(",")]
    
    numTimePointsPerEpoch = args.timePointsPerEpoch
    
    
    print ISProposal(args.numHaps, theta, rhos, popSizes, times, numTimePointsPerEpoch, args.cores)

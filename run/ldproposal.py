#! /usr/bin/env python
from __future__ import print_function
from ldpop import ISProposal, rhos_from_string
import logging, argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Print a table of (approximate) likelihoods through time, use in conjunction with bin/ImportanceSampler.jar")
    parser.add_argument("-n", type = int, metavar = "N", required = True, help = "sample size (table will contain configs with 2N haplotypes)" )
    parser.add_argument("-th", type=float, metavar = "theta", required = True, help = "twice the mutation rate")
    parser.add_argument("-rh", type=str, metavar = "rho", required = True, help = "twice the recomb rate.  Alternatively, one can input num_rho,max_rho to make a grid with num_rh uniformly spaced points from 0 to max_rh, inclusive. (((Alternatively, to create a non-uniform grid, use '-rh r0,step0,r1,step1,r2,...rK'. This creates a grid {r0,r0+step0,r0+2*step0,...,r1,r1+step1,...,rK} similar to ldhelmet. Note that non-uniform grid is incompatible with vanilla ldhat.)))")
    parser.add_argument("-s", type=str, metavar="s0,s1,...,sD", help="coalescent scaled population sizes (s0=present size, sD=ancient size)")
    parser.add_argument("-t", type=str, metavar="t1,...,tD", help="times of size changes from present backwards. Must be increasing positive reals.")
    parser.add_argument("--cores", type=int, default=1, help="Number of parallel processes to use.")
    parser.add_argument("--grdpts", default = 5, metavar = "D", type = int, help = "Number of points per epoch to discretize time.  More points will results in a longer run time, but more efficient importance sampling.")
    parser.add_argument("--log", type=str, metavar="logfile", help="Log extra info to logfile. If logfile='.', logs to STDERR.")     
    args = parser.parse_args()

    if args.log == ".":
        logging.basicConfig(level=logging.INFO)
    elif args.log is not None:
        logging.basicConfig(filename=args.log, level=logging.INFO)
    
    rhos = rhos_from_string(args.rh)
    
    
    assert (args.s is None) == (args.t is None)
    if args.s is None:
        popSizes = [1.]
        times = []
    else:
        popSizes = [float(i) for i in args.s.split(",")]
        times = [float(t) for t in args.t.split(",")]

    assert len(popSizes) == len(times)+1    
    
    print(ISProposal(args.n, args.th, rhos, popSizes, times, args.grdpts, args.cores))

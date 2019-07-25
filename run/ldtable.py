#! /usr/bin/env python
from __future__ import print_function
from ldpop import LookupTable, rhos_from_string
import argparse
import logging

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Print a lookup table, as expected by ldhat or ldhelmet.')
    parser.add_argument('-n', type=int, metavar='N', required=True,
                        help='sample size')
    parser.add_argument('-th', type=float, required=True, metavar='theta',
                        help='twice the mutation rate')
    parser.add_argument('-rh', type=str,
                        required=True, metavar='num_rh,max_rh',
                        help='grid of rhos (twice the recomb rate). The grid '
                             'has num_rh uniformly spaced points from 0 to '
                             'max_rh, inclusive. (((Alternatively, to create '
                             'a non-uniform grid, use '
                             '"-rh r0,step0,r1,step1,r2,...rK". '
                             'This creates a grid '
                             '{r0,r0+step0,r0+2*step0,...,r1,r1+step1,...,rK} '
                             'similar to ldhelmet. Note that non-uniform grid '
                             'is incompatible with vanilla ldhat.)))')
    parser.add_argument('-s', type=str, metavar='s0,s1,...,sD',
                        help='coalescent scaled population sizes '
                             '(s0=present size, sD=ancient size)')
    parser.add_argument('-t', type=str, metavar='t1,...,tD',
                        help='times of size changes from present backwards. '
                             'Must be increasing positive reals.')
    parser.add_argument('--approx', action='store_true',
                        help='use finite moran model. A reasonable '
                             'approximation that is much faster than '
                             'the exact formula. Accuracy of the  '
                             'approximation can be improved by taking n '
                             'larger than needed, and using ldhat/lkgen.c '
                             'to subsample. (Converges to the "exact" '
                             'diffusion as n->infty)')
    parser.add_argument('--cores', type=int, default=1,
                        help='Number of parallel processes to use.')
    parser.add_argument('--log', type=str, metavar='logfile',
                        help='Log extra info to logfile. '
                             'If logfile=".", logs to STDERR.')

    args = parser.parse_args()

    if args.log == '.':
        logging.basicConfig(level=logging.INFO)
    elif args.log is not None:
        logging.basicConfig(filename=args.log, level=logging.INFO)

    exact = not args.approx
    numCores = args.cores

    rhos = rhos_from_string(args.rh)

    popSizes, times = args.s, args.t
    assert (popSizes is None) == (times is None)
    if popSizes is None:
        popSizes = [1.]
        times = []
    else:
        popSizes = [float(_) for _ in popSizes.split(',')]
        times = [float(_) for _ in times.split(',')]

    assert len(popSizes) == len(times)+1

    print(LookupTable(args.n, args.th, rhos, popSizes, times, exact, numCores))

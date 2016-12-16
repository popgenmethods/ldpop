from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
import pytest
import ldpop, sys, re
import pickle as pickle
import numpy as np


example_ldtable_f = 'example_ldtable.txt'
def example_ldtable():
    return ldpop.LookupTable(n=5, theta=.001, rhos=[0,2,4,6,8,10],
                             pop_sizes=[100,.1,1], times=[.5,.58], processes=2)

def test_ldtable():
    # read in stored table
    with open(example_ldtable_f,'r') as f:
        lines = [l.strip() for l in f]
    # skip first 3 lines
    lines = lines[3:]
    # skip blank lines
    lines = [l for l in lines if l]
    # extract configs, likelihoods from each line
    line_re = re.compile(r"\d+ # (\d+ \d+ \d+ \d+) : (.+)")
    lines = [line_re.match(l) for l in lines]
    lines = [(l.group(1), l.group(2)) for l in lines]
    configs, likelihoods = zip(*lines)
    configs = list(configs)
    likelihoods = np.array([[float(l) for l in row.split()]
                            for row in likelihoods])

    # check that example_ldtable() matches the stored table
    ldtable = example_ldtable()
    assert np.allclose(ldtable.table.values, likelihoods)
    assert list(ldtable.table.index) == configs

example_ldproposal_f = 'example_ldproposal.txt'
def example_ldproposal():
    return ldpop.ISProposal(n=5, theta=.001, rhos=[0,2,4,6,8,10],
                            pop_sizes=[100,.1,1], times=[.5,.58],
                            numTimePointsPerEpoch=3)

def test_ldproposal():
    with open(example_ldproposal_f, 'r') as f:
        lines = [l.strip() for l in f]
    # skip first 5 lines
    lines = lines[5:]
    assert lines[0] == "%"
    line_re = re.compile(r"(\d+ \d+ \d+ \d+) : (.+)")
    panel = []
    for line in lines:
        if line.startswith("rho") or line.startswith("config"):
            continue
        elif line == '%' or line == '$':
            panel.append([])
        else:
            line = line_re.match(line)
            config = line.group(1)
            likelihoods = [float(l) for l in line.group(2).split()]
            panel[-1].append(likelihoods)
    assert len(panel[-1]) == 0
    del panel[-1]

    panel = np.array(panel)
    ldproposal = example_ldproposal()
    assert np.allclose(panel, ldproposal.panel.values)

if __name__=="__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "generate":
        with open(example_ldtable_f,'w') as f:
            print(example_ldtable(), file=f)
        with open(example_ldproposal_f,'w') as f:
            print(example_ldproposal(), file=f)
    elif len(sys.argv) == 1:
        pass
    else:
        raise Exception("Unrecognized command line options.")

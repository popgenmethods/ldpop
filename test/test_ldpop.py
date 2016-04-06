import pytest
import ldpop, sys
import cPickle as pickle
import numpy as np


example_ldtable_f = 'example_ldtable.pickle'
def example_ldtable():
    return ldpop.LookupTable(n=5, theta=.001, rhos=[0,2,4,6,8,10],
                             pop_sizes=[100,.1,1], times=[.5,.58], processes=2)

def test_ldtable():
    with open(example_ldtable_f,'r') as f:
        assert np.allclose(example_ldtable().table.values,
                           pickle.load(f).table.values)


example_ldproposal_f = 'example_ldproposal.pickle'
def example_ldproposal():
    return ldpop.ISProposal(n=5, theta=.001, rhos=[0,2,4,6,8,10],
                            pop_sizes=[100,.1,1], times=[.5,.58],
                            numTimePointsPerEpoch=3)

def test_ldproposal():
    with open(example_ldproposal_f,'r') as f:
        assert np.allclose(example_ldproposal().panel.values,
                           pickle.load(f).panel.values)


if __name__=="__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "generate":
        with open(example_ldtable_f,'w') as f:
            pickle.dump(example_ldtable(), f)
        with open(example_ldproposal_f,'w') as f:
            pickle.dump(example_ldproposal(), f)
    elif len(sys.argv) == 1:
        pass
    else:
        raise Exception("Unrecognized command line options.")
    

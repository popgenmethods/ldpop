#python build_ld_hat_table.py --numHaps 5 --theta .001 --popSizes 100,.1,1 --epochTimes .5,.58 --ldHelmetRhos 0,2,10 > test_example.txt
#python build_timed_likelihood_grid.py --numHaps 5 --theta .001 --popSizes 100,.1,1 --epochTimes .5,.58 --ldHelmetRhos 0,2,10 --timePointsPerEpoch 3 > test_example2.txt
#python -m pdb build_ld_hat_table.py --numHaps 5 --theta .001 --popSizes 100,.1,1 --epochTimes .5,.58 --ldHelmetRhos 0,2,10
python ../ldpop/build_ld_hat_table.py --numHaps 5 --theta .001 --popSizes 100,.1,1 --epochTimes .5,.58 --ldHelmetRhos 0,2,10 --cores 2 --log . | diff test_example.txt -
#python build_ld_hat_table.py --numHaps 5 --theta .001 --ldHelmetRhos 0,2,10 --cores 2
#python build_ld_hat_table.py --numHaps 5 --theta .001 --popSizes 100,.1,1 --epochTimes .5,.58 --ldHelmetRhos 0,2,10 --cores 2 --approx | diff test_example.txt -
#python -m pdb build_timed_likelihood_grid.py --numHaps 5 --theta .001 --popSizes 100,.1,1 --epochTimes .5,.58 --ldHelmetRhos 0,2,10 --timePointsPerEpoch 3
#python ../ldpop/build_timed_likelihood_grid.py --numHaps 5 --theta .001 --popSizes 100,.1,1 --epochTimes .5,.58 --ldHelmetRhos 0,2,10 --timePointsPerEpoch 3 --log . | diff test_example2.txt -
#python build_timed_likelihood_grid.py --numHaps 5 --theta .001 --ldHelmetRhos 0,2,10 --timePointsPerEpoch 3

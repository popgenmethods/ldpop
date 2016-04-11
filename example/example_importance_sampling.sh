#this example will generate a proposal distribution, and then use it to sample ten 2-locus ARGs for this configuration: 00:1, 01:2, 10:1, 11:1, 
#with rho = 1.0, and the population history described in the accompanying paper

#first, make a proposal distribution using run/ldproposal.py
../run/ldproposal.py -n 5 -th .008 -s 100,.1,1 -t .5,.58 -rh 1.0 --grdpts 3 --log . > proposal.txt

#now, use the proposal to print ARGs
java -ea -jar ../run/ImportanceSampler.jar --ldproposal proposal.txt --num_samples 10 --config 1,2,1,1 --trees --seed 995539384 > args.txt

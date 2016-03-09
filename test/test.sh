## Some very simple tests to check the output of ldtable.py and ldproposal.py

#../bin/ldtable.py -n 5 -th .001 -s 100,.1,1 -t .5,.58 -rh 6,10 --cores 2 --log . > example_ldtable.txt
../bin/ldtable.py -n 5 -th .001 -s 100,.1,1 -t .5,.58 -rh 6,10 --cores 2 --log . | diff example_ldtable.txt -
#python -m pdb ../bin/ldtable.py -n 5 -th .001 -s 100,.1,1 -t .5,.58 -rh 6,10 --cores 2 --log .


#../bin/ldproposal.py -n 5 -th .001 -s 100,.1,1 -t .5,.58 -rh 0,2,10 --grdpts 3 > example_ldproposal.txt
../bin/ldproposal.py -n 5 -th .001 -s 100,.1,1 -t .5,.58 -rh 0,2,10 --grdpts 3 --log . | diff example_ldproposal.txt -

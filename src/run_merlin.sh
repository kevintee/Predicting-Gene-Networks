# Runs MERLIN and places output in folder results/
# There is a seg fault, but some results are inputted into results.


mkdir -p results;
./gpdream/modules/Merlin/src/merlin -d data/BRCA_reformatted.txt -l data/tfs.txt -o results

# Runs MERLIN and places output in folder results/
# There is a seg fault, but some results are inputted into results.

cd gpdream/modules/Merlin/gsl-1.9;
./configure;
make;
sudo make install;
cd ../src;
make;
mkdir -p results;
./merlin -d example/net1_expression.txt -l example/net1_transcription_factors.tsv -o results/;

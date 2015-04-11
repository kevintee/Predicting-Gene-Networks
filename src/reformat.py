# Reformats data in data/ directory and saves .txt file in same directory with format described below.
import parse
import pdb
"""
	every gene expression data file is in the form:
	Name exp1 exp2 ...
	G1 valG1 valG1 valG1 ...
	G2 valG2 valG2 valG2 ...

	We want to convert that into the following format:
	G1     G2      G3    ...
	valG1  valG2   valG3 ...
	valG1  valG2   valG3 ... 
	...

"""


def main():
	# writes into data/<cancer_type>_reformatted.txt
	parse.reformat_data()

if __name__ == '__main__':
	main()
# Parses TCGA data from `data/` directory

import numpy as np

NUM_GENES = 8499

directory = 'data/'
files = ['BRCA.txt', 'COAD.txt', 'KIRC.txt', 'LUSC.txt', 'OV.txt', 'UCEC.txt']
reg_file = 'tfs.txt'

# X is in the form: Genes vs Tumor Sample
# regulators is a list of regulators
def parse_data():
    X = [[] for x in xrange(NUM_GENES)]
    for fname in files:
        with open(directory + fname) as f:
            for i,line in enumerate(f):
                # Skip title
                if i == 0:
                    continue
                vals = line.strip().split('\t') # Parse tsv
                vals = vals[1:] # Remove name
                vals = [float(x) for x in vals] # Convert to floats
                X[i-1].extend(vals)

    regulators = []
    with open(directory + reg_file) as f:
        regulators = [x.strip() for x in f]

    return np.asarray(X), regulators

def main():
    pass

if __name__ == '__main__':
    main()

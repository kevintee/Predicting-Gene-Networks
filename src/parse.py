# Parses TCGA data from `data/` directory

import numpy as np
import pdb

NUM_GENES = 8499

directory = 'data/'
files = ['BRCA.txt', 'COAD.txt', 'KIRC.txt', 'LUSC.txt', 'OV.txt', 'UCEC.txt']
reg_file = 'tfs.txt'

# X is in the form: Genes vs Tumor Sample
# regulators is a list of regulators
def parse_data():
    gene_data = {} # Dictionary of cancer to matrix
    for fname in files:
        with open(directory + fname) as f:
            X = []
            for i,line in enumerate(f):
                # Skip title
                if i == 0:
                    continue
                vals = line.strip().split('\t') # Parse tsv
                vals = vals[1:] # Remove name
                vals = [float(x) for x in vals] # Convert to floats
                X.append(vals)
            gene_data[fname.strip()[:-4]] = np.asarray(X)

    regulators = []
    with open(directory + reg_file) as f:
        regulators = [x.strip() for x in f]

    return gene_data, regulators


### for example converts BRCA.txt to BRCA_reformatted.txt New format
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

def reformat_data():
    gene_data = {} # Dictionary of cancer to matrix
    for fname in files:
        with open(directory + fname) as f:
            X = []
            gene_names = []
            for i,line in enumerate(f):
                # Skip title
                if i == 0:
                    continue
                vals = line.strip().split('\t') # Parse tsv
                gene_names.append(vals[0])
                vals = vals[1:] # Remove name
                vals = [float(x) for x in vals] # Convert to floats
                X.append(vals)
            label = fname.strip()[:-4]
            gene_data[label] = np.asarray(X)
            first_string = '' + gene_names[0]

            with open(directory + fname[:-4] + '_reformatted.txt', 'wb') as f_out:
                for gene_name in gene_names[1:]:
                    first_string += '\t%s' %(gene_name)
                f_out.write(first_string + '\n')
                next_string = ''
                height, width = gene_data[label].shape
                for i in range(width):
                    exp_string = ''
                    for j in range(height):
                        exp_string += '%s\t' %(str(2.0**(gene_data[label][j, i])))
                    f_out.write(exp_string[:-1] + '\n')

"""
This method parses modules.txt and outputs in a good datastructure described below.
"""

def parse_module(fname):
	f = open(fname, 'rb')
	module_id_to_genes = {}
	gene_to_module_id = {}
	content = f.read()
	f.close()
	content_lst = content.split('\n')

	modules_seen = set()
	for line in content_lst:
		if not line:
			break
		line_split = line.split('\t')
		gene, module_id = line_split[0], line_split[1]
		if not module_id in modules_seen:
			modules_seen.add(module_id)
			module_id_to_genes[module_id] = set()
		module_id_to_genes[module_id].add(gene)
		gene_to_module_id[gene] = module_id

	return module_id_to_genes, gene_to_module_id

def main():
    pass

if __name__ == '__main__':
    main()

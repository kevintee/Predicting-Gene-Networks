import pdb

"""
This file parses modules.txt and outputs in a good datastructure described below.
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

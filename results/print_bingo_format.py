import pdb

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
	module_id_to_genes, _ = parse_module('sbm/UCEC_cluster.txt')
	for cluster_id, genes in module_id_to_genes.iteritems():
		if len(genes) < 5:
			continue
		print 'cluster_%s' %(cluster_id)
		for gene in genes:
			print gene
		print 'batch'


if __name__ == '__main__':
	main()
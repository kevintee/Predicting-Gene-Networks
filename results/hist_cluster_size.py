import pdb
import matplotlib
import matplotlib.pyplot as pyplot

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
	module_id_to_genes, _ = parse_module('fold14/modules.txt')
	data = []
	# module_id_to_genes, _ = parse_module('sbm/UCEC_cluster.txt')
	for cluster_id, genes in module_id_to_genes.iteritems():
		data.append(len(genes))

	pdb.set_trace()
	measurements = data
	graph_minimum = min(data)
	graph_maximum = max(data)
	fig = pyplot.figure()
	ax = fig.add_subplot(1,1,1,)
	n, bins, patches = ax.hist(measurements, bins=50, range=(graph_minimum, graph_maximum), histtype='bar')

	#ax.set_xticklabels([n], rotation='vertical')

	for patch in patches:
	    patch.set_facecolor('r')

	pyplot.title('Distribution of cluster size')
	pyplot.xlabel('Cluster size')
	pyplot.ylabel('')
	pyplot.show()


if __name__ == '__main__':
	main()
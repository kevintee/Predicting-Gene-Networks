import networkx as nx
import pdb


def main():
	g = nx.read_weighted_edgelist('../results/Merlin/prediction_k300.txt',
						          delimiter='\t')
	pdb.set_trace()

if __name__ == '__main__':
	main()
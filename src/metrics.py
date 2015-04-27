# Add your metrics here

import scipy.stats
from parse import parse_sbm_results, parse_chip_tf, parse_all_sbm

key_tf = 'TFAP2A'

def score_gene_weights(test_vals, true_vals):
    intersect = set(test_vals) & set(true_vals)
    test_unique = set(test_vals) - intersect
    true_unique = set(true_vals) - intersect
    #return scipy.stats.hypergeom.cdf

def evaluate_network():
    tfs, genes, sbm_binary_matrix, sbm_module_to_gene, \
        sbm_gene_to_module = parse_sbm_results()

    sbm_results = parse_all_sbm(genes, sbm_binary_matrix, key_tf)
    chip_results = parse_chip_tf(key_tf)
    return score_gene_weights(sbm_results, chip_results)

def main():
    evaluate_network()

if __name__ == '__main__':
    main()

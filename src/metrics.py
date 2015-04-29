# Add your metrics here

import scipy.stats
from itertools import chain
from parse import parse_sbm_results, parse_chip_seq

def score_gene_weights(test_vals, true_vals):
    intersect = set(test_vals) & set(true_vals)
    test_unique = set(test_vals) - intersect
    true_unique = set(true_vals) - intersect
    #return scipy.stats.hypergeom.cdf

# This method makes the 'universe' of both the SBM and ChIP-seq results the same
def remove_unique_genes(sbm_results, chip_results):
    # First get the intersection of the two sets
    overlapping_genes = set(chain(*sbm_results.values()))
    overlapping_genes &= set(chain(*chip_results.values()))

    # Remove all results that aren't in either of the two sets
    sbm_modified = {}
    chip_modified = {}
    keys = set(sbm_results.keys()) & set(chip_results.keys())
    for k,v in sbm_results.items():
        if k in overlapping_genes and k in keys:
            v = list(set(v) & overlapping_genes)
            sbm_modified[k] = v

    for k,v in chip_results.items():
        if k in overlapping_genes and k in keys:
            v = list(set(v) & overlapping_genes)
            chip_modified[k] = v

    return sbm_modified, chip_modified

def evaluate_network():
    tfs, genes, sbm_results, sbm_module_to_gene, \
        sbm_gene_to_module = parse_sbm_results()

    chip_results = parse_chip_seq()
    sbm_results, chip_results = remove_unique_genes(sbm_results, chip_results)
    return score_gene_weights(sbm_results, chip_results)

def main():
    evaluate_network()

if __name__ == '__main__':
    main()

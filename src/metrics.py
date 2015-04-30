# Add your metrics here

import math
import scipy.stats
from itertools import chain
from parse import parse_sbm_results, parse_chip_seq

def score_gene_weights(test_vals, true_vals, genes):
    keys = set(test_vals.keys()) & set(true_vals.keys())
    num_genes = len(genes)
    p_vals = []
    for k in keys:
        test_val = test_vals[k]
        true_val = true_vals[k]
        intersection = len(set(test_val) & set(true_val))
        p_val = scipy.stats.hypergeom.sf(intersection, num_genes,
                len(true_val), len(test_val))
        if not math.isnan(p_val):
            p_vals.append(p_val)
    return p_vals

# This method makes the 'universe' of both the SBM and ChIP-seq results the same
def remove_unique_genes(sbm_results, chip_results):
    # First get the intersection of the two sets
    overlapping_genes = set(chain(*sbm_results.values()))
    overlapping_genes &= set(chain(*chip_results.values()))

    # Remove all results that aren't in either of the two sets
    sbm_modified = {}
    chip_modified = {}
    for k,v in sbm_results.items():
        if k in overlapping_genes:
            v = list(set(v) & overlapping_genes)
            sbm_modified[k] = v

    for k,v in chip_results.items():
        if k in overlapping_genes:
            v = list(set(v) & overlapping_genes)
            chip_modified[k] = v

    return sbm_modified, chip_modified, overlapping_genes

def evaluate_network():
    tfs, genes, sbm_results, sbm_module_to_gene, \
            sbm_gene_to_module = parse_sbm_results()

    chip_results = parse_chip_seq()
    sbm_results, chip_results, genes = \
            remove_unique_genes(sbm_results, chip_results)

    return score_gene_weights(sbm_results, chip_results, genes)

def main():
    print evaluate_network()

if __name__ == '__main__':
    main()

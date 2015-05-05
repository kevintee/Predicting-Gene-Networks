# Add your metrics here

import pdb
import math
import scipy.stats
from itertools import chain
from parse import parse_sbm_results, parse_chip_seq, parse_merlin
from plots import plot_p_vals

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

def bin_results(p_vals):
    bins = [0]*10
    for p_val in p_vals:
        i = 0
        while p_val < 1 and p_val != 0:
            i += 1
            p_val *= 10
        bins[i] += 1
    return bins

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
    merlin_results = parse_merlin()

    # Compare to SBM
    sbm_results, chip_results, genes = \
            remove_unique_genes(sbm_results, chip_results)

    p_vals = score_gene_weights(sbm_results, chip_results, genes)
    plot_p_vals(p_vals, 'Stochastic Block Model')

    # Compare to MERLIN
    merlin_results, chip_results, genes = \
            remove_unique_genes(merlin_results, chip_results)

    p_vals = score_gene_weights(merlin_results, chip_results, genes)
    plot_p_vals(p_vals, 'MERLIN')

def main():
    evaluate_network()

if __name__ == '__main__':
    main()

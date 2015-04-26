# Add your metrics here

from parse import parse_sbm_results, parse_chip_tf, parse_all_sbm

def score_gene_weights(test_vals, true_vals):
    pass

def main():
    key_tf = 'TFAP2A'
    tfs, genes, sbm_binary_matrix, sbm_module_to_gene, \
        sbm_gene_to_module = parse_sbm_results()

    sbm_results = parse_all_sbm(genes, sbm_binary_matrix, key_tf)
    chip_results = parse_chip_tf(key_tf)
    print sbm_results, chip_results

if __name__ == '__main__':
    main()

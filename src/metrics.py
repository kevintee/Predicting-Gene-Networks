# Add your metrics here

from parse import parse_sbm_results, parse_chip_tf

def sample_metric(X):
    pass

def main():
    tfs, genes, sbm_binary_matrix, sbm_module_to_gene, \
        sbm_gene_to_module = parse_sbm_results()

    chip_regulation = parse_chip_tf('TFAP2A')
    print chip_regulation

if __name__ == '__main__':
    main()

# Add your experiments here

from parse import parse_data

def sample_experiment(X):
    pass

def main():
    gene_data, regulators = parse_data()
    X = gene_data['BRCA']

if __name__ == '__main__':
    main()

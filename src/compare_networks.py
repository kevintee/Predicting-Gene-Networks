import networkx as nx
import pdb
import os
import argparse
import sys

""" Compare two networks (without modules)
    
    - Look @ the differences/similarities
    Measures:

    Compare two networks (with modules)


    Measures:
    - Look @ the differences/similarities


"""
def compare_two_networks(g0, g1): # Without module information


    return None


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', dest='input_files', nargs=2, help="enter atleast 2 input files.",
                        default=['../results/Merlin/prediction_k300.txt', '../results/Merlin/prediction_k300.txt'])
    parser.add_argument('-o', dest='output_file', nargs=1, help="enter output file name",
                        default=['log.txt'])
    # parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
    #                     default=sys.stdout)
    args = parser.parse_args()
    print args.input_files
    print args.output_file
    net0_fname = args.input_files[0]
    net1_fname = args.input_files[1]

    g0 = nx.read_weighted_edgelist(net0_fname, delimiter='\t')
    g1 = nx.read_weighted_edgelist(net1_fname, delimiter='\t')

    compare_two_networks(g0, g1)


if __name__ == '__main__':
    main()
import networkx as nx
import pdb
import os
import argparse
import sys


import parse_modules
import graph_io

""" Compare two networks (without modules)

    - Look @ the differences/similarities
    Measures:

    Compare two networks (with modules)

"""
def get_top_n_modules(m_to_gene, n):
    "return top n largest modules"
    d_sorted = sorted(m_to_gene.keys(), key=lambda x: len(m_to_gene[x]))
    d_sorted.reverse()
    return d_sorted[:n]

def compare_modules(m_g0, g_m0, m_g1, g_m1, g0, g1):
    """
        currently saves two images of the biggest module for graph0 and graph1.
    """
    print 'calling compare_modules'
    top10_0 = get_top_n_modules(m_g0, 10)
    top10_1 = get_top_n_modules(m_g1, 10)

    top_module_genes0 = list(m_g0[top10_0[0]])
    top_module_genes1 = list(m_g1[top10_1[0]])

    H0 = g0.subgraph(top_module_genes0)
    H1 = g1.subgraph(top_module_genes1)

    print top_module_genes0
    print top_module_genes1
    print top_module_genes0 == top_module_genes1

    print 'saving the two biggest modoles'

    #todo get this to display the same graph, since edgelist is the same?
    graph_io.save_graph(H0, '../results/Merlin/h0.png')
    graph_io.save_graph(H1, '../results/Merlin/h1.png')


    return None


def compare_two_networks(g0, g1): # Without module information


    return None


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', dest='input_files', nargs=2, help="enter atleast 2 input files.",
                        default=['../results/Merlin/prediction_k300.txt', '../results/Merlin/prediction_k300.txt'])
    parser.add_argument('-m', dest='input_modules', nargs=2, help="enter 2 input files corresponding to modules",
                         default=['../results/Merlin/modules.txt', '../results/Merlin/modules.txt'])
    parser.add_argument('-o', dest='output_file', nargs=1, help="enter output file name",
                        default=['log.txt'])

    args = parser.parse_args()
    print args.input_files
    print args.output_file
    print args.input_modules
    net0_fname = args.input_files[0]
    mod0_fname = args.input_modules[0]
    net1_fname = args.input_files[1]
    mod1_fname = args.input_modules[1]

    g0 = nx.read_weighted_edgelist(net0_fname, delimiter='\t')
    g1 = nx.read_weighted_edgelist(net1_fname, delimiter='\t')

    m_g0, g_m0 = parse_modules.parse_module(mod0_fname)
    m_g1, g_m1 = parse_modules.parse_module(mod1_fname)


    compare_modules(m_g0, g_m0, m_g1, g_m1, g0, g1)
    #pdb.set_trace()
    #compare_two_networks(g0, g1)


if __name__ == '__main__':
    main()

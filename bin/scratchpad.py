#!/usr/bin/env pypy
#
#   Scratchpad for working with implementation from the root folder.
#

from matplotlib.pylab import plot, ylabel, xlabel, savefig, close, title, figtext
from preprocessing.discretization import *
from preprocessing.preprocessors import remove_empty_samples, compute_relative_values, discrete_dataset_cleaning
import os
from plots.mtv_results import read_run_results
from plots.mtv_results import plot_run_results

# from plots.bacteria_histogram import plot_bacteria_hist
# plot_bacteria_hist('../plots/hist/normalized_depth_6/', depth=6)
# exit()

# plot_run_results('../experiments/4/')
# plot_run_results('../experiments/2a/')
# plot_run_results('../experiments/backup/')
# exit()
#
# from plots.faust_result_discretized import plot_faust_relationships
# plot_faust_relationships()
# exit()


def write_dataset_to_experiment(file_name, ds):
    from itemsets import binary_vectors_to_ints
    from utils.files import write_dat_file
    import pickle
    import itemsets
    # Write .dat file
    abundance = abundance_matrix(ds)
    D = binary_vectors_to_ints(abundance)
    write_dat_file(file_name + '.dat', D)

    # Create a header file
    headers = []
    for header in ds[0][2:]:
        vals = header.split('|')
        if len(vals) > 1:
            headers.append('|'.join(vals[-2:]))
        else:
            headers.append(vals[0])

    with open(file_name + '.headers', 'wb') as fd:
        line = ' '.join(headers)
        fd.write(line)

    # save the raw dataset
    with open(file_name + '.pickle', 'wb') as fd:
        pickle.dump(ds, fd)


def run_discretization_all_nodes():
    """
    TODO: Work in progress. Code to use a phylogenetic tree, and get a
    daset for the entire phylogenetic tree
    """
    from preprocessing import parser
    from preprocessing.tree import Tree
    from utils.dataset_helpers import abundance_matrix

    # Get the stool dataset and discretize it
    ds = parser.get_dataset()
    ds = compute_relative_values(ds)
    t = Tree(ds)
    print 'tree nodes before: ', t.count_nodes()
    print 'tree leafs before: ', t.count_leafs()
    ds = t.dataset_for_all_nodes()
    ds = median_discretization(ds)
    ds = discrete_dataset_cleaning(ds, 0.20)

    print 'Final number of attributes: ', len(ds[0]) - 2

    write_dataset_to_experiment('../experiments/1/Stool_maxent_discretized_all_nodes_020', ds)


# run_discretization_all_nodes()


def run_discretization_faust_nodes_and_leafs():
    """
    TODO: This currently result in only 37 attributes, and 32 after cleaning..
    Is this enough? Can we search the tree differently?
    """
    from preprocessing import parser
    from preprocessing.tree import Tree
    from utils.dataset_helpers import abundance_matrix
    from itemsets import binary_vectors_to_ints
    from utils.files import write_dat_file
    from preprocessing import faust_parser
    import itemsets

    # Get the Faust result clade names
    faust_results = faust_parser.results()
    faust_clades = []
    for faust_result in faust_results:
        faust_clades.append(faust_result.clade_1)
        faust_clades.append(faust_result.clade_2)

    # Get the stool dataset and discretize it
    ds = parser.get_dataset()
    ds = compute_relative_values(ds)
    t = Tree(ds)

    ds = t.dataset_for_clades_or_leaf(faust_clades)
    print 'Attributes: ', len(ds[0][2:])
    ds = median_discretization(ds)
    ds = discrete_dataset_cleaning(ds)
    print 'Attributes after cleaning: ', len(ds[0][2:])

    write_dataset_to_experiment('../experiments/5/Stool_maxent_discretized_faust_nodes_and_leafs', ds)

# run_discretization_faust_nodes_and_leafs()

def run_discretization_for_tree_depth(depth):
    """
    Experiment with all nodes at a given depth

    depth 6:
    Cleaning:
    Attributes:  166
    discrete dataset cleaning, removed bacteria:  65
    Attributes after cleaning:  101

    depth 5:
    Attributes:  66
    discrete dataset cleaning, removed bacteria:  21
    Attributes after cleaning:  45
    """
    from preprocessing import parser
    from preprocessing.tree import Tree
    from preprocessing import faust_parser
    import itemsets

    # Get the stool dataset and discretize it
    ds = parser.get_dataset()
    ds = compute_relative_values(ds)
    t = Tree(ds)
    ds = t.dataset_at_depth(depth)

    print 'Attributes: ', len(ds[0][2:])
    ds = median_discretization(ds)
    ds = discrete_dataset_cleaning(ds, 0.10)
    print 'Attributes after cleaning: ', len(ds[0][2:])

    write_dataset_to_experiment('../experiments/4/Stool_maxent_discretized_nodes_depth_6_010', ds)


# run_discretization_for_tree_depth(6)

def load_model():
    from mtv import MTV
    from utils import files
    import itemsets
    headers = files.parse_header_file('../experiments/1/Stool_maxent_discretized_all_nodes.headers')
    C = files.parse_dat_file('../experiments/1/summary.dat')
    D = files.parse_dat_file('../experiments/1/Stool_maxent_discretized_all_nodes.dat')
    # print clade names
    for X in C:
        print itemsets.to_index_list(X, headers)



def build_summary_table():
    # Construct summary table
    pass


def compare_to_faust():
    """
    Go through the faust results and check if our pairs occur
    :return:
    """
    from preprocessing.faust_parser import results
    from utils.files import parse_header_file, parse_dat_file
    from itemsets import to_index_list
    from itertools import combinations
    from preprocessing import parser
    from preprocessing.tree import Tree
    from preprocessing import faust_parser
    import itemsets

    # Get the stool dataset, and construct the
    # phylogenetic tree, so we can look up
    # ancestors of the patterns we have discovered
    ds = parser.get_dataset()
    tree = Tree(ds)

    genus_nodes = 0
    leaf_nodes = 0
    for node in tree.nodes.values():
        if node.depth == 6:
            genus_nodes += 1
        if node.is_leaf():
            leaf_nodes += 1
    print 'nodes: ', len(tree.nodes)
    print 'genus nodes: ', genus_nodes
    print 'leaf nodes: ', leaf_nodes


    headers = parse_header_file('../experiments/4/Stool_maxent_discretized_nodes_depth_6_010.headers')
    summary = parse_dat_file('../experiments/4/summary.dat')
    # Get the indeces of the pair
    bin_indeces = [to_index_list(x) for x in summary]

    # Get summary itemsets by their clade names
    patterns = []
    for indeces in bin_indeces:
        pattern = []
        for index in indeces:
            # co-occurrence
            if index > len(headers):
                index = index - len(headers)
            pattern.append(headers[index])
        patterns.append(pattern)

    faust_results = results()
    faust_results_copy = faust_results

    print 'Total faust results: ', len(faust_results)

    genues_in_faust = []
    for faust_result in faust_results:
        a = faust_result.clade_1
        b = faust_result.clade_2
        if a in tree.nodes and b in tree.nodes:
            node_a = tree.nodes[a]
            node_b = tree.nodes[b]
            if node_a.depth == 6 and node_b.depth == 6:
                genues_in_faust.append((a,b))
    print 'genus in faust: ', len(genues_in_faust)

    # Check every pair in every pattern for existence in the faust results
    print 'In faust:'
    for pattern in patterns:
        for a, b in combinations(pattern, 2):

            # Get nodes for our pattern
            node_a = tree.nodes[a]
            node_b = tree.nodes[b]
            # print 'depth a: ', node_a.depth
            # print 'depth b: ', node_b.depth

            for faust_result in faust_results:
                if faust_result.clade_1 == a or faust_result.clade_2 == a:
                    if faust_result.clade_1 == b or faust_result.clade_2 == b:
                        print 'Genus results:'
                        print 'Faust: %s %s, %d' % (faust_result.clade_1, faust_result.clade_2, faust_result.direction)
                        print 'MTV: %s, %s' % (a.replace('|', '-'), b.replace('|', '-'))

                        if faust_result in faust_results_copy:
                            faust_results_copy.remove(faust_result)
                        if (faust_result.clade_1, faust_result.clade_2) in genues_in_faust:
                            genues_in_faust.remove((faust_result.clade_1, faust_result.clade_2))
                elif faust_result.clade_1 in tree.nodes and faust_result.clade_2 in tree.nodes:
                    faust_a = tree.nodes[faust_result.clade_1]
                    faust_b = tree.nodes[faust_result.clade_2]

                    # Do not look for stuff that is too general
                    # only parent nodes
                    parent_depth = 3
                    if faust_a.depth <= parent_depth or faust_b.depth <= parent_depth:
                        continue

                    # Node a
                    is_descendents = 0
                    res = ''
                    if node_a.is_in_lineage(faust_a)  and node_b.is_in_lineage(faust_b):
                        res += 'a in a, b in b'
                        is_descendents = True
                    if node_a.is_in_lineage(faust_b)  and node_b.is_in_lineage(faust_a):
                        res += 'a in b, b in a'
                        is_descendents = True

                    if is_descendents:
                        if faust_result in faust_results_copy:
                            faust_results_copy.remove(faust_result)
                        print '####'
                        print 'Result in same lineage: ' + res
                        print 'Faust %d' % faust_result.direction
                        print 'Faust: %s, %s' % (faust_a.name, faust_b.name)
                        print 'MTV %s, %s' % (a, b)
                        print '####'


    print 'Faust results not found: ', [(fr.clade_1, fr.clade_2) for fr in faust_results_copy]
    print 'Not found at genus level: ', genues_in_faust

# compare_to_faust()



def clade_table():
    from utils.files import parse_header_file, parse_dat_file
    from itemsets import to_index_list
    headers = parse_header_file('../experiments/4/Stool_maxent_discretized_nodes_depth_6_005.headers')

    summary = parse_dat_file('../experiments/2b/summary.dat')
    # Get the indeces of the pair
    bin_indeces = [to_index_list(x) for x in summary]


    clades = []
    for i, indeces in enumerate(bin_indeces):
        pattern = []
        for index in indeces:
            # co-occurrence
            if index > len(headers):
                index = index - len(headers)

            pattern.append(headers[index])

        clades.append(pattern)
        # print '[' + ', '.join(pattern) + ']'
        print '%', i
        print '\item {\small '
        for index, clade in enumerate(pattern):
            line = '%s, ' % clade.replace('|', '-')
            if clade == pattern[-1]:
                line = line.replace(',', '')
            line = line.replace('_', '\\_')
            print '\t' + line
        print '}'
        print ''
    print 'Clades in itemset: '

    # for clade1, clade2 in clades:
    #     print '%s & %s\\\\' % (clade1, clade2)

# clade_table()


def clade_pair_abundances():
    from utils.files import parse_header_file, parse_dat_file
    from itemsets import to_index_list
    from plots.clade_correlation import plot_clades_relationships
    from itertools import combinations
    headers = parse_header_file('../experiments/1/Stool_maxent_discretized_all_nodes.headers')
    headers = headers + headers

    # parse_dat_file returns a list of ints
    summary = parse_dat_file('../experiments/tmp/summary.dat')
    # Get the indeces of the pair
    bin_indeces = [to_index_list(x) for x in summary]


    clades = []
    for pattern in bin_indeces:
        for binindex1, binindex2 in combinations(pattern, 2):
            clade1 = headers[binindex1]
            clade2 = headers[binindex2]
            clades.append((clade1, clade2))

    plot_clades_relationships(clades, '../experiments/2b/plots/')


# clade_pair_abundances()
def plot_clades_pairwise(clades):
    from plots.clade_correlation import plot_clades_relationships
    from itertools import combinations
    clade_pairs = []

    for clade1, clade2 in combinations(clades, 2):
        clade_pairs.append((clade1, clade2))


    plot_clades_relationships(clade_pairs, '../experiments/2a/plots/')

# clades = ['Bacteria|Bacteroidetes', 'Bacteroidia|Bacteroidales', 'Bacteroidetes|Bacteroidia', 'Bacteroidales|Bacteroidaceae', 'Bacteroidaceae|Bacteroides', 'Tenericutes|Mollicutes']
# clades = ['Veillonellaceae|Phascolarctobacterium', 'Veillonellaceae|Dialister', 'Porphyromonadaceae|Parabacteroides', 'Bacteroidaceae|Bacteroides']
# clades = ['Porphyromonadaceae|Barnesiella', 'Porphyromonadaceae|unclassified']
# clades = ['Erysipelotrichaceae|Turicibacter', 'Veillonellaceae|Veillonella', 'Pasteurellaceae|Haemophilus', 'Pasteurellaceae|unclassified']
# plot_clades_pairwise(clades)

def plot_clades():
    from plots.clade_correlation import plot_clades_relationships

    # Plot all combs
    # lh = ['Ruminococcaceae|Oscillibacter', 'Ruminococcaceae|Anaerotruncus', 'Ruminococcaceae|unclassified', 'Ruminococcaceae|Ruminococcus', 'Ruminococcaceae|Subdoligranulum', 'Incertae_Sedis_XIII|Anaerovorax', 'Coriobacteriaceae|Collinsella']
    # clades = []
    # for c in lh:
    #     clades.append((c, 'Bacteroidaceae|Bacteroides'))


    clades = [('Ruminococcaceae|unclassified', 'Bacteroidaceae|Bacteroides'),
              ('Ruminococcaceae|unclassified', 'Bacteroidia|Bacteroidales'),
              ('Ruminococcaceae|unclassified', 'Bacteroidales|Prevotellaceae'),
              ('Ruminococcaceae|unclassified', 'Bacteroidales|Porphyromonadaceae'),
              ('Ruminococcaceae|unclassified', 'Bacteroidales|Rikenellaceae'),
              ('Ruminococcaceae|unclassified', 'Bacteroidales|unclassified'),
              ('Ruminococcaceae|unclassified', 'Ruminococcaceae|Anaerofilum'),
              ('Ruminococcaceae|unclassified', 'Rikenellaceae|Alistipes')]
    # clades = [('Veillonellaceae|Phascolarctobacterium', 'Veillonellaceae|Dialister')]
    # clades = [['Ruminococcaceae|unclassified', 'Bacteroidaceae|Bacteroides']]
    # clades = [('Bacteroidaceae|Bacteroides', 'Prevotellaceae|unclassified')]
    # clades = [('Alcaligenaceae|Sutterella', 'Alcaligenaceae|Parasutterella')]
    clades = [('Veillonellaceae|Phascolarctobacterium', 'Veillonellaceae|Dialister')]
    clades= [('Alcaligenaceae|Sutterella', 'Alcaligenaceae|Parasutterella')]

    plot_clades_relationships(clades, '../experiments/2b/plots/')

# plot_clades()



# plot_clades()

# from preprocessing import faust_parser
# res = faust_parser.results()
# for r in res:
#     print r.clade_1 + ' + ' + r.clade_1



def format_stats(summary_file):
    """
    Reads a file containing a copy of the MTV output
    :param f:
    """

    heuristics, BIC_scores, independent_models, size_of_c, iteration_time, relation_ships, queries = read_run_results(summary_file)

    print 'heuristics =', heuristics
    print 'BIC_scores =', BIC_scores
    print 'independent_models =', independent_models
    print 'size_of_c =', size_of_c
    print 'iteration_time =', iteration_time
    print 'relation_ships: ', relation_ships
    print
    # for latex table
    print 'LATEX:'
    for index, i in enumerate(heuristics):
        line = '$ %d $ & ' % (index + 1)
        line += '$ %s $ & ' % relation_ships[index]
        line += '$ %.3f $ & ' % queries[index]
        line += '$ %d $ & ' % independent_models[index]
        line += '$ %d $ & ' % size_of_c[index]
        line += '$ %.1f $ &' % iteration_time[index]
        line += '$ %.3f $ & ' % heuristics[index]
        line += '$ %d $ ' % BIC_scores[index]
        print line + '\\\\'

# format_stats('../experiments/4/run_result.txt')

def write_tree():
    from preprocessing.tree import Tree
    from preprocessing.parser import get_dataset
    ds = get_dataset()
    t = Tree(ds)
    xml = t.root.to_xml()
    with open('../../../dev/thesis/data/tree.xml', 'wb') as fd:
       fd.write(xml)

write_tree()

# from preprocessing.tree import Tree
#
# from preprocessing import parser
#
# ds = parser.get_dataset()
#
# t = Tree(ds)
# clades = ['Bacteroidaceae|Bacteroides', 'Bacteroidia|Bacteroidales',
#           'Bacteroidales|Prevotellaceae', 'Bacteroidales|Porphyromonadaceae',
#           'Bacteroidales|Rikenellaceae', 'Bacteroidales|unclassified',
#             'Ruminococcaceae|unclassified', 'Ruminococcaceae|Anaerofilum']
#
# s = t.dot_graph_for_clades(clades)
# print s

def print_appendix_figures(figures=10):
    for i in range(figures):

        print '\\begin{figure}[h]'
        print '\\begin{center}'
        print '\\makebox[\\textwidth][c]{\\includegraphics[scale=0.3, width=1.2\\textwidth]{figures/experiment2a/%d.png}}%%' % i
        print '\\end{center}'
        print '\\caption{Phylogenetic tree for summary result %d in experiment 2a}' % i+1
        print '\\end{figure}\n'

# print_appendix_figures(10)

def print_report_clade_lots():
    from plots.clade_correlation import plot_clades_relationships
    clades= [('Ruminococcaceae|unclassified', 'Bacteroidaceae|Bacteroides'),
             ('Ruminococcaceae|unclassified','Bacteroidia|Bacteroidales'),
             ('Ruminococcaceae|unclassified', 'Rikenellaceae|Alistipes')]

    plot_clades_relationships(clades, '../experiments/4/plots/')

# print_report_clade_lots()

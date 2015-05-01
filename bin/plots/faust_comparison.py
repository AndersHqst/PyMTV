#!/usr/bin/env pypy

#
#   Scripts to compare a model with the Faust results.
#   faust_comparison will print a table formattet for a
#   latex table
#

# Careful. Gives warn('no cffi linalg functions and no _umath_linalg_capi module, expect problems.')
# running pypy
from preprocessing.parser import compute_relative_values
from preprocessing.parser import get_dataset

def load_model():
    from mtv import MTV
    from utils.files import parse_dat_file
    from utils.files import parse_header_file

    D = parse_dat_file('../experiments/4/Stool_maxent_discretized_all_nodes.dat')
    headers = parse_header_file('../experiments/1/Stool_maxent_discretized_all_nodes.headers')
    summary = parse_dat_file('../experiments/1/summary.dat')

    print 'Creating MTV object'
    mtv = MTV(D, summary, s=0.0)
    print 'Done'
    mtv.headers = headers
    return mtv


def faust_comparison(body_site='Stool'):
    from preprocessing.tree import Tree
    from preprocessing import faust_parser
    from scipy.stats import pearsonr, spearmanr
    from utils.correlation import phicoeff
    import itemsets
    from utils.rule_miner import association_rule

    mtv = load_model()

    # For each result in faust
    # translate clade names, to indeces
    #   create mtv query(headers) as a convenience
    #   let it throw a key error if header is unknown
    # crate table row with
    # key -/+ relationship faust pearson, spearman, association rules, in summary, agree? model phi

    faust_results = faust_parser.results(body_site)

    ds = get_dataset('Stool')
    ds = compute_relative_values(ds)

    # Construct a tree to get abundances for faust results
    tree = Tree(ds)
    count = 0
    # header = 'ID & relation ship & faust & pearson & spearman & phi & $X \\rightarrow Y$ & $Y \\rightarrow Y$ & in C\\\\'
    header = 'ID & faust & pearson & spearman & phi & $X \\rightarrow Y$ & $Y \\rightarrow Y$ & in C\\\\'
    print header
    for faust_result in faust_results:
        count += 1
        if abs(faust_result.direction) < 3:
            continue

        # make sure the faust result is in the tree
        # ex Clostridiales|IncertaeSedisXIV is not in the data set
        if not (tree.has_clade(faust_result.clade_1) and tree.has_clade(faust_result.clade_2)):
            continue

        clade_1 = faust_result.clade_1
        clade_2 = faust_result.clade_2

        # Get the nodes in the phylogenetic tree
        from_node = tree.node_for_clade_name(clade_1)
        to_node = tree.node_for_clade_name(clade_2)

        # Get the total abundance for the clades in the tree
        xs = tree.abundance_column_in_subtree(from_node)
        ys = tree.abundance_column_in_subtree(to_node)

        # Get faust correlations
        pearson = pearsonr(xs, ys)[0]
        spearman = spearmanr(xs, ys)[0]

        # Association rules
        # First, get the corresponding itemsets for the clades
        itemset_1 = itemsets.itemset_for_headers([clade_1], mtv.headers)
        itemset_2 = itemsets.itemset_for_headers([clade_2], mtv.headers)
        # Get association rules
        X_Y, Y_X = association_rule(mtv, itemset_1, itemset_2)

        # Get phi from model
        # that is, we need the counts for 00, 01, 10, 11. We will get these by the probabilities
        intersection = mtv.query_headers([clade_2, clade_1])
        X_estiamte = mtv.query_headers([clade_1])
        Y_estiamte = mtv.query_headers([clade_2])
        _00 = len(mtv.D) * (1 - X_estiamte - Y_estiamte + intersection)
        _01 = len(mtv.D) * (X_estiamte - intersection)
        _10 = len(mtv.D) * (Y_estiamte - intersection)
        _11 = len(mtv.D) * (intersection)
        phi = phicoeff(_00, _01, _10, _11)

        line = '$%d$ & ' % faust_result.id

        # line += '%s - %s & ' % (clade_1.replace('---', '|'), clade_2.replace('---', '|'))

        line += '$%d$ & ' % faust_result.direction

        line += '$%.3f$ & ' % pearson

        line += '$%.3f$ & ' % spearman

        line += '$%.3f$ & ' % phi

        line += '$%.3f,%.3f$ & ' % (X_Y.confidence, X_Y.lift)

        line += '$%.3f,%.3f$ & ' % (Y_X.confidence, Y_X.lift)

        in_C = 'N'
        if itemset_1|itemset_2 in mtv.C:
            in_C = 'Y'
        line +=in_C

        print line + '\\\\'

        # 67 faust results in total
        print 'count: ', count

# faust_comparison()
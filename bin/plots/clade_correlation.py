from matplotlib.pylab import plot, scatter, hist, ylabel, xlabel, savefig, close, title, figtext, grid
from preprocessing.parser import *
from preprocessing.discretization import *
from preprocessing.tree import Tree
from scipy.stats import pearsonr, spearmanr
from utils.correlation import phicoeff_lists
from utils.dataset_helpers import pairwise_remove_highest_values





def plot_clades_relationships(clade_pairs, folder, remove_highest=0):

    ds = get_dataset('Stool')
    ds = compute_relative_values(ds)

    # Use numeric Tree. We will discreteze relevant values ourselves
    # and plot the numeric correlations
    tree = Tree(ds)

    same_lineages = 0

    for clade_1, clade_2 in clade_pairs:

        # Get the nodes in the phylogenetic tree
        from_node = tree.node_for_clade_name(clade_1)
        to_node = tree.node_for_clade_name(clade_2)

        if from_node.is_in_lineage(to_node):
            same_lineages += 1

        title_text = 'Relative clade abundances'
        title(title_text)
        # find from-to bacteria abundances
        xs = []
        ys = []

        # Get the total abundance for hte clades in the tree
        abundance_from = tree.abundance_column_in_subtree(from_node)
        abundance_to = tree.abundance_column_in_subtree(to_node)

        x_zeroes = 0
        y_zeroes = 0


        # List the in-sample values for each sample
        for index, _row in enumerate(ds[1:]):
            from_abundance = abundance_from[index]
            to_abundance = abundance_to[index]

            if from_abundance == 0.0:
                x_zeroes += 1

            if to_abundance == 0.0:
                y_zeroes += 1

            # if from_abundance < 0.01 and to_abundance < 0.01:
            xs.append(from_abundance)
            ys.append(to_abundance)
            # else:
            #     xs.append(0.0)
            #     ys.append(0.0)

        if remove_highest > 0.0:
            xs, ys = pairwise_remove_highest_values(remove_highest, xs, ys)
            ys, xs = pairwise_remove_highest_values(remove_highest, ys, xs)

                        # All data was eventyally removed
            if len(xs) == 0:
                print 'Warning. lowest values:%f removed all data from caldes:%s-%s' % (remove_highest, from_node.name, to_node.name)
                continue

        print 'x zeroes: ', x_zeroes
        print 'y zeroes: ', y_zeroes
        grid(True)

        text_x = 0.67
        sample_points = 'Sample points: %d' % len(xs)
        figtext(text_x, 0.85, sample_points, fontsize=10)

        from_clade = from_node.name.replace('|', '-')
        to_clade = to_node.name.replace('|', '-')
        xlabel(from_clade, fontsize=10)
        ylabel(to_clade, fontsize=10)

        disc_x, discrete_xs = discretize_row(xs, maxent_discretization_splitter)
        disc_y, discrete_ys = discretize_row(ys, maxent_discretization_splitter)

        # plot discretization lines
        a, b = [disc_x, disc_x], [-0.001, max(ys)]
        c, d = [-0.001, max(xs)], [disc_y, disc_y]
        plot(a, b, c='r')
        plot(c, d, c='r')

        # Discrete bin sizes
        pairs = zip(discrete_ys, discrete_xs)
        _00 = '00: ' + str(len([x for x in pairs if x == (0,0)]))
        _01 = '01: ' + str(len([x for x in pairs if x == (0,1)]))
        _10 = '10: ' + str(len([x for x in pairs if x == (1,0)]))
        _11 = '11: ' + str(len([x for x in pairs if x == (1,1)]))
        figtext(text_x, 0.82, _00, fontsize=10)
        figtext(text_x, 0.79, _01, fontsize=10)
        figtext(text_x, 0.76, _10, fontsize=10)
        figtext(text_x, 0.73, _11, fontsize=10)
        phi_coeff, phi_r = pearsonr(discrete_xs, discrete_ys)
        phi = 'phi: %f, %f' % (phi_coeff, phi_r)
        figtext(text_x, 0.70, phi, fontsize=10)

        # Depth of nodes in the phylogenetic tree
        from_depth = 'x depth: %d' % from_node.depth
        to_depth = 'y depth: %d' % to_node.depth
        figtext(text_x, 0.67, from_depth, fontsize=10)
        figtext(text_x, 0.64, to_depth, fontsize=10)

        # Discretization values
        median_x = 'splitter x: %f' % disc_x
        median_y = 'splitter y: %f' % disc_y
        figtext(text_x, 0.61, median_x, fontsize=10)
        figtext(text_x, 0.58, median_y, fontsize=10)

        # Same lineage
        same_lineage = 'False'
        if tree.nodes_have_same_lineage(from_node, to_node):
            same_lineage = ' True'
        figtext(text_x, 0.55, 'Same lineage: ' + same_lineage, fontsize=10)

        # Pearson and spearman correlations
        try:
            pearson = pearsonr(xs, ys)
            spearman = spearmanr(xs, ys)

            pearson = 'Pearson: %.3f, %.3f' % (pearson[0], pearson[1])
            spearman = 'Spearman: %.3f, %.3f' % (spearman[0], spearman[1])

            figtext(text_x, 0.52, pearson, fontsize=10)
            figtext(text_x, 0.49, spearman, fontsize=10)
        except Exception, e:
            print e
            print 'clades1: ', from_node.name
            print 'clades2: ', to_node.name
            print 'xs: %s', xs
            print 'ys: %s', ys

        # vals = vals[:-20]
        scatter(xs, ys, s=1, color='#0066FF')
        file_name = folder + from_node.name.replace('|', '-').replace('/', '%') + '---' + to_node.name.replace('|', '-').replace('/', '%')
        # file_name = os.path.join(dir, file_name)
        # print '[RESULT] ', file_name
        # print vals
        savefig(file_name)
        close()

    print 'same lineages: ', same_lineages
    print 'node pairs: ', len(clade_pairs)
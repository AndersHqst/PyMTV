from utils.dataset_helpers import *
from matplotlib.pylab import plot, hist, ylabel, xlabel, savefig, close, title, figtext, grid
from preprocessing.discretization import *
from scipy.stats import tvar
import os
from preprocessing import parser
from preprocessing.tree import Tree
from preprocessing import faust_parser
import itemsets
from preprocessing.preprocessors import remove_empty_samples, compute_relative_values, discrete_dataset_cleaning
from preprocessing.discretization import *

def plot_bacteria_hist(folder, depth=6, mid_quantile=False):
    """
    Saves a histogram of abundances binned
    :return:
    """

    # Get the stool dataset and discretize it
    ds = parser.get_dataset()
    ds = compute_relative_values(ds)
    t = Tree(ds)
    ds = t.dataset_at_depth(depth)

    # Get header names to priint on the plots
    headers = ds[0][2:]

    for index, header  in enumerate(headers):

        node = t.node_for_clade_name(header)
        abundances = t.abundance_column_in_subtree(node)
        abundances = [round(x,3) for x in abundances]

        if mid_quantile:
            abundances.sort()
            abundances = abundances[int(len(abundances)*0.25): -int(len(abundances)*0.25)]

        xlabel('Relative abundance')
        ylabel('Bin size')

        title_text = header.replace('/','-').replace('|', '-')
        title(title_text)
        binwidth = 0.001
        bins, bin_sizes, patches = hist(abundances, bins=np.arange(min(abundances), max(abundances) + binwidth, binwidth), color='#0066FF')

        # Write discretized values
        threshold, discretized_abundances = discretize_row(abundances, maxent_discretization_splitter)
        _0 = '0: ' + str(len([x for x in discretized_abundances if x == 0]))
        _1 = '1: ' + str(len([x for x in discretized_abundances if x == 1]))

        text_x = 0.7

        smaples_text = 'Samples: %d' % len(abundances)
        figtext(text_x, 0.85, smaples_text, fontsize=10)

        threshold_text = 'Splitter: %f' % threshold
        figtext(text_x, 0.82, threshold_text, fontsize=10)
        figtext(text_x, 0.79, _0, fontsize=10)
        figtext(text_x, 0.76, _1, fontsize=10)

        # Draw threshold line
        max_bin = len(abundances)
        if len(bins) != 0:
            max_bin = max(bins)

        a, b = [threshold, threshold], [0, max_bin]
        plot(a, b, c='r')

        grid(True)

        # Write max and avg
        # max_abundance = 'max: %f' % max(abundances)
        # avg_abundance = 'avg: %f' % (sum(abundances) / float(len(abundances)))
        # figtext(text_x, 0.76, max_abundance, fontsize=10)
        # figtext(text_x, 0.73, avg_abundance, fontsize=10)

        # write variance
        # variance = 'var: %f' % tvar(abundances)
        # figtext(text_x, 0.70, variance, fontsize=10)

        # Save fig to folder
        if not (os.path.exists(folder)):
            os.makedirs(folder)
        file_name = os.path.join(folder, title_text)
        print 'Hist: ', file_name
        savefig(file_name)

        close()

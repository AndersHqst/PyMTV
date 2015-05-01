
from utils.dataset_helpers import *

def discretized_dataset(dataset, splitter_method):
    """
    Discretize row of a dataset with a given discretization method
    :param dataset: A dataset
    :param row_method: Discretization method for a rwo, ie abundance vector
    :return:
    """
    import numpy as np

    # Get the abundance matrix and discretize it
    abundances = abundance_matrix(dataset).T
    discrete_matrix = []
    for row in abundances:
        splitter, discrete_row = discretize_row(row, splitter_method)
        discrete_matrix.append(discrete_row)

    # transpose to dataset orientation
    discrete_matrix = np.array(discrete_matrix).T

    # Replace the abundance submatrix
    discretized_dataset = replace_abundance_matrix(dataset, discrete_matrix)

    return discretized_dataset


def discretize_row(row, splitter_method):
    discrete_row = []
    splitter = splitter_method(row)

    for val in row:
        if val <= splitter:
            discrete_row.append(0)
        else:
            discrete_row.append(1)

    return splitter, discrete_row


def median_discretization(dataset):
    from utils.stats import median
    return discretized_dataset(dataset, median)


def maxent_discretization_splitter_dist_threshold(row):
    from math import log

    maxent = 999999999999999999
    log_sums = []
    best_threshold = 0
    for val in row:

        threshold = val
        log_sum = 0
        for i in row:
            log_sum += log(abs(i-threshold)+0.1)

        log_sums.append(log_sum)

        if log_sum < maxent:
            maxent = log_sum
            best_threshold = val

    return best_threshold


def maxent_discretization_splitter_dist(dataset):
    """
    Discretization that maximizes the entropy of the sample point
    distances to the splitter
    :param dataset:
    :return:
    """

    return discretized_dataset(dataset, maxent_discretization_splitter_dist_threshold)


def maxent_discretization_splitter(row):
    """
    Return splitter in row that maximizes the entropy of the binary outcomes
    :param row:
    :return:
    """
    from scipy.stats import entropy
    from utils.stats import median

    highest_entropy = -1
    best_splitter = 0
    l = float(len(row))

    for splitter in row:

        # Get probabilities for the two outcomes w.r.t. the splitter
        a = len([x for x in row if x <= splitter]) / l
        b = len([x for x in row if x > splitter]) / l
        ent = entropy([a, b])

        if ent > highest_entropy:
            highest_entropy = ent
            best_splitter = splitter

    return best_splitter


def maxent_discretization(dataset):
    """
    Binary discretization, which gives the most uniform
    speration of the bins.
    :param dataset:
    :return:
    """
    return discretized_dataset(dataset, maxent_discretization_splitter)
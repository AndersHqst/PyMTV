import os

def pairwise_remove_highest_values(amount, values, pair_list, remove_zero_pairs=False):

    """
    Removes highest percentage of some values
    and removes to corresponding indeces in a paired list
    :param pct:
    :param values:
    :param pair_list:
    :return:
    """
    assert len(values) == len(pair_list)

    # The check for ignoring the lowest values when this is
    # larger than the input, is not strictly correct,
    # but is ok as we would like to not have empty
    # lists of abundance for plotting
    if len(values) == 0 or amount >= len(values):
        return values, pair_list


    pairs = zip(values, pair_list)
    if remove_zero_pairs:
        pairs = filter(lambda (x,y): x != 0.0 and y != 0.0, pairs)
    pairs.sort(lambda x,y: x[0] < y[0] and -1 or 1)

    list_a = []
    list_b = []

    for tup in pairs[:-int(amount)]:
        list_a.append(tup[0])
        list_b.append(tup[1])

    return list_a, list_b

def abundance_matrix(matrix):
    import numpy as np
    """ Return the submatrix of abundances"""
    # From row 1, from column 2
    try:
        return matrix[1:, 2:].astype(np.int)
    except:
        pass
    return matrix[1:, 2:].astype(np.float)


def replace_abundance_matrix(dataset, replacement):
    import numpy as np
    """
    Replace the abundance matrix of a dataset
    :param dataset: Dataset
    :param replacement: Replacement for abundance matrix
    :return:
    """
    clade_names = np.array(dataset[0][2:])
    ds = np.vstack((clade_names, replacement))

    # Attach sample columns, on left side
    left_columns = np.array(dataset)[:,0:2]
    ds = np.hstack((left_columns, ds))

    return ds


def dataset_with_negations(dataset, singletons):
        """
        Returns a dataset with negated attributes added.
        Every row in the dataset, will have bit-flipped version
        of itself, added as a left shift of the number of singletons

        Ex a dataset with binary 1001, 1110 will be returned as 01101001, 00011110
        :param dataset: List of integers representing binary transactions
        :return: Data set with negated attributes
        """
        dataset_with_negations = []
        for X in dataset:
            positive_attributes = len(singletons)
            mask = (2 ** positive_attributes) - 1
            negated_attributes = (mask ^ X) << positive_attributes
            row = X | negated_attributes
            dataset_with_negations.append(row)

        return dataset_with_negations


def is_negated_pattern(itemset, singletons):
    """
    :param itemset: Itemset
    :param singletons: All singleton, positive and negated
    :return: True if the itemset includes a negated attribute
    and tuple of (positive (values, negated_value), or False and (0.0)
    """

    # With cooccurrences, the data has been double, thus
    # dividing by 2 should always be the splitting value
    positive_bits = int(len(singletons) / 2.0)
    negated_attribute = itemset >> positive_bits
    positive_attributes = (2**positive_bits - 1) & itemset

    if negated_attribute != 0:
        return (True, (positive_attributes, negated_attribute<<positive_bits))

    return (False, (0,0))


def is_co_exclusion(output_folder, pattern_index):
    """
    Returns true if a patterns is a co-exclusion pattern
    due to a U<1

    This is not pretty, but currently the only way we can get this info
    from the output files from MTV. A better way would be
    to have a computer-readable file output from MTV with the
    U values. Here we strip it from the
    human-readable summary file
    :param output_folder: Output folder from MTV
    :param pattern_index: Pattern index, 0-indexed
    :return: True if co-exclusion pattern
    """

    prev_line = ''
    for line in open(os.path.join(output_folder, 'summary.txt'), 'rw'):
        if len(line) > 0 and len(prev_line) > 0:
            if line[0] == 'U' and int(prev_line[0]) == pattern_index:
                u = float(line.replace('U', '').replace('=', ''))
                return u < 1
        prev_line = line

    return False
import numpy as np

def compute_relative_values(dataset):
    """
    Returns a dataset with relative abundance counts, relative
    to the sample size.
    :param matrix:
    :return: Relative abundance matrix
    """

    # Keep headers
    relative_matrix = [list(dataset[0])]

    # relative per sample!
    # Skip row with headers
    for row in dataset[1:]:
        abundances = [int(x) for x in row[2:]]
        total = sum(abundances)
        new_row = list(row[:2])
        for value in abundances:
            val = value / float(total)
            if val > 1:
                pass
            new_row.append(value / float(total))
        relative_matrix.append(list(new_row))

    return np.array(relative_matrix)

def remove_empty_samples(dataset):
    # remove 0 samples
    no_zero_samples = []
    for index, row in enumerate(dataset):
        # Header rows
        if index < 2:
            no_zero_samples.append(row)
        else:
            sample_abundances = [int(x) for x in row[2:]]
            if max(sample_abundances) > 0:
                no_zero_samples.append(row)

    return no_zero_samples

def data_cleaning(dataset, threshold=2):
    """
    Returns a dataset where bacteria that to not occur more
    times than the threshold are removed
    :param dataset:
    :return:
    """

    no_zero_samples = remove_empty_samples(dataset)

    # Remove bacteria with abundance count <= 2
    # transpose the zero sample matrix to iterate
    # bacterias as rows
    cleaned_dataset = []
    for index, row in enumerate(np.array(no_zero_samples).T):
        # Header rows
        if index < 2:
            cleaned_dataset.append(row)
        else:
            abundances = [int(x) for x in row[1:]]
            if not max(abundances) <= threshold:
                cleaned_dataset.append(row)

    # Return the result, transposed to the original
    return np.array(cleaned_dataset).T


def discrete_dataset_cleaning(dataset, threshold=0.05):
    """
    Returns a dataset where bacteria that to not occur more
    times than the threshold are removed
    The default threshold of 0.05 would make good sense if this is also the
    support used in MTV.
    :param dataset: Dataset
    :param threshold: minumum fraction of values that should be 1
    :return:
    """

    cleaned_dataset = []
    removed_rows = 0

    for index, row in enumerate(np.array(dataset).T):
        # Header rows
        if index < 2:
            cleaned_dataset.append(row)
        else:
            # Discrete binary values
            abundances = [int(x) for x in row[1:]]
            if  threshold * len(abundances) <= sum(abundances):
                cleaned_dataset.append(row)
            else:
                removed_rows += 1

    print 'discrete dataset cleaning, removed bacteria: ', removed_rows
    # Return the result, transposed to the original
    return np.array(cleaned_dataset).T
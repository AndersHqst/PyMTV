def median(vals):
    """
    Return the median value from a list of numerical values
    :param vals: List of values
    :return:
    """
    from math import floor
    sorted_vals = sorted(vals)
    mid_index = int(floor(len(sorted_vals) / 2.0))
    if len(vals) % 2 == 1:
        return sorted_vals[mid_index]
    else:

        return (sorted_vals[mid_index] + sorted_vals[mid_index+1]) / 2.0


def fraction_splitter(vals, fraction):
    """
    Return the median value from a list of numerical values
    :param vals: List of values
    :return:
    """
    from math import ceil
    sorted_vals = sorted(vals)
    index = int(ceil(len(sorted_vals) * fraction))
    if len(vals) % 2 == 1:
        return sorted_vals[index]
    else:
        return (sorted_vals[index] + sorted_vals[index+1]) / 2.0


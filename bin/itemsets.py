"""
A set of helper functions for working with itemsets represented as integers
"""

from utils.memoisation import memoise
from utils.timer import *

def contains(a, b):
    """ True if a contains b """
    return a & b == b


def singletons(D):
    """
    Finds all singletons present in a dataset D and adds them to a set I.
    All rows in D are or'ed together to create a bit mask of all attributes.
    These are then all read as the 1 bit in the mask being shifted right
    :param I: Set singletons are added to
    :param D: Dataset where singletons will be taken from
    :return:
    """

    # Mask all bits observed in D
    mask = 0
    for d in D:
        mask = mask | d

    # Create set of singletons for all 1 bits
    return singletons_of_itemsets([mask])

def to_index_list(itemset, headers=None):
    """
    Returns a sorted list of the binary indeces of 1's
    found in the itemset
    :param itemset:
    :param headers: if headers are provided,
    indexes will be converted to header names in the returned list
    :return: List of binary indeces i itemset
    """
    l = []
    pos = 0
    while itemset != 0:
        if itemset & 1 == 1:
            l.append(pos)
        itemset = itemset >> 1
        pos += 1

    if not (headers is None):
        attribute_names = []
        for i in l:
            if i < len(headers):
                attribute_names.append(headers[i])
            else:
                attribute_names.append(i)
        return attribute_names

    return l


def to_index_lists(itemsets):
    """
    Use to_index_list on every itemset in the passed in
    itemsets and returns them in a list
    """
    l = []
    for itemset in itemsets:
        l.append(to_index_list(itemset))
    return l


def union_of_itemsets(itemsets):
    """Union of items in itemsets"""
    timer_start('union_of_itemsets')
    result = 0
    for c in itemsets:
        result = c | result
    timer_stop('union_of_itemsets')
    return result


def itemset_for_headers(itemset_headers, headers):

    itemset = 0

    for itemset_header in itemset_headers:
        itemset = itemset | 2**headers.index(itemset_header)

    return itemset


@memoise
def singletons_of_itemset(itemset):
    """
    Returns a set of singleton values in the itemset.
    This is a small modification of the bitCount algorithm
    ex https://wiki.python.org/moin/BitManipulation
    :param itemset:
    :return:
    """
    timer_start('Singletons of itemsets')

    singletons = set()
    while itemset:

        # Save the current itemset
        prev = itemset

        # Mask all bits but the lowest that is set
        itemset &= itemset - 1

        # Get the bit that was removed
        removed = itemset ^ prev

        singletons.add(removed)

    timer_stop('Singletons of itemsets')
    return singletons


def singletons_of_itemsets(itemsets):
    singletons = set()
    for itemset in itemsets:
        for singleton in singletons_of_itemset(itemset):
            singletons.add(singleton)
    return singletons


def iterate_singletons_of_itemset(itemset):
    """
    Iterate the singletons in an itemset.
    This is a small modification of the bitCount algorithm
    ex https://wiki.python.org/moin/BitManipulation
    :param itemset:
    :return:
    """

    while itemset:

        # Save the current itemset
        prev = itemset

        # Mask all bits but the lowest that is set
        itemset &= itemset - 1

        # Return bit that was removed
        yield itemset ^ prev


def binary_vectors_to_ints(binary_matrix):
    values = []
    for row in binary_matrix:
        val = 0
        pos = 0
        for bin_val in row:
            # the cast to int is necessay if the bin_val is
            # a numpy int64 - the binary operand will fail for higher values
            bit = 2 ** pos * int(bin_val)
            val = val | bit
            pos += 1
        values.append(val)
    return values


def itemset_from_binary_indeces(binaries):
    """
    Given a list of binaru indeces, returns an itemset as
    its integer representation
    :param binaries:
    :return:
    """
    val = 0
    for b in binaries:
        val = val | 2 ** b
    return val




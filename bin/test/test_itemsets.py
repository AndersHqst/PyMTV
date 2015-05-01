from unittest import TestCase

import itemsets

class ItemsetTest(TestCase):

    def test_itemsetunion(self):
        assert itemsets.contains(0, 0)
        assert itemsets.contains(1, 0)
        assert itemsets.contains(1, 1)
        assert itemsets.contains(3, 2)
        assert itemsets.contains(0, 2) == False
        assert itemsets.contains(2, 1) == False


    def test_singletons_of_itemset(self):
        singletons = [2**0, 2**1, 2**2, 2**3, 2**42]
        itemset = 0
        for i in singletons:
            itemset |= i
        for i in itemsets.singletons_of_itemset(itemset):
            assert i in singletons

        assert len(itemsets.singletons_of_itemset(itemset)) == len(singletons)


    def test_iterate_singletons_of_itemset(self):
        singletons = [2**0, 2**1, 2**2, 2**3, 2**42]
        itemset = 0
        for i in singletons:
            itemset |= i
        for i in itemsets.iterate_singletons_of_itemset(itemset):
            assert i in singletons

        assert len(itemsets.singletons_of_itemset(itemset)) == len(singletons)


    def test_singletons_of_dataset(self):
        # Build data set from 0,1,2,many + + rando singletons
        from random import randint
        I = [2**0, 2**1, 2**2, 2**3, 2**42]
        for i in range(10):
            I.append(2**randint(0, 100))

        # Construct dataset
        dataset = []
        for i in range(100):
            transaction = 0
            for i in range(randint(0,5)):
                transaction |= I[randint(0, len(I)-1)]
            dataset.append(transaction)

        # assert all singletons found are known
        for i in itemsets.singletons(dataset):
            assert i in I


    def test_to_index_list(self):
        itemset_0 = 2**0
        itemset_1 = 2**1
        itemset_2 = 2**2
        itemset_3 = 2**1 | 2**2

        assert itemsets.to_index_list(itemset_0) == [0]
        assert itemsets.to_index_list(itemset_1) == [1]
        assert itemsets.to_index_list(itemset_2) == [2]
        assert itemsets.to_index_list(itemset_3) == [1, 2]


    def test_union_of_itemsets(self):
        itemset_a = 2**0 | 2**4
        itemset_b = 2**1 | 2**17

        l = [itemset_a, itemset_b]

        union = itemsets.union_of_itemsets(l)

        assert union == itemset_a | itemset_b


    def test_binary_vestore_to_ints(self):

        results = [0, 1, 2, 3, 17]

        # These are in big-endian,
        # we think about transactions as big-endians
        # but represent them as little-endian internally
        i_0 = [0, 0, 0, 0]
        i_1 = [1, 0, 0, 0]
        i_2 = [0, 1, 0, 0]
        i_3 = [1, 1, 0, 0]
        i_17 = [1, 0, 0, 0, 1]

        M = [i_0, i_1, i_2, i_3, i_17]

        for index, i in enumerate(itemsets.binary_vectors_to_ints(M)):
            assert i == results[index]


    def test_itemset_from_binary(self):

        i_1 = [0, 0, 0]
        i_2 = [1]
        i_4 = [2]
        i_6a = [1, 2]
        i_6b = [2, 1]

        assert itemsets.itemset_from_binary_indeces(i_1) == 1
        assert itemsets.itemset_from_binary_indeces(i_2) == 2
        assert itemsets.itemset_from_binary_indeces(i_4) == 4
        assert itemsets.itemset_from_binary_indeces(i_6a) == 6
        assert itemsets.itemset_from_binary_indeces(i_6b) == 6









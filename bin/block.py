from utils.charitems import to_chars

class Block(object):
    """Block for itemsets in C"""
    def __init__(self):
        super(Block, self).__init__()
        self.union_of_itemsets = 0
        self.itemsets = set()
        self.singletons = set()
        self.block_size = 0
        self.cummulative_block_size = 0
        self.block_weight = 0
        self.cummulative_block_weight = 0

    def __str__(self):
        return to_chars(self.union_of_itemsets) + ' blocksize: ' + str(self.block_size) + ' blockweight: ' + str(self.block_weight)

    def __key(self):
        return self.union_of_itemsets

    def __eq__(x, y):
        return x.__key() == y.__key()

    def __hash__(self):
        return hash(self.__key())

    def __lt__(self, other):
        """Less-than to make class sortable
            Defined has the partial order on blocks,
            sets(T_1, C) < sets(T_2, C)
        """
        return self.union_of_itemsets & other.union_of_itemsets == self.union_of_itemsets
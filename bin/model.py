from __future__ import division
from math import log
from itertools import combinations
import sys

from utils.charitems import to_chars
from block import Block
from utils.timer import *
from utils.counter import *
from settings import *
import itemsets
from heuristic import h

sys.setrecursionlimit(1500)


class Model(object):

    def __init__(self, mtv):
        super(Model, self).__init__()

        # Reference to MTV for run parameters
        self.mtv = mtv

        # Model parameters
        self.u0 = 0
        self.U = {}

        # summary
        self.C = []

        # Singletons
        self.I = set()

        self.union_of_C = 0

        # Block w.r.t C
        self.T_c = []

        # Cached queries in FindBestItemSet
        self.query_cache = {}

        # Heurestics from h() for X in C at the time X was added
        self.heurestics = {}

        # BIC scores for X in C at the time X was added
        self.BIC_scrores = {}

        # Total weight of all singletons.
        # We keep this as a property as we
        # do not always need to recompute it
        self.total_weight = 0


    def p(self, T, y):
        """
        Return a model estiamte for y in block T
        :param T: A block
        :param y: Itemset
        :return: Model estimate for y in T
        """
        res = 1.0

        for x in self.C:
            assert not (x in self.I), "Singletons are not in summary calling p()"
            if itemsets.contains(T.union_of_itemsets, x):
                res = res * self.U[x]

        return self.u0 * res * T.block_weight


    def independence_estimate(self, y):
        """
        Return the indendence estiamte for y
        :param y: Itemset
        :return:
        """

        timer_start('independence_estimate')
        independence_estimate = 1.0
        for i in itemsets.iterate_singletons_of_itemset(y):
            independence_estimate *= self.U[i] / (1 + self.U[i])

        counter_inc('Independence estimates')
        timer_stop('independence_estimate')
        return independence_estimate


    def closure(self, y):
        """
        Returns a closure of itemset in C that har subsets of y
        :param y: Itemset to create clojure for
        :return: Clojure over C w.r.t y
        """
        closure = set()

        for itemset in self.C:
            if itemset & y == itemset:
                closure.add(itemset)

        return closure


    def query(self, y, scaling=False):
        """
        Query the probability on an itemset y.
        :param y: Itemset
        :param scaling: When scaling, query will not use any cache,
        and will recompute the total weight of songletons for every query
        :return: Estimate of y
        """

        if y in self.query_cache and not scaling:
            return self.query_cache[y]

        counter_inc('Total queries')

        estimate = 0.0

        # If y does not intersect any block,
        # we do an independence estimate
        if y & self.T_c[0].union_of_itemsets == 0:
            estimate = self.independence_estimate(y)
        else:
            counter_inc('Block queries')

            if scaling:
                self.compute_total_weight()

            T_c = self.compute_block_weights(y)

            timer_start('Compute p')
            for T in T_c:
                estimate += self.p(T, y)
            timer_stop('Compute p')

        if not scaling:
            self.query_cache[y] = estimate

        return estimate


    def compute_blocks(self):
        """
            Compute the set of blocks that C infer
            return: Topologically sorted blocks T_c
        """

        timer_start('Compute blocks')

        T_c = list()
        T_unions = set()

        # iterate the combination sizes in reverse
        a = range(len(self.C)+1)[::-1]
        for i in a:
            choose = i
            for comb in combinations(self.C, choose):
                union = itemsets.union_of_itemsets(comb)

                if not union in T_unions:
                    T_unions.add(union)
                    T = Block()
                    T.union_of_itemsets = union
                    T.singletons = itemsets.singletons_of_itemsets(comb)
                    T.itemsets = set(comb)

                    T_c.append(T)

        timer_stop('Compute blocks')
        return T_c


    def compute_block_sizes(self, T_c):
        for T in T_c:
            T.cummulative_block_size = 2 ** (len(self.I) - len(to_chars(T.union_of_itemsets)))
        for i, Ti in enumerate(T_c):
            Ti.block_size = Ti.cummulative_block_size
            for Tj in T_c[:i]:
                if Ti < Tj:
                    Ti.block_size = Ti.block_size - Tj.block_size
        return T_c


    def block_in_closure(self, T, closure):
        for c in closure:
            if not (c in T.itemsets):
                return False
        return True


    def compute_total_weight(self):
        self.total_weight = 1
        for i in self.I:
            self.total_weight *= (1 + self.U[i])


    def compute_block_weights(self, y=0):
        """
        Returns blocks for T_c + y, or all blocks if no y is passed
        :param y: Optional query parameter
        :return:
        """

        U = self.U
        blocks = []

        closure = self.closure(y)

        timer_start('Cummulative weight')
        for T in self.T_c:

            if self.block_in_closure(T, closure):

                blocks.append(T)

                T.cummulative_block_weight = self.total_weight

                # Remove singletons from y already covered by the block
                mask = y & T.union_of_itemsets
                ys = mask ^ y
                for i in itemsets.iterate_singletons_of_itemset(ys):
                    T.cummulative_block_weight *= U[i] / (1 + U[i])

                for i in T.singletons:
                    T.cummulative_block_weight *= U[i] / (1 + U[i])

        timer_stop('Cummulative weight')

        timer_start('Block weight')
        for i, Ti in enumerate(blocks):
            Ti.block_weight = Ti.cummulative_block_weight
            for Tj in blocks[:i]:
                if Ti < Tj:
                    Ti.block_weight = Ti.block_weight - Tj.block_weight

        timer_stop('Block weight')

        return blocks


    def iterative_scaling(self):

        # reset cache
        self.query_cache = {}

        # Initialize U and u0
        self.u0 = 2 ** -len(self.I)
        _C = self.I.union(self.C)
        for c in _C:
            self.U[c] = 1.0

        self.T_c = self.compute_blocks()

        timer_start('Iterative scaling')

        iterations = 0
        epsilon = 1e-4

        while iterations < 100:

            max_error = 0

            for x in _C:

                estimate = self.query(x, scaling=True)

                # if self.mtv.fr(x) == 0 or estimate == 0:
                #     msg = 'itemset %d has frequency=%f and p=%f. It should not be added to the summary' % (x, self.mtv.fr(x), estimate)
                #     assert False, msg

                fr_x = self.mtv.fr(x)
                if  abs(1 - fr_x) < float_precision:
                    fr_x = 0.9999999999
                elif fr_x < float_precision:
                    fr_x = 0.0000000001
                if  abs(1 - estimate) < float_precision:
                    estimate = 0.9999999999
                elif estimate < float_precision:
                    estimate = 0.0000000001

                self.U[x] = self.U[x] * (fr_x / estimate) * ((1 - estimate) / (1 - fr_x))
                self.u0 = self.u0 * (1 - fr_x) / (1 - estimate)

                max_error = max(max_error, 1 - min(fr_x, estimate) / max(fr_x, estimate))

            iterations += 1
            counter_max('Iterative scaling max iterations', iterations)

            if max_error < epsilon:
                break

        timer_stop('Iterative scaling')


    def score(self):
        # print 'score called'
        try:
            C = self.I.union(self.C)
            U = self.U
            u0 = self.u0
            D = self.mtv.D

            return -1 * (len(D) * (log(u0, 2) + sum([self.mtv.fr(x) * log(U[x], 2) for x in C])))

        except Exception, e:
            print 'Exception in score function, ', e
            # print 'Singletons: ', self.I
            for u in U:
                if U[u] <= 0:
                    'Negative U: ', u

            for x in C:
                if self.mtv.fr(x) <= 0:
                    'Itemset with 0 frequency: ', x

            print 'len of C: ', len(C)
            print 'len of D: ', len(D)
            print 'u0: ', u0

            exit()


    def add_to_summary(self, itemset):
        """
        Added an itemset to C and update current BIC score and heurestic
        for the itemset.
        :param itemset: An itemset to be added to C
        :return:
        """
        heuristic = h(self.mtv.fr(itemset), self.query(itemset))

        # Add X to summary
        self.C.append(itemset)
        self.union_of_C = itemsets.union_of_itemsets(self.C)
        self.heurestics[itemset] = heuristic

        # Update model
        self.iterative_scaling()

        # Compute score
        self.BIC_scrores[itemset] = self.score()


    def is_in_sumamry(self, y):
        """
        :param y: Itemset to look for
        :return: Ture if y is in C
        """
        for x in self.C:
            if y == x:
                return True
        return False


    def __str__(self):

        str = u'Summary: {0:s} '.format(self.C)
        str += u'U: {0:s} '.format(self.U)
        str += u'u0: {0:f} '.format(self.u0)

        return str
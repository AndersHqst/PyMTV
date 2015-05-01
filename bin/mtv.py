from math import log
from time import time

from settings import *
import itemsets
from model import Model
from heuristic import h
from utils.dataset_helpers import dataset_with_negations
from graph import Graph
from utils.timer import *
from utils.counter import *


class MTV(object):

    def __init__(self, D, initial_C=[], k=DEFAULT_K, m=DEFAULT_M, s=DEFAULT_S, z=DEFAULT_Z, v=DEFAULT_V, q=DEFAULT_Q, add_negated=DEFAULT_ADD_NEGATED, greedy=DEFAULT_GREEDY, headers=None):
        super(MTV, self).__init__()

        # Mine up to k itemsets
        self.k = k

        # Maximum itemset size
        self.m = m

        # Support
        self.s = s

        # Constraint on max model size
        self.q = q
        # If q is set, we will black list singletons from models
        # having reached the max size
        self.black_list_singletons = set()

        # Be verbose
        self.v = v

        # Header strings for attributes
        self.headers = headers

        # If set to True, MTV will also produce negated patterns
        self.add_negated = add_negated

        # If true FindBestItemset will be more greedy
        self.greedy = greedy

        # Number of candidate itemsets FindBestItemSet should search for
        # Will result in a list of top-z highest heuristics
        self.z = z

        # Dataset, is it right to always remove empty rows?
        tmp = []
        for i in D:
            if i != 0:
                tmp.append(i)
        self.D = tmp

        # Singletons
        self.I = itemsets.singletons(self.D)

        if self.add_negated:
            self.D  = dataset_with_negations(self.D, self.I)
            self.I = itemsets.singletons(self.D)


        # Cached frequency counts in D
        self.fr_cache = {}

        # Global summary
        self.C = list()

        self.union_of_C = itemsets.union_of_itemsets(self.C)

        self.BIC_scores = {}
        self.heuristics = {}

        # Create a model for holding all singletons
        # Singletons not used by model in the graph
        # will be in this model
        self.singleton_model = Model(self)
        self.singleton_model.I = self.I.copy()
        self.singleton_model.iterative_scaling()

        # Cached queries
        self.query_cache = {}

        # List to track history of disjoint components
        self.independent_components = []

        # List to track history of C size
        self.summary_sizes = []

        # List to track how large a search space was searched
        self.search_space = []

        # List to track history of timings of a loop in mtv
        self.loop_times = []

        # Initialize Graph of independent models
        self.graph = Graph()
        self.__init_graph(initial_C)

        self.BIC_scores['initial_score'] = self.score()


    def mtv(self):
        """
        Run the mtv algorithm
        """

        timer_stopwatch('run')

        # Added 0 loop_times for seeded itemsets
        for X in self.C:
            self.loop_times.append(0)

        # Run until we have converged
        while not self.finished():

            start = time()

            X = self.find_best_itemset()

            if not (self.validate_best_itemset(X)):
                break

            self.add_itemset(X)

            self.loop_times.append(time()-start)

            if self.v:
                print 'Found itemset (%.2f secs): %s, BIC-score: %f, model-sizes: %s, searched-nodes: %d' % (timer_stopwatch_time('run'), itemsets.to_index_list(X), self.BIC_scores[X], self.summary_sizes[-1], self.search_space[-1])


    def query(self, y):
        """
        Query using models intersected by y
        """

        timer_start('mtv_query')
        # query intersected models independently
        mask = y
        p = 1.0
        for model in self.graph.model_iterator():

            # Is this an intersected model?
            if y & model.union_of_C != 0:

                # get intersection
                intersection = model.union_of_C & mask

                # remove from mask
                mask = intersection ^ mask

                # query the intersected model
                p *= model.query(intersection)

        # disjoint singletons
        p *= self.singleton_model.query(mask)

        timer_stop('mtv_query')

        return p


    def query_headers(self, itemset_headers):
        """
        Query an itemset by its header names.
        This method will give a ValueError if the
        queried headers are not in the headers property
        of MTV
        :param itemset_headers: A list of headers
        :return: model query of the queried itemset
        """

        # itemset for the provided header names, will throw ValueError
        # if a header name is not in the self.headers property
        itemset = itemsets.itemset_for_headers(itemset_headers, self.headers)

        return self.query(itemset)


    def score(self):

        total_score = self.singleton_model.score()

        for model in self.graph.model_iterator():
            total_score += model.score()

        total_score += 0.5 * len(self.C) * log(len(self.D), 2)

        return total_score


    def finished(self):
        """
        Return True if the model has converged, or if k is provided, that k itemsets have been found.
        :return:
        """
        if not (self.k is None):
            return self.k <= len(self.C)

        if 1 < len(self.C):
            # If previous score is lower, the model score has increased
            # and we should finish.
            return self.BIC_scores[self.C[-2]] < self.BIC_scores[self.C[-1]]

        return False


    def add_itemset(self, X):
        """
        Add an itemset X to C and update MTV.
        warning: Adding itemsets to C should always be
         done with this methods, or MTV will be left in an
         invalid state.
        :param X: Itemset to be added to C
        :return:
        """
        heuristic = h(self.fr(X), self.query(X))

        # Add X to global summary
        self.C.append(X)
        self.union_of_C = itemsets.union_of_itemsets(self.C)
        self.heuristics[X] = heuristic

        self.update_graph(X)
        # Compute score
        self.BIC_scores[X] = self.score()


    def __init_graph(self, initial_c):
        """
        Init the graph with itemsets in C
        :return:
        """

        # Build graph
        for X in initial_c:
            self.add_itemset(X)

        # initialize independent models
        # and removed singletons from singleton model
        # for model in self.graph.independent_models():
        #     model.iterative_scaling()
        #     self.singleton_model.I -= model.I

        # finally initialize the singleton model
        # self.singleton_model.iterative_scaling()


    def update_graph(self, X):
        """
        Updates the graph with a new itemset X. This will always results in a new
        model being initialized. The new model's C corresponds to a new graph
        component, that may contain a merge of one or more existing
        graph components.
        """
        timer_start('Build independent models')

        new_model, components = self.graph.add_nodes(X, Model(self))
        new_model.iterative_scaling()

        self.update_model_constraints(new_model)

        # Update the singleton model
        self.singleton_model.I -= new_model.I
        self.singleton_model.iterative_scaling()

        timer_stop('Build independent models')

        self.graph_stats(components)


    def cached_itemset_query(self, X):
        """
        Helper function to cache queries.
        Note this must only be used from e.g. FindBestItemSet,
        where the model parameters are not altered between cache hits.
        :param X: Itemset
        :return: Model estimate of X
        """
        timer_start('Cached query')
        estimate = 0.0

        if X in self.query_cache:
            estimate = self.query_cache[X]
        else:
            estimate = self.query(X)
            self.query_cache[X] = estimate

        timer_stop('Cached query')
        return estimate


    def find_best_itemset(self):
        """
        Find best itemset in the sample
        space defined by I.
        Subject to model parameters z, m, s and
        the heuristic function h
        :return: Best itemset Z
        """

        # reset query caches
        self.query_cache = {}

        timer_start('Find best itemset')
        self.search_space.append(0)
        Z = self.find_best_itemset_rec(0, self.I.copy() - self.black_list_singletons, [(0,0)])
        timer_stop('Find best itemset')

        # Edge case, where we only find singletons not exactly described by the model
        # We search the top 10 Zs to see if there was a non singleton itemset
        for z in Z:
            if not (z in self.I) and z != 0:
                return z[0]
        print 'No valid z in Z: ', Z
        return Z[0][0]


    def validate_negated_pattern(self, X, y):
        """
        Return true if y unioned with X is a valid itemset
        when negated attributes are included. The builds on the
        assumption that only 1 negated attribute can be added to
        a pattern, which may not be meaningful in all domains.

        X|y will not be valid if a negated attribute is already in X
        or if the positive counterpart of y, is already in X

        :param X: Itemset
        :param y: Itemset to be unioned with X, if allowed when
        negated attributes are used
        :return: True, if y can be unioned with X
        """

        assert self.add_negated

        # MTV should be setup so half of the attributes
        # positive
        positive_attributes = int(len(self.I)/2.)

        # Check if y is a negated attribute
        if 2**positive_attributes <= y:

            # check no other negated attribute is set
            if X >> positive_attributes != 0:
                return False

            # Check if positive counterpart of y is set
            pos = y >> positive_attributes
            if pos & X == pos:
                return False

        else:
            # Check y negative counter part is not set
            pos = y << positive_attributes
            if pos & X == pos:
                return False

        return True


    def find_best_itemset_rec(self, X, Y, Z, X_length=0, parent_h=0):
        """
        :param X: itemset
        :param Y: remaining itemsets
        :param Z: currently best itemset
        :param s: min support
        :param m: max itemset size
        :param X_length: number of items in X. No pretty, but since X is an int,
                         this is the fastest way to know its length
        :return: Best itemsets Z
        """

        fr_X = self.fr(X)
        if fr_X < self.s or X in self.C:
            return Z

        p_X = self.cached_itemset_query(X)

        h_X = h(fr_X, p_X)

        # Greedy approach
        if self.greedy and h_X < parent_h:
            return Z

        if h_X > Z[-1][1] or len(Z) < self.z:
            # print 'Z improved: ', h_X
            Z.append((X, h_X))

            # Sort by descending  heuristic
            Z.sort(lambda x, y: x[1] < y[1] and 1 or -1)
            if self.z < len(Z):
                Z.pop()


        XY = X | itemsets.union_of_itemsets(Y)
        fr_XY = self.fr(XY)

        p_XY = self.cached_itemset_query(XY)

        b = max(h(fr_X, p_XY), h(fr_XY, p_X))

        if Z[0][0] == 0 or b > Z[-1][1]:

            if self.m == 0 or X_length < self.m:

                while 0 < len(Y):
                    y = Y.pop()

                    self.search_space[-1] += 1

                    # If we are mining negated patterns
                    # we have to check that y can be unioned with X
                    if not self.add_negated or self.validate_negated_pattern(X, y):
                        Z = self.find_best_itemset_rec(X | y, Y.copy(), Z, X_length + 1, parent_h=h_X)

        return Z


    def fr(self, x):
        """
        :param x: Itemset
        :return: Frequency of x in D
        """

        if x in self.fr_cache:
            return self.fr_cache[x]

        p = 0.0
        for xi in self.D:
            if itemsets.contains(xi, x):
                p += 1
        p = p / len(self.D)

        assert p <= 1.0

        self.fr_cache[x] = p

        return p


    def U(self):
        """
        Returns U over all models
        :return: U
        """
        U = {}

        for model in self.graph.model_iterator():
            U.update(model.U)

        return U


    def u0(self):
        """
        Returns u0 over all models
        :return: u0
        """
        u0 = self.singleton_model.u0

        for model in self.graph.model_iterator():
            u0 *= model.u0

        return u0


    def update_model_constraints(self, newest_model):
        """
        If C_i in some model has grow larger than q,
        we remove the singletons in the model from
        the search space
        :param newest_model:
        """
        if not (self.q is None):
            # Blacklist model singletons
            if len(newest_model.C) >= self.q:
                self.black_list_singletons = self.black_list_singletons.union(newest_model.I)


    def validate_best_itemset(self, itemset):
        """
        Returns True if an itemset can be added to C.
        Singletons or the empty set should not be added to C, but this can happen
        in cases where e.g. thresholds for support or minimum itemset size
        are too strict.
        """
        if itemset in self.I:
            print 'X was a singleton! These should not be possible from the heurestic'
            return False

        if itemset == 0:
            print 'Best itemset found was the empty set (0). This could ' \
                      'mean the heurestic could not find any itemset ' \
                      'not already predicted by the model, or above ' \
                      'the provided thresholds. Exiting MTV'
            return False

        return True


    def graph_stats(self, components):
        """
        Log stats on the graph
        :param components: Current Graph Components
        """
        self.independent_components.append(len(components))

        C_sizes = []
        for component in components:
            C_sizes.append(len(component.model.C))
        self.summary_sizes.append(C_sizes)

        counter_max('Independent models', len(components))

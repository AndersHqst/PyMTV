
from itertools import combinations
from itemsets import singletons_of_itemset, union_of_itemsets, to_index_list
from model import Model
import math
from settings import float_precision

class AssociationRule(object):
    def __init__(self):
        super(AssociationRule, self).__init__()

        # The rule is X -> Y, ie given an X, probability of Y
        self.X = 0
        self.Y = 0
        self.confidence = 0
        self.lift = 0

    def __str__(self):
        return '(conf:%f, lift:%f) %s -> %s' % (self.confidence, self.lift, to_index_list(self.X), to_index_list(self.Y))



def cached_query(mtv, itemset, cache):
    p = 0.0

    if itemset in cache:
        return cache[itemset]
    else:
        p = mtv.query(itemset)
        cache[itemset] = p

    return p

def association_rules(mtv, itemsets, use_observed_frequency=False):
    """
    Mines all association rules for a set of itemsets
    subject to the provided model and its dataset.
    Empty itemsets are ignored
    :param model: Model object
    :param itemsets: Itemsets to mine ruels for
     :param use_observed_frequency: If set to True, association rules will
     be based on the actual frequencies of the itemsets in D. Default is
     to use the estimate frequency of the model
    :return:
    """

    # Since we iterate several itemsets,
    # we may visit the same subsets X,Y several times
    # so we track visited nodes to avoid dublicates
    association_rules_set = set()

    cache = {}

    # List of all possible rules
    # that come from all subsets of the
    # itemsets
    rules = []

    for itemset in itemsets:

        singletons = singletons_of_itemset(itemset)

        for k in range(len(singletons)):
            choose_X = k + 1
            for comb in combinations(singletons, choose_X):

                X = union_of_itemsets(comb)

                prob_X = 0
                if use_observed_frequency:
                    prob_X = mtv.fr(X)
                else:
                    prob_X = cached_query(mtv, X, cache)

                Ys = set(singletons) - set(comb)

                for i in range(len(Ys)):
                    choose_Y = i + 1
                    for Y_comb in combinations(Ys, choose_Y):

                        Y = union_of_itemsets(Y_comb)

                        if not (X, Y) in association_rules_set:
                            association_rules_set.add((X, Y))
                            XY = X | Y
                            prob_XY = 0
                            prob_Y = 0
                            if use_observed_frequency:
                                prob_XY = mtv.fr(XY)
                                prob_Y = mtv.fr(Y)
                            else:
                                prob_XY = cached_query(mtv, XY, cache)
                                prob_Y = cached_query(mtv, Y, cache)

                            if prob_X > float_precision and prob_Y > float_precision:

                                conf = prob_XY / prob_X
                                lift = conf / prob_Y

                                rule = AssociationRule()
                                rule.X = X
                                rule.Y = Y
                                rule.confidence = conf
                                rule.lift = lift
                                rules.append(rule)


    # Return sorted list, regading the rules as
    # either association or disassociation rules
    # Association rules, descending prob, lift > 1 for true rules
    association_rules = filter(lambda rule: rule.lift > 1, rules)

    # Sort association rules descendingly, and cut away lower half
    association_rules.sort(lambda ar1, ar2: ar1.lift < ar2.lift and 1 or -1)
    split = int(math.ceil(len(association_rules)/2.))
    association_rules_high_lift = association_rules[:split]
    association_rules_low_lift = association_rules[split:]
    association_rules_high_lift.sort(lambda ar1, ar2: ar1.confidence < ar2.confidence and 1 or -1)
    association_rules_low_lift.sort(lambda ar1, ar2: ar1.confidence < ar2.confidence and 1 or -1)
    association_rules = association_rules_high_lift + association_rules_low_lift

    # Disassociation rules, ascending prob, lift < 1 for true rules
    disassociation_rules = filter(lambda rule: rule.lift < 1, rules)
    # Sort disassociation rules, ascindinglyly, and cut away lower half
    disassociation_rules.sort(lambda ar1, ar2: ar1.lift < ar2.lift and -1 or 1)
    d_split = int(math.ceil(len(disassociation_rules)/2.))
    disassociation_lower_lift = disassociation_rules[:d_split]
    disassociation_higher_lift = disassociation_rules[d_split:]
    disassociation_lower_lift.sort(lambda ar1, ar2: ar1.confidence < ar2.confidence and -1 or 1)
    disassociation_higher_lift.sort(lambda ar1, ar2: ar1.confidence < ar2.confidence and -1 or 1)
    disassociation_rules = disassociation_lower_lift + disassociation_higher_lift

    return association_rules, disassociation_rules


def association_rule(mtv, X, Y):
    """
    Returns association rules for X -> Y, Y -> X
    :param model: Mtv object to query probabilities
    :param X: An itemset
    :param Y: An itemset
    :return: Tuple of association rules (X->Y, Y->X)
    """

    prob_X = mtv.query(X)
    prob_Y = mtv.query(Y)
    prob_XY = mtv.query(X|Y)

    # X -> Y
    confX_Y = prob_XY / prob_X
    liftX_Y = confX_Y / prob_Y

    ruleX_Y = AssociationRule()
    ruleX_Y.X = X
    ruleX_Y.Y = Y
    ruleX_Y.confidence = confX_Y
    ruleX_Y.lift = liftX_Y

    # Y -> X
    confY_X = prob_XY / prob_Y
    liftY_X = confY_X / prob_X

    ruleY_X = AssociationRule()
    ruleY_X.X = Y
    ruleY_X.Y = X
    ruleY_X.confidence = confY_X
    ruleY_X.lift = liftY_X

    return (ruleX_Y, ruleY_X)


#
# Debug code
#
# D= [15, 4, 6, 0, 0, 2, 4, 6, 4, 2, 0, 2, 4, 4, 6, 2, 2, 4, 0, 7, 4, 4, 4, 5, 0, 10, 0, 2, 0, 0, 6, 6, 13, 0, 0, 4, 0, 0, 7, 6, 6, 6, 3, 0, 2, 0, 0, 2, 2, 15, 0, 4, 10, 0, 0, 12, 4, 6, 0, 6, 6, 2, 6, 6, 6, 4, 14, 0, 5, 4, 4, 4, 6, 2, 1, 0, 0, 0, 0, 0, 2, 7, 2, 6, 0, 2, 0, 4, 6, 0, 3, 0, 0, 2, 2, 0, 0, 12, 4, 4, 6, 14, 2, 4, 6, 4, 0, 2, 0, 0, 2, 4, 2, 2, 4, 4, 0, 4, 2, 2, 2, 4, 4, 2, 2, 0, 2, 1, 0, 2, 5, 7, 4, 0, 0, 6, 2, 2, 2, 0, 0, 4, 2, 10, 2, 0, 4, 0, 4, 4, 0, 6, 0, 0, 0, 6, 3, 2, 2, 0, 2, 6, 0, 0, 0, 6, 0, 4, 2, 6, 2, 2, 0, 3, 2, 4, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 4, 0, 4, 4, 4, 4, 2, 0, 0, 4, 1, 2, 6, 2, 4, 6, 4, 0, 0, 0, 4, 0, 0, 0, 2, 0, 0, 4, 4, 0, 0, 4, 4, 0, 0, 4, 0, 10, 0, 0, 4, 0, 0, 0, 0, 6, 4, 0, 0, 6, 4, 0, 6, 12, 0, 4, 0, 2, 0, 0, 4, 2, 6, 4, 4, 12, 4, 2, 2, 4, 4, 4, 0, 0, 2, 2, 2, 6, 4, 0, 0, 4, 4, 0, 0, 4, 8, 4, 6, 6, 0, 2, 4, 4, 0, 6, 6, 4, 0, 4, 6, 0, 0, 0, 2, 4, 0, 4, 0, 4, 0, 4, 6, 4, 0, 2, 2, 0, 4, 0, 0, 6, 8, 2, 2, 0, 0, 4, 6, 6, 2, 2, 0, 0, 2, 6, 0, 0, 2, 0, 0, 0, 0, 4, 4, 0, 9, 6, 4, 4, 6, 0, 0, 0, 6, 4, 2, 4, 2, 7, 0, 0, 6, 4, 6, 4, 4, 4, 0, 0, 0, 0, 4, 6, 6, 0, 8, 6, 4, 4, 15, 6, 6, 0, 4, 4, 4, 6, 4, 4, 2, 2, 0, 6, 6, 2, 11, 6, 0, 4, 6, 6, 4, 0, 6, 13, 0, 6, 6, 0, 10, 9, 6, 7, 0, 12, 6, 6, 6, 5, 0, 2, 2, 7, 0, 5, 0, 1, 2, 9, 2, 4, 0, 0, 6, 0, 0, 0, 0, 6, 6, 2, 2, 0, 0, 2, 2, 0, 4, 0, 0, 4, 4, 4, 0, 4, 0, 6, 4, 0, 2, 2, 2, 2, 0, 6, 2, 14, 4, 0, 0, 4, 2, 6, 0, 6, 2, 0, 6, 0, 4, 4, 0, 0, 12, 4, 6, 0, 6, 6, 6, 4, 6, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 0, 4, 6, 4, 0, 6, 4, 4, 2, 2, 0, 4, 0, 0, 0, 0, 4, 6, 4, 2, 0, 0, 6, 0, 4, 2, 0, 4, 0, 0, 4, 0, 6, 4, 4, 0, 6, 0, 4, 6, 6, 6, 4, 0, 2, 4, 0, 0, 6, 6, 4, 0, 4, 4, 0, 2, 6, 6, 0, 2, 6, 6, 6, 4, 6, 0, 6, 0, 6, 6, 6, 4, 6, 4, 6, 9, 6, 6, 2, 2, 4, 6, 6, 0, 0, 4, 6, 4, 0, 0, 4, 0, 4, 4, 0, 2, 2, 15, 4, 0, 4, 2, 0, 0, 4, 2, 2, 0, 0, 4]
# model = Model(D, s=0.0)
# model.mtv()
# association_rules, disassociation_rules = association_rules(mtv, model.C)
# print 'Association rules'
# for ar in association_rules:
#     print ar
#
# print '\nDisassociation rules'
# for ar in disassociation_rules:
#     print ar
#
# print 'done'
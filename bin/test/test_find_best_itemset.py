from unittest import TestCase

import test_helper
from utils import charitems


class TestFindBestItemset(TestCase):

    def test_find_best_itemset(self):

        model = test_helper.init_simple_model()

        best_itesmset = model.mtv.find_best_itemset()

        # We happen to know what the best itemsets should be in the
        # simple model..
        cd = charitems.to_binary('cd')
        assert cd == best_itesmset
        model.add_to_summary(cd)

        bc = charitems.to_binary('bc')
        best_itesmset = model.mtv.find_best_itemset()
        assert bc == best_itesmset
        model.add_to_summary(bc)

        ab = charitems.to_binary('ab')
        best_itesmset = model.mtv.find_best_itemset()
        assert ab == best_itesmset

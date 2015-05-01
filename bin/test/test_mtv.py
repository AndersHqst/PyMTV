from unittest import TestCase

import test_helper
from settings import float_precision
from utils import charitems


class TestMTV(TestCase):


    def helper_query_simple_singletons(self, mtv):
        """
        Assert query on all singletons are equal to their frequency
        :param model:
        :return:
        """
        # Assert precesion is within one transaction
        pct_one_transaction = (100/float(len(mtv.D))) / 100.
        for i in mtv.I:

            # Debug code to inspect singleton precesions
            # print 'q:%f, fr:%f' % (model.query(i), model.mtv.fr(i))

            assert abs(mtv.query(i) - mtv.fr(i)) < pct_one_transaction


    def test_mtv_query(self):
        # These values assume certain frequencies have
        # been set in the simple model

        mtv = test_helper.init_simple_model().mtv

        self.helper_query_simple_singletons(mtv)

        _0 = charitems.to_binary('')
        ab = charitems.to_binary('ab')
        mtv.add_itemset(ab)
        assert abs(mtv.query(ab) - 0.1) < float_precision
        assert abs(mtv.query(_0) - 1) < float_precision
        assert len(mtv.graph.components) == 1

        self.helper_query_simple_singletons(mtv)

        cd = charitems.to_binary('cd')
        mtv.add_itemset(cd)
        assert abs(mtv.query(cd) - 0.1) < float_precision
        assert abs(mtv.query(ab) - 0.1) < float_precision
        assert abs(mtv.query(_0) - 1) < float_precision
        assert len(mtv.graph.components) == 2

        self.helper_query_simple_singletons(mtv)

        bc = charitems.to_binary('bc')
        mtv.add_itemset(bc)
        assert abs(mtv.query(cd) - 0.1) < float_precision
        assert abs(mtv.query(ab) - 0.1) < float_precision
        assert abs(mtv.query(bc) - 0.1) < float_precision
        assert abs(mtv.query(_0) - 1) < float_precision
        assert len(mtv.graph.components) == 1

        self.helper_query_simple_singletons(mtv)





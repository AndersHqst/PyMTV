from unittest import TestCase

import test_helper
from settings import float_precision
from utils import charitems


class TestModel(TestCase):

    def test_independence_estimate(self):
        # These values assume certain frequencies have
        # been set in the simple model

        model = test_helper.init_simple_model()

        _0 = charitems.to_binary('')
        a = charitems.to_binary('a')
        b = charitems.to_binary('b')
        c = charitems.to_binary('c')
        d = charitems.to_binary('d')
        e = charitems.to_binary('e')
        assert abs(model.independence_estimate(_0) - 1) < float_precision
        assert abs(model.independence_estimate(a) - 0.5) < float_precision
        assert abs(model.independence_estimate(b) - 0.3) < float_precision
        assert abs(model.independence_estimate(c) - 0.2) < float_precision
        assert abs(model.independence_estimate(d) - 0.2) < float_precision
        assert abs(model.independence_estimate(e) - 0.1) < float_precision


    def test_frequency(self):
        # These values assume certain frequencies have
        # been set in the simple model

        model = test_helper.init_simple_model()

        _0 = charitems.to_binary('')
        a = charitems.to_binary('a')
        b = charitems.to_binary('b')
        c = charitems.to_binary('c')
        d = charitems.to_binary('d')
        e = charitems.to_binary('e')
        ab = charitems.to_binary('ab')
        cd = charitems.to_binary('cd')
        assert abs(model.mtv.fr(_0) - 1) < float_precision
        assert abs(model.mtv.fr(a) - 0.5) < float_precision
        assert abs(model.mtv.fr(b) - 0.3) < float_precision
        assert abs(model.mtv.fr(c) - 0.2) < float_precision
        assert abs(model.mtv.fr(d) - 0.2) < float_precision
        assert abs(model.mtv.fr(e) - 0.1) < float_precision
        assert abs(model.mtv.fr(ab) - 0.1) < float_precision
        assert abs(model.mtv.fr(cd) - 0.1) < float_precision


    def helper_query_simple_singletons(self, model):
        """
        Assert query on all singletons are equal to their frequency
        :param model:
        :return:
        """
        # Assert precesion is within one transaction
        pct_one_transaction = (100/float(len(model.mtv.D))) / 100.
        for i in model.I:

            # Debug code to inspect singleton precesions
            # print 'q:%f, fr:%f' % (model.query(i), model.mtv.fr(i))

            assert abs(model.query(i) - model.mtv.fr(i)) < pct_one_transaction


    def test_query(self):
        # These values assume certain frequencies have
        # been set in the simple model

        model = test_helper.init_simple_model()

        self.helper_query_simple_singletons(model)

        _0 = charitems.to_binary('')
        ab = charitems.to_binary('ab')
        model.add_to_summary(ab)
        assert abs(model.query(ab) - 0.1) < float_precision
        assert abs(model.query(_0) - 1) < float_precision

        self.helper_query_simple_singletons(model)

        cd = charitems.to_binary('cd')
        model.add_to_summary(cd)
        assert abs(model.query(cd) - 0.1) < float_precision
        assert abs(model.query(ab) - 0.1) < float_precision
        assert abs(model.query(_0) - 1) < float_precision

        self.helper_query_simple_singletons(model)

        bc = charitems.to_binary('bc')
        model.add_to_summary(bc)
        assert abs(model.query(cd) - 0.1) < float_precision
        assert abs(model.query(ab) - 0.1) < float_precision
        assert abs(model.query(bc) - 0.1) < float_precision
        assert abs(model.query(_0) - 1) < float_precision

        self.helper_query_simple_singletons(model)


    def test_closure(self):
        # These values assume certain frequencies have
        # been set in the simple model

        model = test_helper.init_simple_model()

        # Closure should be empty
        abc = charitems.to_binary('abc')
        empty_cls = model.closure(abc)
        assert len(empty_cls) == 0

        # ab should be in the closure
        ab = charitems.to_binary('ab')
        model.add_to_summary(ab)
        cls = model.closure(abc)
        assert ab in cls
        assert len(cls) == 1

        # Empty closure
        a = charitems.to_binary('a')
        cls = model.closure(a)
        assert len(cls) == 0

        # ab and cd in closure for abcd
        cd = charitems.to_binary('cd')
        model.add_to_summary(cd)
        abcd = charitems.to_binary('abcd')
        closure_two = model.closure(abcd)
        assert ab in closure_two
        assert cd in closure_two
        assert len(closure_two) == 2

        # nothing in closure
        bc = charitems.to_binary('bc')
        closure_empty2 = model.closure(bc)
        assert len(closure_empty2) == 0








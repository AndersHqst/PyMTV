from unittest import TestCase

import test_helper
from utils import charitems


class TestBlocks(TestCase):

    def test_compute_blocks(self):
        # print 'TEST BLOCKS'
        test_model = test_helper.init_simple_model()

        ab = charitems.to_binary('ab')
        cd = charitems.to_binary('cd')
        cde = charitems.to_binary('cde')

        test_model.add_to_summary(ab)

        assert len(test_model.T_c) == 2 # ab and the empty block
        assert ab in test_model.T_c[0].itemsets
        assert len(test_model.T_c[-1].itemsets) == 0

        test_model.add_to_summary(cd)
        assert len(test_model.T_c) == 4
        assert ab in test_model.T_c[0].itemsets
        assert cd in test_model.T_c[0].itemsets
        assert ab in test_model.T_c[1].itemsets or cd in test_model.T_c[1].itemsets
        assert ab in test_model.T_c[2].itemsets or cd in test_model.T_c[2].itemsets
        assert len(test_model.T_c[-1].itemsets) == 0

        # Debug code to inspect blocks
        # for i, block in enumerate(test_model.T_c):
        #     print 'block: ', i
        #     for itemset in block.itemsets:
        #         print charitems.to_chars(itemset)

        # Adding {cde} restricts that blocks with {cde, not cd} must be impossible
        test_model.add_to_summary(cde)
        assert len(test_model.T_c) == 6
        assert ab in test_model.T_c[0].itemsets
        assert cd in test_model.T_c[0].itemsets
        assert cde in test_model.T_c[0].itemsets
        for block in test_model.T_c:
            if cde in block.itemsets:
                assert cd in block.itemsets
        assert len(test_model.T_c[-1].itemsets) == 0







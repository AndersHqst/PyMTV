from unittest import TestCase

from utils import charitems


class TestCharItems(TestCase):

    def test_to_binary(self):
        assert charitems.to_binary('a') == 1
        assert charitems.to_binary('ab') == 3
        assert charitems.to_binary('b') == 2
        assert charitems.to_binary('') == 0
        assert charitems.to_binary('ad') == 9
        assert charitems.to_binary('da') == 9


    def test_to_chars(self):
        assert charitems.to_chars(0) == ''
        assert charitems.to_chars(1) == 'a'
        assert charitems.to_chars(2) == 'b'
        assert charitems.to_chars(3) == 'ab'
        assert charitems.to_chars(9) == 'ad'

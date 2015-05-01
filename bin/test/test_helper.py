import mtv
from utils import charitems


def init_simple_model():
    """

    IMPORTANT: do not modify the model constructed here
    as certain test methods make assumptions about the data.
    To create other init model helpers, copy this method..


    Initializes and returns a simple model
    with attributes {a,b,c,d,e}

    The simple model is constructed w.r.t. a MTV object
    that holds some synthetetic data D

    fr(a) = 0.5
    fr(b) = 0.3
    fr(c) = 0.2
    fr(d) = 0.2
    fr(e) = 0.1
    :return: Model
    """

    D = [
            charitems.to_binary('a'),
            charitems.to_binary('a'),
            charitems.to_binary('a'),
            charitems.to_binary('a'),
            charitems.to_binary('b'),
            charitems.to_binary('ab'),
            charitems.to_binary('bc'),
            charitems.to_binary('d'),
            charitems.to_binary('cd'),
            charitems.to_binary('e')
        ]

    m = mtv.MTV(D)

    return m.singleton_model

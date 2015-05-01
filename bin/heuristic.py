from settings import float_precision
from math import log

def h(x, y):
    """
    Heurestic for scoring an itemset
    :param x:
    :param y:
    :return:
    """
    # if x > y:
    #     return 0

    # x == 0 and y == 1
    if x < float_precision and abs(1.0 - y) < float_precision:
        return 0.0

    if y < float_precision:
        y = 0.000000000001
        # y =0.1

    # x == y
    if abs(x - y) < float_precision:
        return 0.0

    # y >= 1
    if abs(1.0 - y) < float_precision or y > 1.0:
        y = 0.9999999999
        # y=0.9

    # x == 0
    if x < float_precision:
        try:
            return -log(1.0 - y, 2)
        except Exception, e:
            print e
            print 'y: ', y
            print 'x: ', x
            exit()

    # x == 1.0
    if abs(1.0 - x) < float_precision:
        try:
            return -log(y, 2)
        except Exception, e:
            print e
            print 'y: ', y
            print 'x: ', x
            exit()


    return x * log(x / y, 2) + (1.0 - x) * log((1.0 - x) / ( 1.0 - y), 2)
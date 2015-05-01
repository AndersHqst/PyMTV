import math

def phicoeff_lists(x, y):
    """
    a, b, c, d are frequency counts for the various paired levels of dichotomous variables.
    2014 ahkj: Script found at http://adorio-research.org/wordpress/?p=10780. Modified
    to take two lists of binary values
        |     X
     Y  |  0     1
    ---------------
     0  |  a     b
     1  |  c     d

    :param x: List of binary values
    :param y: List of binary values
    :return:
    """

    xys = zip(x, y)

    a = len([0 for (x,y) in xys if x == 0 and y == 0])
    b = len([0 for (x,y) in xys if x == 1 and y == 0])
    c = len([0 for (x,y) in xys if x == 0 and y == 1])
    d = len([0 for (x,y) in xys if x == 1 and y == 1])

    ab = (a + b)
    if ab == 0:
       ab = 1

    cd = (c + d)
    if cd == 0:
       cd = 1

    ac = (a + c)
    if ac == 0:
       ac = 1

    bd = (b + d)
    if bd == 0:
       bd = 1


    return (a*d - b * c) / math.sqrt(ab * cd * ac * bd)

def phicoeff(a, b, c, d):
    """
    Script found at http://adorio-research.org/wordpress/?p=10780. Modified
    to take two lists of binary values
        |     X
     Y  |  0     1
    ---------------
     0  |  a     b
     1  |  c     d

    :param a:
    :param b:
    :param c:
    :param d:
    :return:
    """

    ab = (a + b)
    if ab == 0:
       ab = 1

    cd = (c + d)
    if cd == 0:
       cd = 1

    ac = (a + c)
    if ac == 0:
       ac = 1

    bd = (b + d)
    if bd == 0:
       bd = 1

    return (a*d - b * c) / math.sqrt(ab * cd * ac * bd)


def phi_correlation_in_model(mtv, X, Y):

    # that is, we need the counts for 00, 01, 10, 11. We will get these by the probabilities
    intersection = mtv.query(X|Y)
    X_estiamte = mtv.query(X)
    Y_estiamte = mtv.query(Y)
    _00 = int(len(mtv.D) * (1 - X_estiamte - Y_estiamte + intersection))
    _01 = int(len(mtv.D) * (X_estiamte - intersection))
    _10 = int(len(mtv.D) * (Y_estiamte - intersection))
    _11 = int(len(mtv.D) * (intersection))
    return phicoeff(_00, _01, _10, _11)

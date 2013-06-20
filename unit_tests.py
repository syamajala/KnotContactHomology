#!/usr/bin/sage -ipython
from sage.all import *
from braids import *

braid = Braid(2, [1, 1, 1])


def test_phi_b_ext():

    a = (1, 0)
    r = [[(2, 0), (2, 1), (1, 0)], (2, 1), [(1, 0), (1, 2), [(2, 0), (2, 1), (1, 0)]]]
    t = braid.phi_b_ext(*a)

    print "test 1:"
    "a1s"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "FAILED"

    print

    a = (2, 0)
    r = [(1, 0), (1, 2), [(2, 0), (2, 1), (1, 0)]]
    t = braid.phi_b_ext(*a)

    print "test 2:"
    "a2s"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "FAILED"

    print


def test_expand():

    a = [(0, 2), (0, 1), (1, 2)]
    r = [(-1, [(0, 2)]), (-1, [(0, 1), (1, 2)])]
    t = braid.expand(a)

    print "test 1:"
    print "-as2 - as1a12"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "FAILED"

    print

    a = [[(0, 2), (0, 1), (1, 2)], (2, 1), (1,  0)]
    r = [(1, [(0, 2)]), (1, [(0, 1), (1, 2)]), (-1, [(2, 1), (1, 0)])]
    t = braid.expand(a)

    print "test 2:"
    print "-(-as2 - as1*a12) - a21*a1s"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "FAILED"

    print 

    a = [(1, 0), (2, 1), [(0, 1), (1, 2), (2, 0)]]
    r = [(-1, [(1, 0)]), (1, [(2, 1), (0, 1)]), (1, [(2, 1), (1, 2), (2, 0)])]
    t = braid.expand(a)

    print "test 3:"
    print "-a1s - a21*(-as1 - a12*a2s)"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "FAILED"

    print

    a = [[(2, 0), (2, 1), (1, 0)], (2, 1), [(1, 0), (1, 2), [(2, 0), (2, 1), (1, 0)]]]
    r = [(1, [(2, 0)]), (1, [(2, 1), (1, 0)]), (1, [(2, 1), (1, 0)]), (-1, [(2, 1), (1, 2), (2, 0)]), (-1, [(2, 1), (1, 2), (2, 1), (1, 0)])]
    t = braid.expand(a)

    print "test 4:"
    print "-(-a2s - a21*a1s) - a21*(-a1s - a12*(-a2s - a21*a1s))"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "FAILED"

    print

    a = [(0, 1), [(0, 2), (0, 1), (1, 2)], (2, 1)]
    r = [(-1, [(0, 1)]), (1, [(0, 2), (2, 1)]), (1, [(0, 1), (1, 2), (2, 1)])]
    t = braid.expand(a)

    print "test 5:"
    print "-as1 - (-as2 - as1*a12)*a21"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "FAILED"

    print

    a = [(1, 0), [(1, 2), (2, 0), (1, 0)], [(2, 1), (1, 0), (2, 0)]]
    r = [(-1, [(1, 0)]), (-1, [(1, 2), (2, 1)]),
         (-1, [(1, 2), (1, 0), (2, 0)]), (-1, [(2, 0), (1, 0), (2, 1)]),
         (-1, [(2, 0), (1, 0), (1, 0), (2, 0)])]
    t = braid.expand(a)

    print "test 6:"
    print "-as1 - (-a12 + a2s*a1s)(-a21 - a1s*a2s)"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "FAILED"

    print

    a = [(1, 0), [(1, 2), [(0, 2), (0, 1), (1, 2)], (1, 0)],
         [(2, 1), (1, 0), (2, 0)]]
    r = [(-1, [(1, 0)]), (-1, [(1,  2), (2, 1)]), 
         (-1, [(1, 2), (1, 0), (2, 0)]), 
         (1, [(0, 2), (1, 0), (2, 1)]), (1, [(0, 2), (1, 0), (1, 0), (2, 0)]),
         (1, [(0, 1), (1, 2), (1, 0), (2, 1)]), 
         (1, [(0, 1), (1, 2), (1, 0), (1, 0), (2, 0)])]
    t = braid.expand(a)

    print "test 7:"
    print "-a1s - (-a12 - (-as2 - as1*a12)*a1s)*(-a12 - a1s*a2s)"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "FAILED"

    print

    print "test 9:"

    a = [(1, 0), [(1, 2), (2, 0), (1, 0)], [(2, 1), [(0, 2), (0, 1), (1, 2)], 
                                            (2, 0)]]
    r = [(-1, [(1, 0)]),
         (-1, [(1, 2), (2, 1)]),
         (1, [(1, 2), (0, 2), (2, 0)]),
         (1, [(1, 2), (0, 1), (1, 2), (2, 0)]),
         (-1, [(2, 0), (1, 0), (2, 1)]),
         (1, [(2, 0), (1, 0), (0, 2), (2, 0)]),
         (1, [(2, 0), (1, 0), (0, 1), (1, 2), (2, 0)])]
    t = braid.expand(a)

    print "-a1s - (-a12 - a2s*a1s)(-a21 - (-as2 - as1*a12)*a2s)"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "FAILED"

    print


test_phi_b_ext()
test_expand()


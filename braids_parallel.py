#!/usr/bin/sage -python
from collections import Counter
from itertools import count
from sage.all import *
from IPython.parallel import *
from braids import Braid, satellite


def setup_view():

    global c
    c = Client(profile='sage')
    global dview
    dview = c[:]
    dview.execute("import sys")
    dview.execute("import os")
    dview.execute("sys.path.append(os.path.expanduser('~/dev/KnotContactHomology'))")
    dview.execute("from braids_parallel import *")


def linearize(args):
    braid, i, j, l, constants = args
    s = l

    while not braid.linearized(s):
        r = 0

        for k, v in s._FreeAlgebraElement__monomial_coefficients.iteritems():
            e = 0
            els = k._element_list

            if len(els) > 1:
                g1 = braid.alg.gen(els[0][0])
                if els[0][1] != 1:
                    g2 = braid.alg.gen(els[0][0])
                    e = -2*g1 - 2*g2 - 4
                    e = e*(g2**(els[0][1]-2))
                    for g, p in els[1:]:
                        e = e*(braid.alg.gen(g)**p)
                else:
                    g2 = braid.alg.gen(els[1][0])
                    e = -2*g1 - 2*g2 - 4
                    if els[1][1] != 1:
                        e = e*(g2**(els[1][1]-1))
                    for g, p in els[2:]:
                        e = e*(braid.alg.gen(g)**p)

            elif els != [] and els[0][1] != 1:
                e = -4 - 4*braid.alg.gen(els[0][0])
                e = e*(braid.alg.gen(els[0][0])**(els[0][1]-2))
            elif els != []:
                e = braid.alg.gen(els[0][0])
            else:
                e = 1

            r = r + v*e

        s = r

    if not constants:

        l = 0
        one = braid.alg.monoid().one_element()

        for k, v in s._FreeAlgebraElement__monomial_coefficients.iteritems():
                if k != one:
                    for g, p in k._element_list:
                        l = l + v*(braid.alg.gen(g)**p)
        return (i, j, l)

    else:
        return (i, j, s)


def linearizem(self, m, constants=False):

    dview.push(dict(braid=self))
    l = CartesianProduct(range(0, self.s), range(0, self.s))
    r = Reference('braid')
    args = []

    for e in l:
        i, j = e
        args.append((r, i, j, m[i, j], constants))

    a = dview.map_sync(linearize, args)
    lm = matrix(self.alg, self.s, self.s)

    for e in a:
        i, j, r = e
        lm[i, j] = r

    return lm

Braid.linearizem = linearizem

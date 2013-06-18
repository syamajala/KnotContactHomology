#!/usr/bin/sage -python
from sage.all import *


class Braid:

    def __init__(self, strands, word):
        self.s = strands
        self.w = word
        self.alg = FreeAlgebra(ZZ, (self.s+1)*(self.s+1), self.aij(star=True))
        self.alg.inject_variables()

        m = []
        k = 0
        for i in range(0, self.s+1):
            n = []
            for j in range(0, self.s+1):
                n.append(self.alg.gen(k))
                k = k + 1
            m.append(n)

        self.astar = matrix(self.alg, self.s+1, self.s+1, m)
        self.a = matrix(self.alg, self.s, self.s, self.astar[1:self.s+1, 1:self.s+1])

    def aij(self, star=False):
        a = ''

        if not star:
            start = 1
        else:
            start = 0

        for i in range(start, self.s+1):
            for j in range(start, self.s+1):
                if star:
                    if i == 0:
                        a = a + 'as' + str(j)
                    elif j == 0:
                        a = a + 'a' + str(i) + 's'
                    else:
                        a = a + 'a'+str(i)+str(j)
                else:
                    a = a + 'a'+str(i)+str(j)
                a = a + ','

        return a[:len(a)-1]

    def phi_sk(self, z, x, y):
        """
        phi_sk returns phi_sigma_z(a_xy) 
        """

        k = z - 1
        i = x - 1
        j = y - 1

        if i == -1 and j == -1:
            return [(-1, -1)]

        if (i == k and j == k) or (i == k+1 and j == k+1):
            return [(-1, -1)]
        elif i == k + 1:
            if j != k and j != k + 1:
                return [(k, j)]
            else:
                return [(k, k+1)]
        elif j == k + 1:
            if i != k and i != k + 1:
                return [(i, k)]
            else:
                return [(k+1, k)]
        elif i == k and j != k and j != k + 1:
            return [(k+1, j), (k+1, k), (k, j)]
        elif j == k and i != k and i != k + 1:
            return [(i, k+1), (i, k), (k, k+1)]
        elif i != k and i != k+1 and j != k and j != k+1:
            return [(i, j)]

    def phi_ext_sk(self, z, x, y):
        """
        phi_ext_sk returns phi^ext_sigma_z(a_xy) 
        """

        k = z
        i = x
        j = y 

        if i == -1 and j == -1:
            return [(-1, -1)]

        if (i == k and j == k) or (i == k+1 and j == k+1) or (i == 0 and j == 0):
            return [(-1, -1)]
        elif i == k + 1:
            if j != k and j != k + 1:
                return [(k, j)]
            else:
                return [(k, k+1)]
        elif j == k + 1:
            if i != k and i != k + 1:
                return [(i, k)]
            else:
                return [(k+1, k)]
        elif i == k and j != k and j != k + 1:
            return [(k+1, j), (k+1, k), (k, j)]
        elif j == k and i != k and i != k + 1:
            return [(i, k+1), (i, k), (k, k+1)]
        elif i != k and i != k+1 and j != k and j != k+1:
            return [(i, j)]

    def phi_b_ext(self, i, j):
        """
        phi_b_ext returns phi_b_ext(a_ij)
        """
        a = [[]]

        a[0].append(self.phi_ext_sk(self.w[0], i, j))

        idx = 0
        for k in self.w[1:]:
            a.append([])
            for i in a[idx]:
                a[-1].append(self.phi_b_help(k, i, []))
            idx = idx + 1

        return a[-1][0]

    def phi_b_help(self, w, l, r):

        if l == [(-1, -1)]:
            return l
        elif l == []:
            return r
        else:
            if type(l[0]) == list:
                a = self.phi_b_help(w, l[0], [])
            else:
                a = self.phi_ext_sk(w, *l[0])
            if len(a) == 1:
                r.extend(a)
            else:
                r.append(a)
            return self.phi_b_help(w, l[1:], r)

    def phi_b_ext_a(self):
        """
        phi_b_ext_a returns phi_b_ext for all a_ij
        """

        a = []

        for i in range(0, self.s+1):
            for j in range(0, self.s+1):
                a.append(((i, j), self.phi_b_ext(i, j)))

        return dict(a)

    def expand(self, l, r, s):

        if type(l[0]) == list:
            r.extend(self.expand(l[0], [], -1*s))
        else:
            r.append((-1*s, [l[0]]))

        t = []
        if type(l[1]) == list:
            a = self.expand(l[1], [], -1*s)
            if type(l[2]) == list:
                pass
            else:
                for i in a:
                    x, y = i
                    y.append(l[2])
                    t.append((x, y))
                r.extend(t)
        else:
            if type(l[2]) == list:
                a = self.expand(l[2], [], s)
                for i in a:
                    x, y = i
                    if type(y) == list:
                        new_l1 = [l[1]]
                        new_l1.extend(y)
                        t.append((-1*s*x, new_l1))
                    else:
                        t.append((-1*s*x, [l[1], y]))
            else:
                t.append((-1*s, [l[1]]))
                t[-1][1].append(l[2])
            r.extend(t)

        return r

    def phi_ls(self, i):
        s = i-1
        m = identity_matrix(self.alg, self.s)
        m[s:s+2, s:s+2] = [[-1*self.a[s+1, s], -1], [1, 0]]
        return m

    def phi_rs(self, i):
        s = i-1
        m = identity_matrix(self.alg, self.s)
        m[s:s+2, s:s+2] = [[-1*self.a[s, s+1], 1], [-1, 0]]
        return m

braid = Braid(2, [1, 1, 1])


def test_expand():

    a = [(0, 2), (0, 1), (1, 2)]
    r = [(-1, [(0, 2)]), (-1, [(0, 1), (1, 2)])]
    t = braid.expand(a, [], 1)

    print "test 1:"
    print "-as2 - as1a12"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "test 1 failed"

    print

    a = [[(0, 2), (0, 1), (1, 2)], (2, 1), (1,  0)]
    r = [(1, [(0, 2)]), (1, [(0, 1), (1, 2)]), (-1, [(2, 1), (1, 0)])]
    t = braid.expand(a, [], 1)

    print "test 2:"
    print "-(-as2 - as1*a12) - a21*a1s"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "test 2 failed"

    print 

    a = [(1, 0), (2, 1), [(0, 1), (1, 2), (2, 0)]]
    r = [(-1, [(1, 0)]), (1, [(2, 1), (0, 1)]), (1, [(2, 1), (1, 2), (2, 0)])]
    t = braid.expand(a, [], 1)

    print "test 3:"
    print "-a1s - a21*(-as1 - a12*a2s)"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "test 3 failed"

    print

    a = [[(2, 0), (2, 1), (1, 0)], (2, 1), [(1, 0), (1, 2), [(2, 0), (2, 1), (1, 0)]]]
    r = [(1, [(2, 0)]), (1, [(2, 1), (1, 0)]), (1, [(2, 1), (1, 0)]), (-1, [(2, 1), (1, 2), (2, 0)]), (-1, [(2, 1), (1, 2), (2, 1), (1, 0)])]
    t = braid.expand(a, [], 1)

    print "test 4:"
    print "-(-a2s - a21*a1s) - a21*(-a1s - a12*(-a2s - a21*a1s))"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "test 4 failed"

    print

    a = [(0, 1), [(0, 2), (0, 1), (1, 2)], (2, 1)]
    r = [(-1, [(0, 1)]), (1, [(0, 2), (2, 1)]), (1, [(0, 1), (1, 2), (2, 1)])]
    t = braid.expand(a, [], 1)

    print "test 5:"
    print "-as1 - (-as2 - as1*a12)*a21"
    print "a = " + str(a)
    print "r = " + str(r)
    print "t = " + str(t)

    try:
        assert t == r
        print "PASS"
    except AssertionError:
        print "test 5 failed"

    print



braid = Braid(4, 13)
a = [[1, 0, 0, 0], [0, a14 + a13*a34, -1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
m1 = matrix(braid.alg, a)
b = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -a41, -1], [0, 0, 1, 0]]
m2 = matrix(braid.alg, b)
c = [[a31+a32*a21, -1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
m3 = matrix(braid.alg, c)
d = [[1, 0, 0, 0], [0, -a32, -1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
m4 = matrix(braid.alg, d)
i = m1*m2*m3*m4

a = [[-1*a13-a12*a23 + (-1*a14 - a12*a24)*(a43-a41*(-1*a13 - a12*a23) + a42*a23 + a42*(-1*a23-a21*(-a13 - a12*a23)) - (-1*a41 - a42*a21)*(-a13-a12*a23)), -a14-a12*a24, 1, 0],
     [-1*a32-(-1*a31 + a32*a21*a12) + (-1*a31-a32*a21)*a12, -a23-a21*(-a13-a12*a23), 0, 1],
     [1, 0, 0, 0],
     [0, 1, 0, 0]]
m1 = matrix(braid.alg, a)     

m1*i


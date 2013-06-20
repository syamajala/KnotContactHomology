#!/usr/bin/sage -python
from sage.all import *


class Braid:

    def __init__(self, strands, word):
        self.s = strands
        self.w = word
        g = self.gens(star=True)
        self.alg = FreeAlgebra(ZZ, len(g.split(',')), g)
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

        m = []
        for i in range(0, self.s):
            n = []
            for j in range(0, self.s):
                n.append(self.alg.gen(k))
                k = k + 1
            m.append(n)

        self.b = matrix(self.alg, self.s, self.s, m)

        m = []
        for i in range(0, self.s):
            n = []
            for j in range(0, self.s):
                n.append(self.alg.gen(k))
                k = k + 1
            m.append(n)
        
        self.c = matrix(self.alg, self.s, self.s, m)

    def gens(self, star=False):
        a = ''
        b = ''
        c = ''
        d = ''
        e = ''

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

                if i != 0 and j != 0:
                    b = b + 'b' + str(i) + str(j) + ','
                    c = c + 'c' + str(i) + str(j) + ','
                    d = d + 'd' + str(i) + str(j) + ','

        for i in range(1, self.s+1):
            e = e + 'e' + str(i) + str(i) + ','

        g = a + b + c + d + e[:-1]
        return g

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
        a = [self.phi_ext_sk(self.w[0], i, j)]

        for k in self.w[1:]:
            a.append(self.phi_b_help(k, a[-1], []))

        if len(a[-1]) == 1:
            return a[-1][0]
        else:
            return a[-1]

    def phi_b_help(self, w, l, r):
        """
        w is the braid generator, l is the term we want to call phi_ext_sk on
        r is the result
        """

        if l == []:
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

    def expand(self, l):
        
        return self.expand_help(l, [], 1)

    def expand_help(self, l, r, s):

        """
        Basically what we want to do is walk through a list like 
        l = [(0, 2), (0, 1), (1, 2)] and turn it into something like
        r = [(-1, [(0, 2)]), (-1, [(0, 1), (1, 2)])], so there are a
        few different cases. 

        if l[0] is a list then we should expand that, if not then it must be 
        (1, [(0, 2)]) by definition of phi_b_ext.

        if l[1] is a list, expand that, then if l[2] is a list as well expand
        it and foil, otherwise multiply each element of l[1] by l[2],

        similarly, if l[2] is a list, multiply each element of l[2] by l[1],
        otherwise just multiply l[1] and l[2]
        """

        if type(l[0]) == list:
            r.extend(self.expand_help(l[0], [], -1*s))
        else:
            r.append((-1*s, [l[0]]))

        t = []
        if type(l[1]) == list:
            a = self.expand_help(l[1], [], -1*s)
            if type(l[2]) == list:
                b = self.expand_help(l[2], [], s)
                t = []
                for i in a:
                    s_a, els_a = i
                    for j in b:
                        n_a = list(els_a)
                        s_b, els_b = j
                        n_a.extend(els_b)
                        t.append((s_a*s_b, n_a))
                r.extend(t)
            else:
                for i in a:
                    x, y = i
                    y.append(l[2])
                    t.append((x, y))
                r.extend(t)
        else:
            if type(l[2]) == list:
                a = self.expand_help(l[2], [], s)
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

    def phi_l_ij(self, i, j):

        a = self.phi_b_ext(i, 0)
        e = self.expand(a)

        r = []
        for i in e:
            s, ae = i
            if (j, 0) in ae:
                idx = ae.index((j, 0))
                r.append((s, ae[:idx] + ae[idx+1:]))

        return r

    def phi_lb(self):

        a = []

        for i in range(1, self.s+1):
            a.append([])
            for j in range(1, self.s+1):
                a[-1].append(self.astar_el(self.phi_l_ij(i, j)))

        return matrix(self.alg, a)

    def phi_r_ij(self, i, j):

        a = self.phi_b_ext(0, i)
        e = self.expand(a)

        r = []
        for i in e:
            s, ae = i
            if (0, j) in ae:
                idx = ae.index((0, j))
                r.append((s, ae[:idx] + ae[idx+1:]))

        return r
        
    def phi_rb(self):
        
        a = []

        for i in range(1, self.s+1):
            a.append([])
            for j in range(1, self.s+1):
                a[-1].append(self.astar_el(self.phi_r_ij(i, j)))

        return matrix(self.alg, a)

    def astar_el(self, l):
        a = 0

        for e in l:
            i, j = e

            if j == []:
                a = a + i
            else:
                t = 1
                for ee in j:
                    x, y = ee
                    t = t*self.astar[x][y]
                a = a + i*t

        return a

braid = Braid(2, [1, 1, 1])
braid.phi_lb()

#w = [2, 1, 3, 2]*3
#w.append(1)
#braid = Braid(4, w)
#print braid.phi_lb()

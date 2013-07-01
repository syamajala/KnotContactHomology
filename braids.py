#!/usr/bin/sage -python
from sage.all import *
from collections import deque


class Braid:

    def __init__(self, strands, word, diag=[]):
        self.s = strands
        self.w = word
        self.br = PolynomialRing(ZZ, 'l, li, nm, mi')
        self.br.inject_variables()
        g = self.gens(star=True)
        self.alg = FreeAlgebra(self.br, len(g.split(',')), g)
        self.alg.inject_variables()

        g = iter(range(len(self.alg.gens())))

        m = [[self.alg.gen(g.next()) for j in range(0, self.s+1)] 
             for i in range(0, self.s+1)]

        self.astar = matrix(self.alg, self.s+1, self.s+1, m)
        
        if diag:
            d = iter(diag)
            m = []
            for i in range(1, self.s+1):
                n = []
                for j in range(1, self.s+1):
                    if i == j:
                        n.append(d.next())
                    else:
                        n.append(self.astar[i, j])
                m.append(n)
            self.a = matrix(self.alg, self.s, self.s, m)
        else:
            self.a = matrix(self.alg, self.s, self.s, 
                            self.astar[1:self.s+1, 1:self.s+1])

        m = [[self.alg.gen(g.next()) for j in range(0, self.s)] 
             for i in range(0, self.s)]
        self.b = matrix(self.alg, self.s, self.s, m)

        m = [[self.alg.gen(g.next()) for j in range(0, self.s)]
             for i in range(0, self.s)]
        self.c = matrix(self.alg, self.s, self.s, m)

        m = [[self.alg.gen(g.next()) for j in range(0, self.s)]
             for i in range(0, self.s)]
        self.d = matrix(self.alg, self.s, self.s, m)

        m = identity_matrix(self.alg, self.s)
        for i in range(0, self.s):
            m[i, i] = self.alg.gen(g.next())

        self.e = matrix(self.alg, self.s, self.s, m)

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
                return [[1, [k, j]]]
            else:
                return [[1, [k, k+1]]]
        elif j == k + 1:
            if i != k and i != k + 1:
                return [[1, [i, k]]]
            else:
                return [[1, [k+1, k]]]
        elif i == k and j != k and j != k + 1:
            return [[-1, [k+1, j]], [-1, [k+1, k, k, j]]]
        elif j == k and i != k and i != k + 1:
            return [[-1, [i, k+1]], [-1, [i, k, k, k+1]]]
        elif i != k and i != k+1 and j != k and j != k+1:
            return [[1, [i, j]]]

    def phi_b_ext(self, i, j):
        """
        phi_b_ext returns phi_b_ext(a_ij)
        """
        a = deque(self.phi_ext_sk(self.w[0], i, j))
        head = deque([a.popleft()])

        for i in self.w[1:]:
            for n in range(len(head)):
                s, e = head.popleft()
                r = []
                for j in range(0, len(e), 2):
                    r.append(self.phi_ext_sk(i, e[j], e[j+1]))
                head.extend(self.expand(s, r))

        for i in self.w[1:]:
            for n in range(len(a)):
                s, e = a.popleft()
                r = []
                for j in range(0, len(e), 2):
                    r.append(self.phi_ext_sk(i, e[j], e[j+1]))
                a.extend(self.expand(s, r))

        return list(head) + list(a)

    def expand(self, os, l):
                
        a = []

        if len(l) == 1:
            for i in l[0]:
                s, e = i
                a.append([s*os, e])
            return a

        prev_len = 0
        for i in l:
            if len(i) == 1 and len(a) == 0:
                s, e = i[0]
                a.append([s*os, e])
            elif len(i) == 1 and len(a) != 0:
                if prev_len == 2:
                    s, e = i[0]
                    old_s1, old_e1 = a.pop()
                    old_s2, old_e2 = a.pop()
                    old_e1.extend(e)
                    old_e2.extend(e)
                    a.append([old_s2*s, old_e2])
                    a.append([old_s1*s, old_e1])
                else:
                    s, e = i[0]
                    old_s, old_e = a.pop()
                    old_e.extend(e)
                    a.append([old_s*s, old_e])
            elif len(i) == 2 and len(a) == 0:
                prev_len = 2
                for j in i:
                    s, e = j
                    a.append([s*os, e])
            elif len(i) == 2 and len(a) != 0:
                old_s, old_e = a.pop()
                for j in i:
                    new_e = list(old_e)
                    new_e.extend(j[1])
                    a.append([old_s*j[0], new_e])

        return a

    def simplify(self, l):
        a = []
        skip = []

        for i in l:
            c = l.count(i)
            if c > 1 and i not in skip:
                s, e = i
                a.append([s*c, e])
                skip.append(i)
            elif i in skip:
                pass
            else:
                a.append(i)

        return a

    def astar_el(self, l):
        a = 0
        for i in l:
            s, e = i
            t = 1
            for j in range(0, len(e), 2):
                t = t*self.astar[e[j], e[j+1]]
            a = a + s*t 

        return a

    def phi_l_b_ij(self, p, q, simplify=True):

        l = self.phi_b_ext(p, 0)
        if simplify:
            l = self.simplify(l)
        a = []

        for i in l:
            s, e = i
            for j in range(0, len(e), 2):
                if e[j] == q and e[j+1] == 0:
                    a.append([s, e[:j]+e[j+2:]])

        return a

    def phi_l_b(self):

        a = []
        
        for i in range(1, self.s+1):
            a.append([])
            for j in range(1, self.s+1):
                a[-1].append(self.astar_el(self.phi_l_b_ij(i, j, 
                                                           simplify=False)))

        return matrix(self.alg, a)

    def phi_r_b_ij(self, p, q, simplify=True):
        
        if p == q:
            l = self.phi_l_b_ij(p, q, simplify)
        else:
            l = self.phi_l_b_ij(q, p, simplify)

        for i in l:
            s, e = i
            e.reverse()

        return l

    def phi_r_b(self):

        a = []
        
        for i in range(1, self.s+1):
            a.append([])
            for j in range(1, self.s+1):
                a[-1].append(self.astar_el(self.phi_r_b_ij(i, j, 
                                                           simplify=False)))

        return matrix(self.alg, a)
    
    def diffb(self):

        i = identity_matrix(self.s)
        p = self.phi_l_b()

        return (i - p)*self.a

    def diffc(self):

        i = identity_matrix(self.s)
        p = self.phi_r_b()

        return self.a*(i - p)

    def diffd(self):

        i = identity_matrix(self.s)
        pl = self.phi_l_b()
        pr = self.phi_r_b()

        return (self.b*(i - pr)) - ((i - pl)*self.c)

    def diffe(self):
        
        return self.b + (self.phi_l_b()*self.c)

braid = Braid(2, [1, 1, 1], [-2, -2])
print
print braid.diffb()
print
print braid.diffc()
print
print braid.diffd()
print
print braid.diffe()

# w = [2, 1, 3, 2]*3
# w.append(1)
# braid2 = Braid(4, w)
# print braid2.phi_l_b()

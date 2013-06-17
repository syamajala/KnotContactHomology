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

    def phi_s(self, z, x, y):
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

    def phi_ext(self, z, x, y):
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

    def phi_b(self):
        a = [[]]

        for i in range(0, self.s+1):
            for j in range(0, self.s+1):
                a[0].append(self.phi_ext(self.w[0], i, j))

        idx = 0
        for k in self.w[1:]:
            a.append([])
            for i in a[idx]:
                a[-1].append(self.phi_b_help(k, i, []))
            idx = idx + 1
        
        return a

    def phi_b_help(self, w, l, r):

        if l == [(-1, -1)]:
            return l
        elif l == []:
            return r
        else:
            if type(l[0]) == list:
                a = self.phi_b_help(w, l[0], [])
            else:
                a = self.phi_ext(w, *l[0])
            if len(a) == 1:
                r.extend(a)
            else:
                r.append(a)
            return self.phi_b_help(w, l[1:], r)

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

    def phi_lb(self):
        m = self.phi_ls(self.w[0])

        for i in self.w[1:]:
            m = m*self.phi_ls(i)

        return m

    def phi_rb(self):
        m = self.phi_rs(self.w[0])
        
        for i in self.w[1:]:
            m = m*self.phi_rs(i)

        return m

#b = Braid(4, 13)
# print
# print b.a
# print

# for k in range(1, 4):
#     print "sigma_" + str(k)
#     for i in range(1, 5):
#         for j in range(1, 5):
#             print (i, j),
#             print "=",
#             print (b.phi_s(k, i, j))

braid = Braid(2, [1, 1, 1])
# print
# print braid.astar
# print

# for i in range(0, 3):
#     for j in range(0, 3):
#         print (i, j),
#         print "=",
#         print (braid.phi_ext(1, i, j))

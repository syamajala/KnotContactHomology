#!/usr/bin/sage -python
from sage.all import *


class TorusKnot:

    def __init__(self, crossings):
        self.n = crossings
        self.crossings = [(1, 2, self.n)] + \
                         [(i, i+1, i-1) for i in range(2, self.n)] + \
                         [(self.n, 1, self.n-1)]
        self.e1 = 1
        self.l = var('l')
        self.u = var('u')
        self.aij = matrix(SR, self.n, self.n, var(self.aij()))

    def aij(self):
        a = ''
        for i in range(1, self.n+1):
            for j in range(1, self.n+1):
                a = a + 'a'+str(i)+str(j)+','
        return a[:len(a)-1]

    def x(self, p, q):
        m = [[0]*self.n for __ in range(self.n)]
        m[p-1][q-1] = 1
        return matrix(m, Sparse=True)

    def l(self, a):
        return self.crossings[a][1]

    def r(self, a):
        return self.crossings[a][2]

    def o(self, a):
        return self.crossings[a][0]
        
    def psi_l1(self):
        m = self.l**(-self.e1)*self.x(1, self.r(1))

        for i in range(1, self.n):
            m = m + self.x(i+1, self.r(i))

        return m

    def psi_l2(self):
        m = matrix(SR, 3, 3)

        for i in range(self.n):
            m = m + u*self.x(i, self.l(i)) - \ 
            self.aij[self.l(i-1)][self.o(i-1)]*self.x(i, self.o(i))
            
        return m

t = TorusKnot(3)
print
print(t.psi_l2())
print
print(t.aij())

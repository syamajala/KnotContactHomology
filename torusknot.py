#!/usr/bin/sage -python
from sage.all import *


class TorusKnot:

    def __init__(self, crossings):
        self.n = crossings
        self.crossings = [(0, 1, self.n-1)] + \
                         [(i, i+1, i-1) for i in range(1, self.n-1)] + \
                         [(self.n-1, 0, self.n-2)]

        self.e1 = 1
        self.lmb = var('l')
        self.mu = var('u')
        self.aij = matrix(SR, self.n, self.n, var(self.aij()))

    def aij(self):
        a = ''
        for i in range(1, self.n+1):
            for j in range(1, self.n+1):
                a = a + 'a'+str(i)+str(j)+','
        return a[:len(a)-1]

    def x(self, p, q):
        m = matrix(self.n, self.n, Sparse=True)
        m[p,q] = 1
        return m

    def l(self, a):
        return self.crossings[a][1]

    def r(self, a):
        return self.crossings[a][2]

    def o(self, a):
        return self.crossings[a][0]

    def psi_l1(self):
        m = self.lmb**(-self.e1)*self.x(0, self.r(0))

        for i in range(1, self.n):
            m = m + self.x(i, self.r(i))

        return m

    def psi_l2(self):
        m = matrix(SR, self.n, self.n)

        for i in range(self.n):
            m = m + self.mu*self.x(i, self.l(i)) - \ 
            self.aij[self.l(i)][self.o(i)]*self.x(i, self.o(i))

        return m

    def psi_l(self):
        return self.psi_l1() + self.psi_l2()

    def psi_r1(self):
        m = self.lmb**(self.e1)*self.mu*self.x(self.r(0), 0)

        for i in range(1, self.n):
            m = m + self.mu*self.x(self.r(i), i)

        return m

    def psi_r2(self):
        m = matrix(SR, self.n, self.n)
         
        for i in range(self.n):
            m = m + self.x(self.l(i), i) - \
            self.aij[self.o(i)][self.l(i)]*self.x(self.o(i), i)

        return m

    def psi_r(self):
        return self.psi_r1() + self.psi_r2()

t = TorusKnot(3)

print
print t.psi_l()
print
print t.psi_r()
print
print t.aij()

t = TorusKnot(5)

print
print t.psi_l()
print
print t.psi_r()
print
print t.aij()


#!/usr/bin/sage -python
from sage.all import *


class Braid:

    def __init__(self, strands, crossings):
        self.s = crossings
        self.n = strands
        self.aij = matrix(SR, self.n, self.n, var(self.aij()))
    
    def aij(self):
        a = ''
        for i in range(1, self.n+1):
            for j in range(1, self.n+1):
                a = a + 'a'+str(i)+str(j)+','
        return a[:len(a)-1]

    def phi(self, z, x, y):
        k = z - 1
        i = x - 1
        j = y - 1

        if i == k + 1:
            if j != k and j != k + 1:
                return self.aij[k][j]
            else:
                return self.aij[k, k+1]
        elif j == k + 1:
            if i != k and i != k + 1:
                return self.aij[i][k]
            else:
                return self.aij[k+1][k]
        elif i == k and j != k and j != k + 1:
            return -1*self.aij[k+1][j]-self.aij[k+1][k]*self.aij[k][j]
        elif j == k and i != k and i != k + 1:
            return -1*self.aij[i][k+1]-self.aij[i][k]*self.aij[k][k+1]
        elif i != k and i != k+1 and j != k and j != k+1:
            return self.aij[i][j]
            

b = Braid(4, 13)
print
print b.aij()
print

for i in range(1, 5):
    for j in range(1, 5):
        print 
        print (i, j),
        print "=",
        print (b.phi(2, i, j))

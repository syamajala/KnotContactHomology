#!/usr/bin/sage -python
from sage.all import *


class Braid:

    def __init__(self, strands, crossings):
        self.s = strands
        self.c = crossings
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

        # for i in range(start, self.s+1):
        #     for j in range(start, self.s+1):
        #         a = a +'a'+str(i)+str(j)+','

        return a[:len(a)-1]

    def phi_s(self, z, x, y):
        k = z - 1
        i = x - 1
        j = y - 1

        if (i == k and j == k) or (i == k+1 and j == k+1):
            return 0
        elif i == k + 1:
            if j != k and j != k + 1:
                return self.a[k][j]
            else:
                return self.a[k, k+1]
        elif j == k + 1:
            if i != k and i != k + 1:
                return self.a[i][k]
            else:
                return self.a[k+1][k]
        elif i == k and j != k and j != k + 1:
            return -1*self.a[k+1][j]-self.a[k+1][k]*self.a[k][j]
        elif j == k and i != k and i != k + 1:
            return -1*self.a[i][k+1]-self.a[i][k]*self.a[k][k+1]
        elif i != k and i != k+1 and j != k and j != k+1:
            return self.a[i][j]

    def phi_ext(self, z, x, y):
        k = z
        i = x
        j = y 

        if (i == k and j == k) or (i == k+1 and j == k+1) or (i == 0 and j == 0):
            return 0
        elif i == k + 1:
            if j != k and j != k + 1:
                return self.astar[k][j]
            else:
                return self.astar[k, k+1]
        elif j == k + 1:
            if i != k and i != k + 1:
                return self.astar[i][k]
            else:
                return self.astar[k+1][k]
        elif i == k and j != k and j != k + 1:
            return -1*self.astar[k+1][j]-self.astar[k+1][k]*self.astar[k][j]
        elif j == k and i != k and i != k + 1:
            return -1*self.astar[i][k+1]-self.astar[i][k]*self.astar[k][k+1]
        elif i != k and i != k+1 and j != k and j != k+1:
            return self.astar[i][j]

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

braid = Braid(2, 3)
print
print braid.astar
print

for i in range(0, 3):
    for j in range(0, 3):
        print (i, j),
        print "=",
        print (braid.phi_ext(1, i, j))

#!/usr/bin/sage -python
from collections import Counter
from itertools import count
from sage.all import *
from sage.algebras.free_algebra_element import FreeAlgebraElement


class Braid:

    def __init__(self, strands, word, polyring=False, diag=0, linear=True):
        self.s = strands
        self.w = word
        self.polyring = polyring
        self.linear = linear
        g = self.gens(star=True)

        if self.polyring:
            r = PolynomialRing(ZZ, 'lmb, lmbinv, mu, muinv, U, Uinv')
            ideal = r.ideal(r.gen(0)*r.gen(1)-1, r.gen(2)*r.gen(3)-1, 
                            r.gen(4)*r.gen(5)-1)
            self.br = QuotientRing(r, ideal, 'lmb, lmbinv, mu, muinv, U, Uinv')
            self.br.inject_variables()
            self.alg = FreeAlgebra(self.br, len(g.split(',')), g)
        else:
            self.alg = FreeAlgebra(ZZ, len(g.split(',')), g)
        self.alg.inject_variables()

        g = iter(xrange(len(self.alg.gens())))

        m = [[self.alg.gen(g.next()) for j in xrange(0, self.s+1)] 
             for i in xrange(0, self.s+1)]

        self.astar = matrix(self.alg, self.s+1, self.s+1, m)

        m = [[self.alg.gen(g.next()) for j in xrange(0, self.s)] 
             for i in xrange(0, self.s)]
        self.b = matrix(self.alg, self.s, self.s, m)
        
        m = [[self.alg.gen(g.next()) for j in xrange(0, self.s)]
             for i in xrange(0, self.s)]
        self.c = matrix(self.alg, self.s, self.s, m)

        m = [[self.alg.gen(g.next()) for j in xrange(0, self.s)]
             for i in xrange(0, self.s)]
        self.d = matrix(self.alg, self.s, self.s, m)

        if diag and not self.polyring:
            m = []
            for i in xrange(1, self.s+1):
                n = []
                for j in xrange(1, self.s+1):
                    if i == j:
                        n.append(diag)
                    else:
                        n.append(self.astar[i, j])
                m.append(n)
            self.a = matrix(self.alg, self.s, self.s, m)

            m = identity_matrix(self.alg, self.s)
            for i in xrange(0, self.s):
                m[i, i] = self.alg.gen(g.next())

            self.e = matrix(self.alg, self.s, self.s, m)

        elif self.polyring:
            m = []
            for i in xrange(1, self.s+1):
                n = []
                for j in xrange(1, self.s+1):
                    if i < j:
                        n.append(-1*mu*self.astar[i, j])
                    elif i == j:
                        n.append(mu-1)
                    elif i > j:
                        n.append(self.astar[i, j])
                m.append(n)
            self.a = matrix(self.alg, self.s, self.s, m)
            
            m = []
            for i in xrange(1, self.s+1):
                n = []
                for j in xrange(1, self.s+1):
                    if i < j:
                        n.append(-1*mu*U*self.astar[i, j])
                    elif i == j:
                        n.append(mu*U-1)
                    elif i > j:
                        n.append(self.astar[i, j])
                m.append(n)

            self.ahat = matrix(self.alg, self.s, self.s, m)

            m = []
            for i in xrange(0, self.s):
                n = []
                for j in xrange(0, self.s):
                    if i > j:
                        n.append(self.b[i, j])
                    elif i == j:
                        n.append(0)
                    elif i < j:
                        n.append(-mu*self.b[i, j])
                m.append(n)

            self.b = matrix(self.alg, self.s, self.s, m)

            m = []
            for i in xrange(0, self.s):
                n = []
                for j in xrange(0, self.s):
                    if i > j:
                        n.append(self.b[i, j])
                    elif i == j:
                        n.append(0)
                    elif i < j:
                        n.append(U*self.b[i, j])
                m.append(n)

            self.bhat = matrix(self.alg, self.s, self.s, m)

            m = [[self.alg.gen(g.next()) for j in xrange(0, self.s)]
                 for i in xrange(0, self.s)]
            self.e = matrix(self.alg, self.s, self.s, m)

            m = [[self.alg.gen(g.next()) for j in xrange(0, self.s)]
                 for i in xrange(0, self.s)]
            self.f = matrix(self.alg, self.s, self.s, m)

            self.writhe = 0
            for k, v in Counter(self.w).iteritems():
                self.writhe = self.writhe + v

            self.cl = identity_matrix(self.alg, self.s)
            self.cl[0, 0] = lmb*((-1*muinv)**self.writhe)*(Uinv**((self.writhe-self.s+1)/2))
            self.cli = identity_matrix(self.alg, self.s)
            self.cli[0, 0] = lmbinv*((-1*mu)**self.writhe)*(U**((self.writhe-self.s+1)/2))

        else:
            self.a = matrix(self.alg, self.s, self.s, 
                            self.astar[1:self.s+1, 1:self.s+1])

        self.phi_l_bm = None
        self.phi_r_bm = None
        self.phi_l_b_linm = None
        self.phi_r_b_linm = None

        self.gdicta = []
        self.gdictb = []
        self.gdictc = []
        self.gdictd = []
        self.gdicte = []
        self.gdictf = []

        for e in list(self.alg.gens()):
            s = str(e)
            i = s[1]
            if i == 's':
                i = 0
            else:
                i = int(i)
            j = s[2]
            if j == 's':
                j = 0
            else:
                j = int(j)

            if s[0] == 'a':
                self.gdicta.append((e, (i, j)))
            elif s[0] == 'b':
                self.gdictb.append((e, (i, j)))
            elif s[0] == 'c':
                self.gdictc.append((e, (i, j)))
            elif s[0] == 'd':
                self.gdictd.append((e, (i, j)))
            elif s[0] == 'e':
                self.gdicte.append((e, (i, j)))
            elif s[0] == 'f':
                self.gdictf.append((e, (i, j)))
        
        self.gdicta = dict(self.gdicta)
        self.gdictb = dict(self.gdictb)
        self.gdictc = dict(self.gdictc)
        self.gdictd = dict(self.gdictd)
        self.gdicte = dict(self.gdicte)
        self.gdictf = dict(self.gdictf)
        
        self.i = identity_matrix(self.alg, self.s)

    def gens(self, star=False):
        a = ''
        b = ''
        c = ''
        d = ''
        e = ''
        f = ''

        if not star:
            start = 1
        else:
            start = 0

        for i in xrange(start, self.s+1):
            for j in xrange(start, self.s+1):
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

        g = a + b + c + d

        if self.polyring:
            for i in xrange(1, self.s+1):
                for j in xrange(1, self.s+1):
                    e = e + 'e' + str(i) + str(j) + ','
            for i in xrange(1, self.s+1):
                for j in xrange(1, self.s+1):
                    f = f + 'f' + str(i) + str(j) + ','
            g = g + e + f[:-1]
        else:
            for i in xrange(1, self.s+1):
                e = e + 'e' + str(i) + str(i) + ','
            g = g + e[:-1]

        return g

    def phi_ext_sk(self, k, a):
        """
        phi_ext_sk returns phi^ext_sigma_z(a_xy) 
        """

        i, j = self.gdicta[a]

        if k > 0:

            if i == k + 1:
                if j != k and j != k + 1:
                    return self.astar[k, j]
                else:
                    return self.astar[k, k+1]
            elif j == k + 1:
                if i != k and i != k + 1:
                    return self.astar[i, k]
                else:
                    return self.astar[k+1, k]
            elif i == k and j != k and j != k + 1:
                return -1*self.astar[k+1, j]-self.astar[k+1, k]* \
                    self.astar[k, j]
            elif j == k and i != k and i != k + 1:
                return -1*self.astar[i, k+1]-self.astar[i, k]* \
                    self.astar[k, k+1]
            elif i != k and i != k+1 and j != k and j != k+1:
                return self.astar[i, j]

        elif k < 0:

            k = -1*k

            if i == k + 1:
                if j == k:
                    return self.astar[k, k+1]
                else: 
                    return -1*self.astar[k, j]-self.astar[k, k+1]* \
                        self.astar[k+1, j]
            elif j == k + 1:
                if i == k:
                    return self.astar[k+1, k]
                else:
                    return -1*self.astar[i, k]-self.astar[i, k+1]* \
                        self.astar[k+1, k]
            elif i == k and j != k and j != k + 1:
                return self.astar[k+1, j]
            elif j == k and i != k and i != k + 1:
                return self.astar[i, k+1]
            elif i != k and i != k+1 and j != k and j != k+1:
                return self.astar[i, j]

    def phi_b_ext(self, a):
        """
        phi_b_ext returns phi_b_ext(a_ij)
        """

        e = self.phi_ext_sk(self.w[0], a)

        for i in self.w[1:]:
            new_e = 0
            for k, v in e._FreeAlgebraElement__monomial_coefficients.iteritems():
                new_g = 1
                for g, p in k:
                    new_g = new_g*self.phi_ext_sk(i, g)
                    if p != 1:
                        print "POW!!!!"
                new_g = v*new_g
                new_e = new_e + new_g
            e = new_e

        return e

    def linearized(self, l):

        if type(l) == sage.rings.integer.Integer:
            return True

        for k, v in l._FreeAlgebraElement__monomial_coefficients.iteritems():
            if (len(k) == 0):
                pass
            elif len(k) == 1:
                if k._element_list[0][1] != 1:
                    return False
            else:
                return False

        return True

    def linearize(self, l):
        s = l

        while not self.linearized(s):

            r = 0

            for k, v in s._FreeAlgebraElement__monomial_coefficients.iteritems():
                e = 0
                els = k._element_list

                if len(els) > 1:
                    g1 = self.alg.gen(els[0][0])
                    if els[0][1] != 1:
                        g2 = self.alg.gen(els[0][0])
                        e = -2*g1 - 2*g2 - 4
                        e = e*(g2**(els[0][1]-2))
                        for g, p in els[1:]:
                            e = e*(self.alg.gen(g)**p)
                    else:
                        g2 = self.alg.gen(els[1][0])
                        e = -2*g1 - 2*g2 - 4
                        for g, p in els[2:]:
                            e = e*(self.alg.gen(g)**p)

                elif els != [] and els[0][1] != 1:
                    e = -4 - 4*self.alg.gen(els[0][0])
                    e = e*(self.alg.gen(els[0][0])**(els[0][1]-2))
                elif els != []:
                    e = self.alg.gen(els[0][0])
                else:
                    e = 1

                r = r + v*e

            s = r

        l = 0
        one = self.alg.monoid().one_element()

        for k, v in s._FreeAlgebraElement__monomial_coefficients.iteritems():
                if k != one:
                    l = l + v*FreeAlgebraElement(self.alg, k)
        return l

    def linearizem(self, m):

        a = [[self.linearize(m[i, j]) for j in xrange(0, self.s)] for i in xrange(0, self.s)]

        return matrix(self.alg, self.s, self.s, a)

    def phi_l_b_ij(self, i, j):

        a = self.phi_b_ext(self.astar[i, 0])
        t = self.astar[j, 0]
        r = 0
        add = False

        for k, v in a._FreeAlgebraElement__monomial_coefficients.iteritems():
            e = 1
            for g, p in k:
                if g != t:
                    e = e*g
                elif g == t:
                    add = True
            if add:
                if type(e) != int:
                    r = r + v*FreeAlgebraElement(self.alg, e)
                else:
                    r = r + v*e
                add = False

        return r

    def phi_l_b(self):

        if self.phi_l_bm:
            return self.phi_l_bm

        a = [[self.phi_l_b_ij(i, j) for j in xrange(1, self.s+1)]
             for i in xrange(1, self.s+1)]
        self.phi_l_bm = matrix(self.alg, a)

        return self.phi_l_bm

    def phi_r_b_ij(self, i, j):

        a = self.phi_b_ext(self.astar[0, j])
        t = self.astar[0, i]
        r = 0
        add = False

        for k, v in a._FreeAlgebraElement__monomial_coefficients.iteritems():
            e = 1
            for g, p in k:
                if g != t:
                    e = e*g
                elif g == t:
                    add = True
            if add:
                if type(e) != int:
                    r = r + v*FreeAlgebraElement(self.alg, e)
                else:
                    r = r + v*e
                add = False

        return r

    def phi_r_b(self):

        if self.phi_r_bm:
            return self.phi_r_bm

        a = [[self.phi_r_b_ij(i, j) for j in xrange(1, self.s+1)] 
             for i in xrange(1, self.s+1)]
        self.phi_r_bm = matrix(self.alg, a)

        return self.phi_r_bm

    def diffa(self):

        return matrix(self.alg, self.s, self.s)
    
    @profile
    def diffb(self):

        p = self.phi_l_b()

        if self.linear:
            m = (self.i - p)*self.a
            return self.linearizem(m)

        if self.polyring:
            phi_b = p*self.a*self.phi_r_b()
            return self.a - (self.cl*phi_b*cli)

        return (self.i - p)*self.a

    def diffc(self):

        p = self.phi_r_b()

        if self.polyring:
            return self.ahat - (self.cl*p*self.a)

        if self.linear:
            m = self.a*(self.i - p)
            return self.linearizem(m)

        return self.a*(self.i - p)

    def diffd(self):

        pr = self.phi_r_b()

        if self.polyring:
            return self.a - (self.ahat*pr*self.cli)

        pl = self.phi_l_b()
        return (self.b*(self.i - pr)) - ((self.i - pl)*self.c)

    def diffe(self):
        
        if self.polyring:
            return self.bhat - self.c - (self.cl*self.phi_l_b()*self.d)

        return self.b + (self.phi_l_b()*self.c)

    def difff(self):
        
        if self.polyring:
            return self.b - self.d - (self.c*self.phi_r_b()*self.cli)

        return matrix(self.alg, self.s, self.s)

    def zero_homology(self):

        db = self.diffb()
        dc = self.diffc()

        cdict = []
        g = iter(count(0))
        for e in list(self.alg.monoid().gens()):
            s = str(e)
            i = s[1]
            j = s[2]
            
            if i == 's' or i == '0' or j == 's' or j == '0' or i == j:
                continue
            if s[0] != 'a':
                break
            else:
                cdict.append((e, g.next()))
        
        cdict = dict(cdict)
        r = self.s*self.s*2
        c = g.next()
        p = matrix(braid.alg, r, c)

        row = 0

        for e in CartesianProduct(xrange(0, self.s), xrange(0, self.s)):
            i, j = e
            diff =  db[i, j]
            for k, v in diff._FreeAlgebraElement__monomial_coefficients.iteritems():
                p[row, cdict[k]] = v
            row = row + 1

        for e in CartesianProduct(xrange(0, self.s), xrange(0, self.s)):
            i, j = e
            diff =  dc[i, j]
            for k, v in diff._FreeAlgebraElement__monomial_coefficients.iteritems():
                p[row, cdict[k]] = v
            row = row + 1

        return p

# braid = Braid(2, [1, 1, 1], polyring=False, linear=True, diag=-2)
# p = braid.zero_homology()

w = [2, 1, 3, 2]*3
w.append(1)
braid = Braid(4, w, polyring=False, diag=-2, linear=True)

braid.diffb()
# print
# print braid.diffb()
# print
# print braid.diffc()
# print
# print braid.diffd()
# print
# print braid.diffe()


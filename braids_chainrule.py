#!/usr/bin/sage -python
from collections import Counter
from itertools import count
from sage.all import *


class Braid:

    def __init__(self, strands, word, polyring=False, linear=True):
        self.s = strands
        self.w = word
        self.w.reverse()
        self.polyring = polyring
        self.linear = linear
        g = self.gens(star=True)

        if self.polyring and not linear:
            r = PolynomialRing(ZZ, 'lmb, lmbinv, mu, muinv, U, Uinv')
            ideal = r.ideal(r.gen(0)*r.gen(1)-1, r.gen(2)*r.gen(3)-1, 
                            r.gen(4)*r.gen(5)-1)
            self.br = QuotientRing(r, ideal, 'lmb, lmbinv, mu, muinv, U, Uinv')
            self.br.inject_variables()
            self.alg = FreeAlgebra(self.br, len(g.split(',')), g)
        elif self.polyring and linear:
            r = PolynomialRing(ZZ, 'mu, muinv')
            ideal = r.ideal(r.gen(0)*r.gen(1)-1)
            self.br = QuotientRing(r, ideal, 'mu, muinv')
            self.br.inject_variables()
            self.alg = FreeAlgebra(self.br, len(g.split(',')), g)
        else:
            self.alg = FreeAlgebra(ZZ, len(g.split(',')), g)
        self.alg.inject_variables()

        g = iter(range(len(self.alg.gens())))

        m = [[self.alg.gen(g.next()) for j in range(0, self.s+1)] 
             for i in range(0, self.s+1)]

        self.astar = matrix(self.alg, self.s+1, self.s+1, m)

        m = [[self.alg.gen(g.next()) for j in range(0, self.s)] 
             for i in range(0, self.s)]
        self.b = matrix(self.alg, self.s, self.s, m)
        
        m = [[self.alg.gen(g.next()) for j in range(0, self.s)]
             for i in range(0, self.s)]
        self.c = matrix(self.alg, self.s, self.s, m)

        m = [[self.alg.gen(g.next()) for j in range(0, self.s)]
             for i in range(0, self.s)]
        self.d = matrix(self.alg, self.s, self.s, m)

        if not self.polyring and linear:
            m = []
            for i in range(1, self.s+1):
                n = []
                for j in range(1, self.s+1):
                    if i == j:
                        n.append(-2)
                    else:
                        n.append(self.astar[i, j])
                m.append(n)
            self.a = matrix(self.alg, self.s, self.s, m)

            m = identity_matrix(self.alg, self.s)
            for i in range(0, self.s):
                m[i, i] = self.alg.gen(g.next())

            self.e = matrix(self.alg, self.s, self.s, m)

        elif self.polyring:
            m = []
            for i in range(1, self.s+1):
                n = []
                for j in range(1, self.s+1):
                    if i < j:
                        n.append(-1*mu*self.astar[i, j])
                    elif i == j:
                        n.append(mu-1)
                    elif i > j:
                        n.append(self.astar[i, j])
                m.append(n)
            self.a = matrix(self.alg, self.s, self.s, m)

            m = []
            for i in range(0, self.s):
                n = []
                for j in range(0, self.s):
                    if i > j:
                        n.append(self.b[i, j])
                    elif i == j:
                        n.append(0)
                    elif i < j:
                        n.append(-mu*self.b[i, j])
                m.append(n)

            self.b = matrix(self.alg, self.s, self.s, m)

            self.cl = identity_matrix(self.alg, self.s)
            self.cli = identity_matrix(self.alg, self.s)

            self.writhe = 0
            for k, v in Counter(self.w).iteritems():
                if k > 0:
                    self.writhe = self.writhe + v
                elif k < 0:
                    self.writhe = self.writhe - v

            if not self.linear:
                m = []
                for i in range(1, self.s+1):
                    n = []
                    for j in range(1, self.s+1):
                        if i < j:
                            n.append(-1*mu*U*self.astar[i, j])
                        elif i == j:
                            n.append(mu*U-1)
                        elif i > j:
                            n.append(self.astar[i, j])
                    m.append(n)

                self.ahat = matrix(self.alg, self.s, self.s, m)

                m = []
                for i in range(0, self.s):
                    n = []
                    for j in range(0, self.s):
                        if i > j:
                            n.append(self.b[i, j])
                        elif i == j:
                            n.append(0)
                        elif i < j:
                            n.append(U*self.b[i, j])
                    m.append(n)

                self.bhat = matrix(self.alg, self.s, self.s, m)

                if (self.writhe-self.s+1)/2 < 0:
                    self.cl[0, 0] = lmb*((-1*muinv)**self.writhe)*(U**(-1*(self.writhe-self.s+1)/2))
                    self.cli[0, 0] = lmbinv*((-1*mu)**self.writhe)*(Uinv**(-1*(self.writhe-self.s+1)/2))
                else:
                    self.cl[0, 0] = lmb*((-1*muinv)**self.writhe)*(Uinv**((self.writhe-self.s+1)/2))
                    self.cli[0, 0] = lmbinv*((-1*mu)**self.writhe)*(U**((self.writhe-self.s+1)/2))

            elif self.linear:

                m = []
                for i in range(1, self.s+1):
                    n = []
                    for j in range(1, self.s+1):
                        if i < j:
                            n.append(-1*mu*self.astar[i, j])
                        elif i == j:
                            n.append(mu-1)
                        elif i > j:
                            n.append(self.astar[i, j])
                    m.append(n)

                self.ahat = matrix(self.alg, self.s, self.s, m)

                m = []
                for i in range(0, self.s):
                    n = []
                    for j in range(0, self.s):
                        if i > j:
                            n.append(self.b[i, j])
                        elif i == j:
                            n.append(0)
                        elif i < j:
                            n.append(-mu*self.b[i, j])
                    m.append(n)

                self.bhat = matrix(self.alg, self.s, self.s, m)

                self.cl[0, 0] = (-1*muinv)**self.writhe
                self.cli[0, 0] = (-1*mu)**self.writhe

            m = [[self.alg.gen(g.next()) for j in range(0, self.s)]
                 for i in range(0, self.s)]
            self.e = matrix(self.alg, self.s, self.s, m)

            m = [[self.alg.gen(g.next()) for j in range(0, self.s)]
                 for i in range(0, self.s)]
            self.f = matrix(self.alg, self.s, self.s, m)

        self.phi_l_bm = None
        self.phi_r_bm = None

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

        g = a + b + c + d

        if self.polyring:
            for i in range(1, self.s+1):
                for j in range(1, self.s+1):
                    e = e + 'e' + str(i) + str(j) + ','
            for i in range(1, self.s+1):
                for j in range(1, self.s+1):
                    f = f + 'f' + str(i) + str(j) + ','
            g = g + e + f[:-1]
        else:
            for i in range(1, self.s+1):
                e = e + 'e' + str(i) + str(i) + ','
            g = g + e[:-1]

        return g

    def phi_ext_sk(self, k, a):
        """
        phi_ext_sk returns phi^ext_sigma_z(a_xy) 
        """
        try:
            i, j = self.gdicta[a]
        except KeyError:
            return a

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
                    g1 = (self.phi_ext_sk(i, g)**p)
                    new_g = new_g*g1
                new_e = new_e + v*new_g
            e = new_e

        return e

    def phi_b(self, w, a):
        """
        phi_b returns phi_b(a_ij)
        """

        e = self.phi_ext_sk(w[0], a)

        for i in w[1:]:
            new_e = 0
            for k, v in e._FreeAlgebraElement__monomial_coefficients.iteritems():
                new_g = 1
                for g, p in k:
                    g1 = (self.phi_ext_sk(i, g)**p)
                    new_g = new_g*g1
                new_e = new_e + v*new_g
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

    def linearize(self, l, constants=False):
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
                        if els[1][1] != 1:
                            e = e*(g2**(els[1][1]-1))
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

        if not constants:

            l = 0
            one = self.alg.monoid().one_element()

            for k, v in s._FreeAlgebraElement__monomial_coefficients.iteritems():
                    if k != one:
                        for g, p in k._element_list:
                            l = l + v*(self.alg.gen(g)**p)
            return l

        else:
            return s

    def linearizem(self, m, constants=False):

        a = [[self.linearize(m[i, j], constants) for j in range(0, self.s)] for i in range(0, self.s)]

        return matrix(self.alg, self.s, self.s, a)

    def phi_l_sk(self, k):
        m = copy(self.i)
        
        if k > 0:
            m[k-1:k+1, k-1:k+1] = [[-1*self.a[k, k-1], -1], [1, 0]]
        elif k < 0:
            k = k*-1
            m[k-1:k+1, k-1:k+1] = [[0, 1], [-1, -1*self.a[k-1, k]]]

        return m

    def phi_l_b_help(self, w, r):
        b1 = w[:-1]
        b2 = w[-1]

        if b1 != []:
            p = self.phi_l_sk(b2)
            rb1 = list(b1)
            rb1.reverse()

            if b2 > 0:
                i = b2
            elif b2 < 0:
                i = (b2*-1)+1
                
            p[i-1, i-1] = -1*self.phi_b(rb1, -1*p[i-1, i-1])
            return self.phi_l_b_help(b1, r+[p])

        return r+[self.phi_l_sk(b2)]

    def phi_l_b(self):

        return reduce(lambda x, y: x*y, self.phi_l_b_help(self.w, []))

    def phi_r_sk(self, k):
        m = copy(self.i)
        
        if k > 0:
            m[k-1:k+1, k-1:k+1] = [[-1*self.a[k-1, k], 1], [-1, 0]]
        elif k < 0:
            k = k*-1
            m[k-1:k+1, k-1:k+1] = [[0, -1], [1, -1*self.a[k, k-1]]]
        return m

    def phi_r_b_help(self, w, r):
        b1 = w[:-1]
        b2 = w[-1]

        if b1 != []:
            p = self.phi_r_sk(b2)
            rb1 = list(b1)
            rb1.reverse()
            if b2 > 0:
                i = b2
            elif b2 < 0:
                i = (b2*-1)+1

            p[i-1, i-1] = -1*self.phi_b(rb1, -1*p[i-1, i-1])
            return self.phi_r_b_help(b1, r+[p])

        return r+[self.phi_r_sk(b2)]

    def phi_r_b(self):
        r = self.phi_r_b_help(self.w, [])
        r.reverse()

        return reduce(lambda x, y: x*y, r)

    def diffa(self):

        return matrix(self.alg, self.s, self.s)
    
    def diffb(self):

        p = self.phi_l_b()

        if self.linear and not self.polyring:
            m = (self.i - p)*self.a
            return self.linearizem(m)
        elif self.linear and self.polyring:
            phi_b = p*self.a*self.phi_r_b()
            m = self.a - (self.cl*phi_b*self.cli)
            return self.linearizem(m)
        elif self.polyring:
            phi_b = p*self.a*self.phi_r_b()
            return self.a - (self.cl*phi_b*self.cli)
        else:
            return (self.i - p)*self.a

    def diffc(self):

        p = self.phi_r_b()

        if self.linear and not self.polyring:
            m = self.a*(self.i - p)
            return self.linearizem(m)
        elif self.linear and self.polyring:
            m = self.ahat - (self.cl*p*self.a)
            return self.linearizem(m)
        elif self.polyring and not self.linear:
            return self.ahat - (self.cl*p*self.a)
        else:
            return self.a*(self.i - p)

    def diffd(self):

        pr = self.phi_r_b()

        if self.linear and not self.polyring:
            pl = self.phi_l_b()
            m = (self.b*(self.i - pr)) - ((self.i - pl)*self.c)
            return self.linearizem(m)
        elif self.linear and self.polyring:
            m = self.a - (self.ahat*pr*self.cli)
            return self.linearizem(m)
        elif self.polyring and not self.linear:
            return self.a - (self.ahat*pr*self.cli)
        else:
            return (self.b*(self.i - pr)) - ((self.i - pl)*self.c)

    def diffe(self):
        
        pl = self.phi_l_b()

        if self.linear and not self.polyring:
            m =  self.b + (pl*self.c)
            return self.linearizem(m)
        elif self.linear and self.polyring:
            m = self.bhat - self.c - (self.cl*pl*self.d)
            return self.linearizem(m)
        elif self.polyring and not self.linear:
            return self.bhat - self.c - (self.cl*pl*self.d)
        else:
            return self.b + (self.phi_l_b()*self.c)

    def difff(self):
        pr = self.phi_r_b()

        if self.polyring and self.linear:
            m = self.b - self.d - (self.c*pr*self.cli)
            return self.linearizem(m)
        elif self.polyring and not self.linear:
            return self.b - self.d - (self.c*pr*self.cli)
        else:
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
        
        if not self.polyring:
            p = matrix(r, c)
        elif self.polyring:
            p = matrix(self.br, r, c)

        row = 0

        for e in CartesianProduct(range(0, self.s), range(0, self.s)):
            i, j = e
            diff =  db[i, j]
            for k, v in diff._FreeAlgebraElement__monomial_coefficients.iteritems():
                p[row, cdict[k]] = v
            row = row + 1

        for e in CartesianProduct(range(0, self.s), range(0, self.s)):
            i, j = e
            diff =  dc[i, j]
            for k, v in diff._FreeAlgebraElement__monomial_coefficients.iteritems():
                p[row, cdict[k]] = v
            row = row + 1

        if self.polyring:
            return p
        else:
            return p.smith_form()[0]

    def first_homology(self):

        db = self.diffb()
        dc = self.diffc()
        dd = self.diffd()
        de = self.diffe()

        cdict = []
        g = iter(counpppt(0))
        for e in list(self.alg.monoid().gens())[(self.s+1)**2:]:

            s = str(e)
            i = s[1]
            j = s[2]

            if s[0] != 'b' and s[0] != 'c':
                break
            else:
                cdict.append((e, g.next()))

        cdict = dict(cdict)
        r = (self.s**2)+(self.s**2)-2
        c = g.next()
        p = matrix(r, c)

        row = 0

        for e in CartesianProduct(range(0, self.s), range(0, self.s)):
            i, j = e
            diff =  db[i, j]
            for k, v in diff._FreeAlgebraElement__monomial_coefficients.iteritems():
                p[row, cdict[k]] = v
            row = row + 1

        for e in CartesianProduct(range(0, self.s), range(0, self.s)):
            i, j = e
            diff =  dc[i, j]
            for k, v in diff._FreeAlgebraElement__monomial_coefficients.iteritems():
                p[row, cdict[k]] = v
            row = row + 1

        return p


def satellite(w):
    r = []

    for i in w:
        a = [-1*(2*i+1), -1*(2*i-1), 2*i, 2*i+1, 2*i-1, 2*i]
        r.extend(a)

    return r

braid = Braid(4, satellite([1]*5), polyring=False, linear=True)
p = braid.zero_homology()
print p.str()
#braid = Braid(4, [-3, -3, 2, 3, 1, 2], polyring=False, linear=True)
# braid = Braid(2, [1, 1, 1], polyring=False, linear=True)
# braid = Braid(2, [1, 1, 1], polyring=True, linear=False)
# braid = Braid(3, [1, -2, 1, -2], polyring = False, linear=True)
# braid = Braid(4, [1, -2, 3, -2, -1, 2, 3], polyring = False, linear=True)
# p = braid.zero_homology()
# print p.str()

# w = [2, 1, 3, 2]*3
# w.append(1)
# braid = Braid(4, w, polyring=False, linear=True)
# p = braid.zero_homology()
# print p.str()

# braid = Braid(4, satellite([1]), polyring=False, linear=True)
# p = braid.phi_l_b()


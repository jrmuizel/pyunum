"""
Copyright (c) Taylor & Francis Group, 2014.
Copyright (c) Jeff Muizelaar, 2015.

Permission is hereby granted, free of charge, to any person obtaining \
a copy of this software and associated documentation files (the \
"Software"), to deal in the Software without restriction, including \
without limitation the rights to use, copy, modify, merge, publish, \
distribute, sublicense, and/or sell copies of the Software, and to \
permit persons to whom the Software is furnished to do so, subject to \
the following conditions:

The above copyright notice and this permission notice shall be \
included in all copies or substantial portions of the software.

THE SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, \
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF \
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND \
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS \
BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER LIABILITY, WHETHER IN AN \
ACTION OR CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN \
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE \
SOFTWARE.
"""

from __future__ import division

import math
open = True
closed = False
from fractions import *
import unum_config
e = unum_config.e
f = unum_config.f
esizesize = e
fsizesize = f
esizemax = 2**esizesize
fsizemax = 2**fsizesize
utagsize = 1 + f + e
maxubits = 1 + esizemax + fsizemax + utagsize
ubitmask = 1 << (utagsize - 1)
fsizemask = (1 << f) - 1
esizemask = (ubitmask - 1) - fsizemask
efsizemask = esizemask | fsizemask
utagmask = ubitmask | efsizemask
ulpu = 1 << utagsize
smallsubnormalu = efsizemask + ulpu
smallnormalu = efsizemask + (1 << (maxubits - 1 - esizemax))
signbigu = 1 << maxubits - 1
posinfu = signbigu - 1 - ubitmask
maxrealu = posinfu - ulpu
minrealu = maxrealu + signbigu
neginfu = posinfu + signbigu
negbigu = neginfu - ulpu
qNaNu = posinfu + ubitmask
sNaNu = neginfu + ubitmask
NaN = float("nan")
if utagsize == 1:
    negopeninfu = 0xd # 1101
else:
    negopeninfu = 0xf << (utagsize - 1)

if utagsize == 1:
    posopeninfu = 0x5 # 0101
else:
    posopeninfu = 0x7 << (utagsize - 1)
negopenzerou = 0x9 << (utagsize - 1)

maxreal = Fraction(2**(2**(esizemax-1)) * (2**fsizemax - 1), (2**(fsizemax-1)))
smallsubnormal = Fraction(2**2, 2**(2**(esizemax-1) + fsizemax))

def fractionalPart(x):
    return x - math.trunc(x)


def integerPart(x):
    return math.trunc(x)

def flatten(container):
    for i in container:
        if isinstance(i, list) or isinstance(i, tuple):
            for j in flatten(i):
                yield j
            else:
                yield i

def transpose(l):
    return tuple([tuple(i) for i in zip(*l)])

# XXX: this is a bad implementation, but at least it avoids division
def binomial(n, k):
    if k == 0:
        return 1
    if n == k:
        return 1
    return binomial(n-1, k-1) + binomial(n-1, k)

def log(x, base):
    return math.log(x, base)

def floor(x):
    if math.isinf(x):
        return x
    return int(math.floor(x))

def ceil(x):
    if math.isinf(x):
        return x
    return int(math.ceil(x))


def unumQ(x):
    if isinstance(x, (int, long)):
        if x >= 0 and x <= sNaNu:
            return True
    return False

def uQ(x):
    return unumQ(x) or uboundQ(x)

def unum2g(u):
    if unumQ(u):
        if u == qNaNu or u == sNaNu:
            return ((NaN, NaN), (open, open))
        x = u2f(exact(u))
        y = u2f(exact(u) + ulpu)
        if exQ(u):
            return ((x, x), (closed, closed))
        elif u == (bigu(u) + ubitmask):
            return ((big(u), float('inf')), (open, open))
        elif u == signmask(u) + bigu(u) + ubitmask:
            return ((float('-inf'), -big(u)), (open, open))
        elif sign(u) == 1:
            return ((y, x), (open, open))
        else:
            return ((x, y), (open, open))

def floatQ(x):
    if math.isnan(x) or math.isinf(x) or isinstance(x, (Fraction, int, long, float)):
        return True
    return False

def esizeminus1(u):
    if unumQ(u):
        return (u & esizemask) >> fsizesize

def esize(u):
    if unumQ(u):
        return 1 + esizeminus1(u)

def fsizeminus1(u):
    if unumQ(u):
        return u & fsizemask

def fsize(u):
    if unumQ(u):
        return 1 + fsizeminus1(u)

def numbits(u):
    if unumQ(u):
        return 1 + esize(u) + fsize(u) + utagsize

def signmask(u):
    if unumQ(u):
        return 1 << (numbits(u) - 1)

def hiddenmask(u):
    if unumQ(u):
        return 1 << (fsize(u) + utagsize)
def fracmask(u):
    if unumQ(u):
        return ((1 << fsize(u)) - 1) << utagsize

def expomask(u):
    if unumQ(u):
        return ((1 << esize(u)) - 1) << (fsize(u) + utagsize)
def floatmask(u):
    if unumQ(u):
        return signmask(u) + expomask(u) + fracmask(u)

def bias(u):
    if unumQ(u):
        return 2**(esizeminus1(u)) - 1
def sign(u):
    if unumQ(u):
        return (u & signmask(u) > 0)
def expo(u):
    if unumQ(u):
        return (u & expomask(u)) >> (utagsize + fsize(u))
def hidden(u):
    if unumQ(u):
        return expo(u) > 0
def frac(u):
    if unumQ(u):
        return (u & fracmask(u)) >> utagsize

def inexQ(u):
    if unumQ(u):
        return (u & ubitmask) > 0

def exQ(u):
    if unumQ(u):
        return (u & ubitmask) == 0

def exact(u):
    if unumQ(u):
        if inexQ(u):
            return u ^ ubitmask
        else:
            return u

def expovalue(u):
    if unumQ(u):
        return expo(u) - bias(u) + 1 - hidden(u)

def u2f(u):
    if unumQ(u) and exQ(u):
        if u == posinfu:
            return float("inf")
        if u == neginfu:
            return -float("inf")
        k = ((-1)**sign(u))*(Fraction(2,1)**expovalue(u))*(hidden(u) + Fraction(frac(u), 2**fsize(u)))
        assert not isinstance(k, float)
        return k

# Biggest unum possible with identical utag contents.
def bigu(u):
    if unumQ(u):
        return expomask(u) + fracmask(u) + (efsizemask & u) - ulpu * ((u & efsizemask) == efsizemask)

# Biggest numerical value representable with identical utag contents.
def big(u):
    if unumQ(u):
        return u2f(bigu(u))

def gQ(x):
    if isinstance(x, tuple):
        if len(x) == 2:
            if isinstance(x[0], tuple):
                if len(x[0]) == 2:
                    if isinstance(x[1], tuple):
                        if len(x[1]) == 2:
                            if isinstance(x[1][0], bool) and isinstance(x[1][1], bool) and floatQ(x[0][0]) and floatQ(x[0][1]):
                                if math.isnan(x[0][0]) or math.isnan(x[0][1]):
                                    return True
                                if (x[0][0] == x[0][1] and not x[1][0] and not x[1][1]) or x[0][0] < x[0][1]:
                                    return True
    return False

def uboundQ(x):
    if isinstance(x, tuple):
        if len(x) == 1 or len(x) == 2:
            xL = x[0]
            xR = x[-1]
            if unumQ(xL) and unumQ(xR):
                gL, gR = (unum2g(xL), unum2g(xR))
                if len(x) == 1 or (xL == qNaNu or xL == sNaNu or xR == qNaNu or xR == sNaNu) or \
                        ((gL[0][0] < gR[0][1]) or (gL[0][0] == gR[0][1] and (exQ(xL) and exQ(xR)))):
                            return True
    return False

# should this look at negative inf?
def scale(x):
    if floatQ(x) and x != float("inf") and not math.isnan(x):
        if x == 0:
            return 0
        else:
            return floor(log(abs(x),2))

def ne(x):
    if floatQ(x) and x != float("inf") and not math.isnan(x):
        if x == 0 and scale(x) == 1:
            return 1
        else:
            return int(math.ceil(log(1+abs(scale(x)-1),2))) + 1

def nbits(u):
    if uQ(u):
        if unumQ(u):
            ub = (u,)
        else:
            ub = u
        if len(ub) == 2:
            return 1 + numbits(ub[0]) + numbits(ub[1])
        else:
            return 1 + numbits(ub[0])

def ubound2g(ub):
    if uboundQ(ub):
        uL = ub[0]
        uR = ub[-1]
        if uL == qNaNu or uL == sNaNu or uR == qNaNu or uR == sNaNu:
            return ((float("nan"), float("nan")), (open, open))
        gL, gR = (unum2g(uL), unum2g(uR))
        return ((gL[0][0], gR[0][1]), (gL[1][0], gR[1][1]))

def u2g(u):
    if uQ(u):
        if unumQ(u):
            return unum2g(u)
        else:
            return ubound2g(u)
    #raise TypeError(u)

def x2u(x):
    if floatQ(x):
        # Exceptional nonnumeric values:
        if math.isnan(x):
            return qNaNu
        elif x == float("inf"):
            return posinfu
        elif x == float("-inf"):
            return neginfu
        # Magnitudes too large to represent:
        elif abs(x) > maxreal:
            return maxrealu + ubitmask + (signbigu if x < 0 else 0)
        # Zero is a special case. The smallest unum for it is just 0:
        elif x == 0:
            return 0
        # Magnitudes too small to represent become "inexact zero" with
        # the maximum exponent and fraction field sizes:
        elif abs(x) < smallsubnormal:
            return utagmask + (signbigu if x < 0 else 0)
        # For subnormal numbers, divide by the ULP value to get the
        # fractional part. The while loop strips off the trailing bits.
        elif abs(x) < u2f(smallnormalu):
            y = abs(x) / smallsubnormal
            y = (signbigu if x < 0 else 0) + efsizemask + (ubitmask if y != math.floor(y) else 0) + (floor(y) << utagsize)
            # XXX is this right?
            #assert(False)
            while ((3 << (utagsize - 1)) & y) == 0:
                y = (y - (efsizemask & y))//2 + (efsizemask & y) - 1
            return y
        # All remaining cases are in the normalized range.
        else:
            n = 0
            y = abs(x) / 2**(scale(x))
            n = 0
            while math.floor(y) != y and n < fsizemax:
                n += 1
                y *= 2
            if y == math.floor(y): # then the value is representable
                # exactly. Fill in the fields from right to left:
                # Size of fraction field,
                # fits in the rightmost fsizesize bits...
                y1 = n - (n > 0)
                # Size of exponent field minus 1 fits in the esizesize bits...
                y1 += ((ne(x) - 1) << fsizesize)
                # Significant bits after hidden bit fits left of the unum tag bits...
                y1 += 0 if n == 0 else ((floor(y) - 2**scale(y)) << utagsize)
                # Value of exponent bits, adjusted for bias...
                y1 += (scale(x) + 2**(ne(x) - 1) - 1) << (utagsize + n + (n == 0))
                # If negative, add the sign bit
                y1 += (1 << (utagsize + n + (n == 0) + ne(x))) if x < 0 else 0
                
                # if a number is more concise as a subnormal, make it one
                z = None
                try:
                    z = log(1 - log(abs(x), 2), 2)
                except:
                    pass
                if z and z.is_integer() and z >= 0:
                    z = long(z)
                    return (z << fsizesize) + ulpu + (x < 0) * signmask(z << fsizesize)
                else:
                    return y1
            else:
                # inexact. Use all available fraction bits
                z = math.ceil(abs(x) / 2**(scale(x) - fsizemax)) * 2 ** (scale(x)-fsizemax)
                n = max(ne(x), ne(z))
                # All bits on for the fraction size, since we're using the maximum
                y1 = fsizemask
                # Store the exponent size minus 1 in the exponent size field
                y1 += (n - 1) << fsizesize
                # Back off by one ULP and make it inexact
                y1 += ubitmask - ulpu
                # Fraction bits are the ones to the left of the
                # binary point after removing hidden bit and scaling
                y1 += floor((z / 2**scale(z) - 1) * 2**fsizemax) << utagsize
                # Exponent value goes in the exponent field
                y1 += (scale(z) + 2**(n-1) - 1) << (utagsize + fsizemax)

                if x < 0:
                    y1 += signmask(y1)
                return y1

# View a float as a decimal, using as many digits as needed to be exact.
def autoN(x):
    if math.isnan(x) or x == 0 or x == float('inf'):
        return str(x)
    if x < 0:
        return "-" + autoN(-x)
    y = log(x.denominator, 2)
    if y == 0:
        return str(x).zfill(1 + floor(log(x, 10)))
    if isinstance(x, Fraction) and y == floor(y):
        y = x - floor(x)
        z = floor(log(y.denominator, 2))
        return str(floor(x)) + "." + str(y*10**z).zfill(z)
    return "?"

def view(g):
    if gQ(g) or unumQ(g) or uboundQ(g):
        ((L, R), (LQ, RQ)) = g if gQ(g) else u2g(g)
        if math.isnan(L) or math.isnan(R):
            return "NaN"
        if L == R and not LQ and not RQ:
            return autoN(L)
        if L < R:
            return ("(" if LQ else "[") + autoN(L) + ", " + autoN(R) + (")" if RQ else "]")
        return "NaN"

def plusg(x, y):
    ((xlo, xhi), (xlob, xhib)) = x
    ((ylo, yhi), (ylob, yhib)) = y
    if math.isnan(xlo) or math.isnan(xhi) or math.isnan(ylo) or math.isnan(yhi):
        return ((NaN, NaN), (open, open))
    if xlo == float('-inf') and not xlob:
        (sumleft, openleft) = (NaN, open) if ylo == float("inf") and not ylob else (float("-inf"), closed)
    elif ylo == float('-inf') and not ylob:
        (sumleft, openleft) = (NaN, open) if xlo == float("inf") and not xlob else (float("-inf"), closed)
    elif (xlo == float('inf') and not xlob) or (ylo == float('inf') and not ylob):
        (sumleft, openleft) = (float('inf'), closed)
    elif xlo == float('-inf'):
        (sumleft, openleft) = (float('inf'), closed) if ylo == float('inf') and not ylob else (float('-inf'), open)
    elif ylo == float('-inf'):
        (sumleft, openleft) = (float('inf'), closed) if xlo == float('inf') and not xlob else (float('-inf'), open)
    else:
        (sumleft, openleft) = ((xlo + ylo), (xlob or ylob))
   
    if xhi == float('-inf') and not xhib:
        (sumright, openright) = (NaN, open) if yhi == float("inf") and not yhib else (float("-inf"), closed)
    elif yhi == float('-inf') and not yhib:
        (sumright, openright) = (NaN, open) if xhi == float("inf") and not xhib else (float("-inf"), closed)
    elif (xhi == float('inf') and not xhib) or (yhi == float('inf') and not yhib):
        (sumright, openright) = (float('inf'), closed)
    elif xhi == float('-inf'):
        (sumright, openright) = (float('inf'), closed) if yhi == float('inf') and not yhib else (float('-inf'), open)
    elif yhi == float('-inf'):
        (sumright, openright) = (float('inf'), closed) if xhi == float('inf') and not xhib else (float('-inf'), open)
    else:
        (sumright, openright) = ((xhi + yhi), (xhib or yhib))
    return ((sumleft, sumright), (openleft, openright))

def plusu(u, v):
    if uQ(u) and uQ(v):
        w = g2u(plusg(u2g(u), u2g(v)))
        global ubitsmoved, numbersmoved
        ubitsmoved += nbits(u) + nbits(v) + nbits(w)
        numbersmoved += 3
        return w

def negateg(x):
    return ((-x[0][1], -x[0][0]), (x[1][1], x[1][0]))

def minusg(x, y):
    return plusg(x, negateg(y))

def minusu(u, v):
    w = g2u(minusg(u2g(u), u2g(v)))
    global ubitsmoved, numbersmoved
    ubitsmoved += nbits(u) + nbits(v) + nbits(w)
    numbersmoved += 3
    return w

def ubright(xright):
    open = xright[1]
    u = x2u(xright[0])
    x = xright[0]
    if x == float('inf'):
        if open:
            return posopeninfu
        else:
            return posinfu
    if x == 0 and open:
        return negopenzerou
    if u2f(u) == x:
        if open:
            return (u - (ulpu * (x >= 0))) | ubitmask
        else:
            return u
    else:
        return u | (open * ubitmask)

# Find the left half of a ubound (numerical value and open-closed)
def ubleft(xleft):
    open = xleft[1]
    u = x2u(xleft[0])
    x = xleft[0]
    if x == float('-inf'):
        if open:
            return negopeninfu
        else:
            return neginfu
    if u2f(u) == x:
        if open:
            return (u - (ulpu * (x < 0))) | ubitmask
        else:
            return u
    else:
        return u | (open * ubitmask)

def negateu(u):
    if uQ(u):
        if uboundQ(u):
            if len(u) == 1:
                return (x2u(0) if u2g(u[0]) == u2g(0) else signmask(u[0]) ^ u[0],)
            else:
                return (x2u(0) if u2g(u[1]) == u2g(0) else signmask(u[1]) ^ u[1],
                    x2u(0) if u2g(u[0]) == u2g(0) else signmask(u[0]) ^ u[0])
        else:
            return x2u(0) if u2g(u) == u2g(0) else signmask(u) ^ u


def unifypos(ub):
    if uboundQ(ub):
        u = ub[0]
        v = ub[-1]
        # First do trivial case where endpoints express the same value
        if u2g(u) == u2g(v):
            return g2u(u2g(u))
        # Cannot unify if the interval includes exact 0, 1, 2, or 3
        if nnequQ(ub, x2u(0)) or nnequQ(ub, x2u(1)) or nnequQ(ub, x2u(2)) or nnequQ(ub, x2u(3)):
            return ub
        # Refine the endpoints for the tightest possible unification.
        u = promote(x2u(u2g(u)[0][0]), efsizemask)[0] + (ubitmask if inexQ(u) else -ubitmask)
        v = promote(x2u(u2g(v)[0][1]), efsizemask)[0] - (ubitmask if inexQ(v) else -ubitmask)
        if u == v:
            return (u,)
        # If upper bound is open inf and lower bound > maxreal, special handling is needed
        if u2g(v)[0][1] == float("inf") and u2g(v)[1][1]:
            if ltuQ((maxrealu,), (u,)):
                return (maxrealu + ubitmask,)
            # Defmote the left bound until the upper bound is open inf
            while u2g(u)[0][1] < float('inf'):
                if esize(u) > 1:
                    u = demotee(u)
                else:
                    u = demotef(u)
            return (u,)
        # While demoting exponents is possible and still leaves unums within ubound, demote both exponents
        while u != v and (((u2g(demotee(u))[0][0] < u2g(demotee(v))[0][0] and \
            u2g(demotee(u))[0][1] < u2g(demotee(v))[0][1]  < float('inf')))) and esize(u) > 1:
            u = demotee(u)
            v = demotee(v)
        while u != v and frac(v) != frac(u) and fsize(u) > 1:
            u = demotef(u)
            v = demotef(v)
        if u != v and ((floatmask(u) + ubitmask) | u) == ubitmask and ltuQ((v,), (x2u(1),)):
            n = min(esizemax, floor(log(1 - log(u2g(v + ubitmask if exQ[v] else 0)[0,1], 2), 2)))
            return (x2u(2**(-2**n+1)) - ubitmask,)
        else:
            if u == v:
                return (u,)
            else:
                return ub

def unify(ub):
    if uboundQ(ub):
        u = ub[0]
        v = ub[-1]
        if u == qNaNu or u == sNaNu or v == qNaNu or v == sNaNu:
            return (qNaNu,)
        if u == posinfu and v == posinfu:
            return (posinfu,)
        if u == neginfu and v == neginfu:
            return (neginfu,)
        if u == neginfu or u == posinfu or v == neginfu or v == posinfu or \
                ltuQ((u,), (x2u(0),)) and not ltuQ((v,), (x2u(0),)):
                    return ub
        if ltuQ((u,), (x2u(0),)) and ltuQ((v,), (x2u(0),)):
            return negateu(unifypos(negateu(ub)))
        if u2g(u) == u2g(v):
            return (min(u, v),)
        return unifypos(ub)

def g2u(g):
    if gQ(g):
        ulo = x2u(g[0][0])
        uhi = x2u(g[0][1])
        blo = g[1][0]
        bhi = g[1][1]
        # XXX why not use the named version from above?
        if ulo == qNaNu or uhi == qNaNu or g[0][0] > g[0][1] or (g[0][0] == g[0][1] and (blo or bhi)):
            return qNaNu
        elif ulo == uhi and not (blo != bhi):
            return x2u(g[0][0])
        else:
            u1 = ubleft((g[0][0], blo))
            u2 = ubright((g[0][1], bhi))
            if u2g(unify((u1, u2))) == u2g((u1, u2)):
                return unify((u1, u2))
            else:
                return (u1, u2)
    #raise TypeError(g)

# Test if interval g is strictly less than interval h.
def ltgQ(g, h):
    if gQ(g) and gQ(h):
        if math.isnan(g[0][0]) or math.isnan(g[0][1]) or math.isnan(h[0][0]) or math.isnan(h[0][1]):
            return False
        return g[0][1] < h[0][0] or (g[0][1] == h[0][0] and (g[1][1] or h[1][0]))

# Test if ubound or unum u is strictly less than ubound or unum v.
def ltuQ(u, v):
    if uQ(u) and uQ(v):
        return ltgQ(u2g(u), u2g(v))

# Test if interval g is strictly greater than interval h.
def gtgQ(g, h):
    if gQ(g) and gQ(h):
        if math.isnan(g[0][0]) or math.isnan(g[0][1]) or math.isnan(h[0][0]) or math.isnan(h[0][1]):
            return False
        return g[0][0] > h[0][1] or (g[0][0] == h[0][1] and (g[1][0] or h[1][1]))

def gtuQ(u, v):
    if uQ(u) and uQ(v):
        return gtgQ(u2g(u), u2g(v))

# Test if interval g is not nowhere equal to inveral h
def nneqgQ(g, h):
    if gQ(g) and gQ(h):
        if math.isnan(g[0][0]) or math.isnan(g[0][1]) or math.isnan(h[0][0]) or math.isnan(h[0][1]):
            return False
        return not (ltgQ(g, h) or gtgQ(g, h))

# Test if ubound or unum u is not nowhere equal to ubound or unum v
def nnequQ(u, v):
    if uQ(u) and uQ(v):
        return nneqgQ(u2g(u), u2g(v))

# Test if interval g is identical to interval h
def samegQ(g, h):
    if gQ(g) and gQ(h):
        return g == h

# Test if ubound or unum u value is identical to ubound or unum v value.
def sameuQ(u, v):
    if uQ(u) and uQ(v):
        return samegQ(u2g(u), u2g(v))

def intersectg(g, h):
    k = intersectgi(g,h)
    #print k
    return k
def intersectgi(g, h):
    #print 'intersect', g,h
    if gQ(g) and gQ(h):
        glo, ghi = g[0]
        glob, ghib = g[1]
        hlo, hhi = h[0]
        hlob, hhib = h[1]
        if math.isnan(glo) or math.isnan(ghi) or math.isnan(hlo) or math.isnan(hhi):
            return ((float("nan"), float("nan")), (open, open))
        if glo < hlo or (glo == hlo and hlob):
            # left end of g is left of end of h. Three sub-cases to test.
            if ghi < hlo or (ghi == hlo and (ghib or hlob)):
                return ((float("nan"), float("nan")), (open, open))
            # g is completely left of h.
            if ghi < hhi or (ghi == hhi and (ghib or not hhib)):
                return ((hlo, ghi), (hlob, ghib))
            # right part of g overlaps left part of h.
            return ((hlo, hhi), (hlob, hhib)) # h is entire inside g.
        if glo < hhi or (glo == hhi and (not glob and not hhib)):
            # left end of g is inside h. Two sub-cases to test.
            if ghi < hhi or (ghi == hhi and (ghib or not hhib)):
                # g is entirely inside h.
                return ((glo, ghi), (glob, ghib))
            
            # left end of g overlaps right end of h.
            return ((glo, hhi), (glob, hhib))
        # g is entirely to the right of h.
        return ((NaN, NaN), (open, open))


# Add a zero bit to the fraction length of an exact unum, if possible.
def promotef(u):
    if unumQ(u) and exQ(u):
        if fsize(u) < fsizemax:
            return 2 * (floatmask(u) & u) + (utagmask & u) + 1
        return u

# Increase the length of the exponent field of an exact unum, if possible.
def promotee(u):
    if unumQ(u) and exQ(u):
        e = expo(u)
        es = esize(u)
        f = frac(u)
        fs = fsize(u)
        s = signmask(u) & u
        ut = (utagmask & u) + fsizemax
        # If already maximum exponent size, do nothing. This also handles NaN and inf values
        if es == esizemax:
            return u
        # Take care of u = 0 case, ignoring the sign bit. It's simply the new utag.
        if e == 0 and f == 0:
            return ut
        # If normal (nonzero exponent), slide sign bit left, add 2**(es-1), increment esize.
        if e > 0:
            return 2 * s + (e + 2**(es-1)) * hiddenmask(u) + ((hiddenmask(u) - 1) & u) + fsizemax
        # Subnormal. Room to shift and stay subnormal?
        if fs - (floor(log(f, 2)) + 1) >= 2**(es-1):
            return 2 * s + frac(u) * (2**(2**(es-1))) * ulpu + ut
        # Subnormal becomes normal. Trickiest case.
        # The fraction slides left such that the lefmost 1 becomes the hidden bit
        nsigbits = floor(log(f,2)) + 1
        return 2*s + (2**(es-1) + 1 - fs + nsigbits)*(hiddenmask(u)) + \
                (f - 2**nsigbits)*2**(fs-nsigbits+1)*ulpu + (utagmask & u) + fsizemask



# Promote a pair of exact unums to the same esize and fsize
def promote(u, v):
    if unumQ(u) and unumQ(v) and exQ(u) and exQ(v):
        eu = esize(u)
        ev = esize(v)
        fu = fsize(u)
        fv = fsize(v)
        ut = u
        vt = v
        while eu < ev:
            ut = promotee(ut)
            eu += 1
        while ev < eu:
            vt = promotee(vt)
            ev += 1
        while fu < fv:
            ut = promotef(ut)
            fu += 1
        while fv < fu:
            vt = promotef(vt)
            fv += 1
        return (ut, vt)

# Demote the fraction of a unum if possible,
# even if it makes it inexact.
def demotef(u):
    if unumQ(u):
        # Cannot make the fraction any smaller
        if fsize(u) == 1 or u == posinfu or u == neginfu or u == qNaNu or u == sNaNu:
            return u
        # Else shift fraction right one bit.
        return ((u & floatmask(u)) // 2) | ((utagmask & u) - 1)

def demotee(u):
    if unumQ(u):
        es = esize(u)
        mask = int(signmask(u) / 4)
        fm = floatmask(u)
        ut = u & utagmask
        s = signmask(u) & u
        f = frac(u)
        left2 = Fraction((u & (3*mask)), mask)
        # Cannot make the exponent any smaller:
        if es == 1 or u == posinfu or u == neginfu or u == qNaNu or u == sNaNu:
            return u
        if expo(u) == 0:
            f = Fraction(f, (2**(2**(es-2))))
            ibit = ubitmask if fractionalPart(f) > 0 else 0
            # XXX: is the division ok here?
            return ibit | (s // 2 + integerPart(f) * ulpu + ut - fsizemax)
        # If the left two exponent bits are 00
        # (but it's normal, since we fell through the previous test),
        # result switches to subnormal. The exponent after the first
        # two bits joins the fraction like a fixed-point number,
        # before shifting the fraction to the right. The
        # new exponent is zero, of course.
        if left2 == 0:
            f = Fraction(2**fsize(u) + f, Fraction(2**(2**(es-2)+1), 2**expo(u)))
            ibit = ubitmask if fractionalPart(f) > 0 else 0
            assert unumQ( ibit | (s // 2 + integerPart(f) * ulpu + ut - fsizemax))
            return ibit | (s // 2 + integerPart(f) * ulpu + ut - fsizemax)
        # If the left two exponent bits are 01 or 10,
        # squeeze out the second bit; if that leaves a subnormal exponent,
        # shift the hidden bit and fraction bits right
        if left2 <= 2:
            e = int(((expomask(u) - 3 * mask) & u) + (u & (2*mask)) / 2)
            if e == 0:
                f = Fraction(2**fsize(u) + f, 2)
            ibit = ubitmask if fractionalPart(f) > 0 else 0
            assert unumQ(ibit | (int(s / 2) + e + integerPart(f) * ulpu + ut - fsizemax))
            return ibit | (int(s / 2) + e + integerPart(f) * ulpu + ut - fsizemax)
        # If the first two exponent bits are 11,
        # always get an unbounded unum, all 1s for fraction:
        assert unumQ( int(((u & signmask(u)) + (fm - signmask(u))) / 2) | ut - fsizemax)
        return int(((u & signmask(u)) + (fm - signmask(u))) / 2) | ut - fsizemax
    raise TypeError(u)

def timesposleft(x, y):
    xb = x[1]
    yb = y[1]
    x = x[0]
    y = y[0]
    if (x, xb) == (0, closed):
        if (y, yb) == (float('inf'), closed):
            return (float('nan'), open)
        else:
            return (0, closed)
    if (y, yb) == (0, closed):
        if (x, xb) == (float('inf'), closed):
            return (float('nan'), open)
        else:
            return (0, closed)
    if (x, xb) == (0, open):
        if (y, yb) == (float('inf'), closed):
            return (float('inf'), closed)
        else:
            return (0, open)
    if (y, yb) == (0, open):
        if (x, xb) == (float('inf'), closed):
            return (float('inf'), closed)
        else:
            return (0, open)
    if (x, xb) == (float('inf'), closed) or (y, yb) == (float('inf'), closed):
        return (float('inf'), closed)
    return (x*y, xb or yb)

def timesposright(x, y):
    xb = x[1]
    yb = y[1]
    x = x[0]
    y = y[0]
    if (x, xb) == (float('inf'), closed):
        if (y, yb) == (0, closed):
            return (float('nan'), open)
        else:
            return (float('inf'), closed)
    if (y, yb) == (float('inf'), closed):
        if (x, xb) == (0, closed):
            return (float('nan'), open)
        else:
            return (float('inf'), closed)
    if (x, xb) == (float('inf'), open):
        if (y, yb) == (0, closed):
            return (0, closed)
        else:
            return (float('inf'), open)
    if (y, yb) == (float('inf'), open):
        if (x, xb) == (0, closed):
            return (0, closed)
        else:
            return (float('inf'), open)
    if (x, xb) == (0, closed) or (y, yb) == (0, closed):
        return (0, closed)
    return (x*y, xb or yb)

def neg(x):
    return (-x[0], x[1])

def unionfix(end1, end2):
    return sorted(set(list(end1) + list(end2)))

def timesg(x, y):
    xlo, xhi = x[0]
    xlob, xhib = x[1]
    ylo, yhi = y[0]
    ylob, yhib = y[1]
    lcan = []
    rcan = []
    # If any value is NaN, the result is also NaN.
    if math.isnan(xlo) or math.isnan(xhi) or math.isnan(ylo) or math.isnan(yhi):
        return ((float('nan'), float('nan')), (open, open))
    # Lower left corner is in upper right quadrant, facing uphill:
    if xlo >= 0 and ylo >= 0:
        lcan = unionfix(lcan, (timesposleft((xlo, xlob), (ylo, ylob)),))
    # Upper right corner is in lower left quadrant, facing uphill:
    if (xhi < 0 or (xhi == 0 and xhib)) and (yhi < 0 or (yhi == 0 and yhib)):
        lcan = unionfix(lcan, (timesposleft((-xhi, xhib), (-yhi, yhib)),))
    # Upper left corner is in upper left quadrant, facing uphill:
    if (xlo < 0 or (xlo == 0 and not xlob)) and (yhi > 0 or (yhi == 0 and not yhib)):
        lcan = unionfix(lcan, (neg(timesposright((-xlo, xlob), (yhi, yhib))), ))
    # Lower right corner is in lower right quadrant, facing uphill:
    if (xhi > 0 or (xhi == 0 and not xhib)) and (ylo < 0 or (ylo == 0 and not ylob)):
        lcan = unionfix(lcan, (neg(timesposright((xhi, xhib), (-ylo, ylob))),))
    # Upper right corner is in upper right quadrant, facing downhill:
    if (xhi > 0 or (xhi == 0 and not xhib)) and (yhi > 0 or (yhi == 0 and not yhib)):
        rcan = unionfix(rcan, (timesposright((xhi, xhib), (yhi, yhib)),))
    # Lower left corner is in lower left quadrant, facing downhill:
    if (xlo < 0 or (xlo == 0 and not xlob)) and (ylo < 0 or (ylo == 0 and not ylob)):
        rcan = unionfix(rcan, (timesposright((-xlo, xlob), (-ylo, ylob)),))
    # Lower right corner is in upper left quadrant, facing downhill:
    if (xhi < 0 or (xhi == 0 and xhib)) and ylo >= 0:
        rcan = unionfix(rcan, (neg(timesposright((-xhi, xhib), (ylo, ylob))),))
    # Upper left corner is in lower right quadrant, facing downhill:
    if xlo >= 0 and (yhi < 0 or (yhi == 0 and yhib)):
        rcan = unionfix(rcan, (neg(timesposright((xlo, xlob), (-yhi, yhib))),))

    if any(isinstance(can, float) and math.isnan(can) for can in flatten(lcan)) or \
       any(isinstance(can, float) and math.isnan(can) for can in flatten(rcan)):
           (timesleft, timesright) = (float("nan"), float("nan"))
           (openleft, openright) = (open, open)
    (timesleft, timesright) = (lcan[0][0], rcan[-1][0])
    (openleft, openright) = (lcan[0][1], rcan[-1][1])
    if len(lcan) > 1:
        if lcan[0][0] == lcan[1][0] and (not lcan[0][1] or not lcan[1][1]):
            openleft = closed
    if len(rcan) > 1:
        if rcan[-1][0] == rcan[-2][0] and (not rcan[-1][1] or not rcan[-2][1]):
            openright = closed
    return ((timesleft, timesright), (openleft, openright))

def timesu(u, v):
    if uQ(u) and uQ(v):
        w = g2u(timesg(u2g(u), u2g(v)))
        global ubitsmoved, numbersmoved
        ubitsmoved += nbits(u) + nbits(v) + nbits(w)
        numbersmoved += 3
        return w

def divideposleft(x, y):
    xb = x[1]
    yb = y[1]
    x = x[0]
    y = y[0]
    if (y, yb) == (0, closed):
        return (float('nan'), open)

    if (x, xb) == (float('inf'), closed):
        if (y, yb) == (float('inf'), closed):
            return (float('nan'), open)
        else:
            return (float('inf'), closed)
    if (x, xb) == (0, closed) or (y, yb) == (float('inf'), closed):
        return (0, closed)
    if (x, xb) == (float('inf'), open) or (y, yb) == (0, open):
            return (float('inf'), open)
    return (x/y, xb or yb)


def divideposright(x, y):
    xb = x[1]
    yb = y[1]
    x = x[0]
    y = y[0]
    if (y, yb) == (0, closed):
        return (float('nan'), open)

    if (x, xb) == (float('inf'), closed):
        if (y, yb) == (float('inf'), closed):
            return (float('nan'), open)
        else:
            return (float('inf'), closed)
    if (x, xb) == (0, closed) or (y, yb) == (float('inf'), closed):
        return (0, closed)
    if (x, xb) == (float('inf'), open) or (y, yb) == (0, open):
            return (float('inf'), open)
    return (x/y, xb or yb)



def divideg(x, y):
    xlo, xhi = x[0]
    xlob, xhib = x[1]
    ylo, yhi = y[0]
    ylob, yhib = y[1]
    lcan = []
    rcan = []
    # If any value is NaN, the result is also NaN.
    if math.isnan(xlo) or math.isnan(xhi) or math.isnan(ylo) or math.isnan(yhi) \
            or ((ylo < 0 or (ylo == 0 and not ylob)) and (yhi > 0 or (yhi == 0 and not yhib))):
        return ((float('nan'), float('nan')), (open, open))
    # Upper left corner is in upper right quadrant, facing uphill:
    if xlo >= 0 and (yhi > 0 or (yhi == 0 and not yhib)):
        lcan = unionfix(lcan, (divideposleft((xlo, xlob), (yhi, yhib)),))
    # Lower right corner is in lower left quadrant, facing uphill:
    if (xhi < 0 or (xhi == 0 and xhib)) and (ylo < 0 or (ylo == 0 and not ylob)):
        lcan = unionfix(lcan, (divideposleft((-xhi, xhib), (-ylo, ylob)),))
    # Lower left corner is in upper left quadrant, facing uphill:
    if (xlo < 0 or (xlo == 0 and not xlob)) and ylo >= 0:
        lcan = unionfix(lcan, (neg(divideposright((-xlo, xlob), (ylo, ylob))),))
    # Upper right corner is in lower right quadrant, facing uphill:
    if (xhi > 0 or (xhi == 0 and not xhib)) and (yhi < 0 or (yhi == 0 and yhib)):
        lcan = unionfix(lcan, (neg(divideposright((xhi, xhib), (-yhi, yhib))),))
    
    # Lower right corner is in upper right quadrant, facing downhill:
    if (xhi > 0 or (xhi == 0 and not xhib)) and ylo >= 0:
        rcan = unionfix(rcan, (divideposright((xhi, xhib), (ylo, ylob)),))
    # Upper left corner is in lower left quadrant, facing downhill:
    if (xlo < 0 or (xlo == 0 and not xlob)) and (yhi < 0 or (yhi == 0 and yhib)):
        rcan = unionfix(rcan, (divideposright((-xlo, xlob), (-yhi, yhib)),))
    # Upper right corner is in upper left quadrant, facing downhill:
    if (xhi < 0 or (xhi == 0 and xhib)) and (yhi > 0 or (yhi == 0 and not yhib)):
        rcan = unionfix(rcan, (neg(divideposleft((-xhi, xhib), (yhi, yhib))),))
    # Lower left corner is in lower right quadrant, facing downhill:
    if xlo >= 0 and (ylo < 0 or (ylo == 0 and not ylob)):
        rcan = unionfix(rcan, (neg(divideposleft((xlo, xlob), (-ylo, ylob))),))
 
    if any(isinstance(can, float) and math.isnan(can) for can in flatten(lcan)) or \
       any(isinstance(can, float) and math.isnan(can) for can in flatten(rcan)):
           (divleft, divright) = (float("nan"), float("nan"))
           (openleft, openright) = (open, open)
    (divleft, divright) = (lcan[0][0], rcan[-1][0])
    (openleft, openright) = (lcan[0][1], rcan[-1][1])
    if len(lcan) > 1:
        if lcan[0][0] == lcan[1][0] and (not lcan[0][1] or not lcan[1][1]):
            openleft = closed
    if len(rcan) > 1:
        if rcan[-1][0] == rcan[-2][0] and (not rcan[-1][1] or not rcan[-2][1]):
            openright = closed
    return ((divleft, divright), (openleft, openright))

def divideu(u, v):
    if uQ(u) and uQ(v):
        w = g2u(divideg(u2g(u), u2g(v)))
        global ubitsmoved, numbersmoved
        ubitsmoved += nbits(u) + nbits(v) + nbits(w)
        numbersmoved += 3
        return w

def squareg(g):
    if gQ(g):
        g1 = g[0][0]
        g2 = g[0][1]
        b1 = g[1][0]
        b2 = g[1][1]
        #XXX: what is indeterminate?
        if math.isnan(g1) or math.isnan(g2):
            return f2g(float('nan'))
        t1 = g1*g1
        t2 = g2*g2
        tset = tuple(sorted(((t1, b1), (t2, b2))))
        # See if 0 is in the range
        if (g1 < 0 and g2 > 0) or (g1 > 0 and g2 < 0) or (g1 == 0 and not b1) \
                or (g2 == 0 and not b2):
            if t1 == t2:
                return ((0, t1), (closed, b1 and b2))
            else:
                return ((0, tset[1][0]), (closed, test[1][1]))
        return transpose(tset)
# Square in the u-layer, with tallying of bits and number.
def squareu(u):
    if uQ(u):
        i = nbits(u)
        v = g2u(squareg(u2g(u)))
        global ubitsmoved, numbersmoved
        ubitsmoved += nbits(u) + nbits(v)
        numbersmoved += 2
        return v

# We don't have a way to express the exact
# square root of a Fraction so we approximate
# it up the required precision instead.
# This just uses newton's method for approximation.
def sqrtFraction(x):
    x0 = x
    while True:
        adjust = Fraction(x0*x0 - x, 2*x0)
        x1 = x0 - adjust
        if adjust < Fraction(1,2):
            break
        x0 = x1
    # check if x is a square number
    x0int = x0.limit_denominator(1)
    if x0int*x0int == x:
        return x0int
    # otherwise keep improving our approximation until we run out of precision
    while True:
        adjust = Fraction(x0*x0 - x, 2*x0)
        x1 = x0 - adjust
        if x2u(x0) == x2u(x1):
            break
        x0 = x1
    return x1

# BUG: the original seems to have an unnecessary i = nbits[u] at the begining
def sqrtu(u):
    if uQ(u):
        g = u2g(u)
        if math.isnan(g[0][0]) or math.isnan(g[0][1]) or g[0][0] < 0:
            result = f2g(NaN)
        else:
            result = ((sqrtFraction(g[0][0]), sqrtFraction(g[0][1])), (g[1][0], g[1][1]))
        v = g2u(result)
        global ubitsmoved, numbersmoved
        ubitsmoved += nbits(u) + nbits(v)
        numbersmoved += 2
        return v

# The "left" power function tables for general intervals.
def powposleft(x, y):
    (x, xb) = x
    (y, yb) = y
    if x >= 1 and y >= 0:
        if (x, xb) == (1, closed):
            if (y, yb) == (float('inf'), closed):
                return (NaN, open)
            else:
                return (1, closed)
        elif (y, yb) == (0, closed):
            if (x, xb) == (float('inf'), closed):
                return (NaN, open)
            else:
                return (1, closed)
        elif (x, xb) == (1, open):
            if (y, yb) == (float('inf'), closed):
                return (float('inf'), closed)
            else:
                return (1, open)
        elif (y, yb) == (0, open):
            if (x, xb) == (float('inf'), closed):
                return (float('inf'), closed)
            else:
                return (1, open)
        elif (x, xb) == (float('inf'), closed) or y(y, yb) == (float('inf'), closed):
            return ('float', closed)
        else:
            return (x**y, xb or yb)

# XXX: BUG the comment is wrong (says left) in upstream
# The "right" power function tables for general intervals.
def powposright(x, y):
    (x, xb) = x
    (y, yb) = y
    if x >= 1 and y >= 0:
        if (x, xb) == (float('inf'), closed):
            if (y, yb) == (0, closed):
                return (NaN, open)
            else:
                return (float('inf'), closed)
        elif (y, yb) == (float('inf'), closed):
            if (x, xb) == (1, closed):
                return (NaN, open)
            else:
                return (float('inf'), closed)
        elif (x, xb) == (float('inf'), open):
            if (y, yb) == (0, closed):
                return (1, closed)
            else:
                return (float('inf'), open)
        elif (y, yb) == (float('inf'), open):
            if (x, xb) == (1, closed):
                return (1, closed)
            else:
                return (float('inf'), open)
        elif (x, xb) == (1, closed) or y(y, yb) == (0, closed):
            return (1, closed)
        else:
            return (x**y, xb or yb)

# Reciprocal helper function for the power function.
def rec(x):
    (x, xb) = x
    if math.isnan(x):
        return NaN
    elif x == 0:
        return float('inf')
    else:
        return 1/x

def powg(x, y):
    ((xlo, xhi), (xlob, xhib)) = x
    ((ylo, yhi), (ylob, yhib)) = y
    NaNg = ((NaN, NaN), (open, open))
    lcan = ()
    rcan = ()
    # If any value is NaN, the result is also NaN.
    if math.isnan(xlo) or math.isnan(xhi) or math.isnan(yhi) or math.inan(yhi):
        return NaNg
    # Do not allow exact zero to a negative or zero power,
    # unless the negative power is in an exact even integer.

# Fused multiply-add in the g-layer.
def fmag(ag, bg, cg):
    if gQ(ag) and gQ(bg) and gQ(cg):
        return plusg(timesg(ag, bg), cg)

def famg(ag, bg, cg):
    if gQ(ag) and gQ(bg) and gQ(cg):
        return timesg(plusg(ag, bg), cg)

def splitub(ub):
    ((g1, g2), (b1, b2)) = u2g(ub)
    if math.isnan(g1) and math.isnan(g2) and b1 and b2:
        return ub
    elif not b1 and not b2:
        if g1 == g2:
            # Cannot split exact single values
            return [ub]
        else:
            # else cleave off exact endpoints
            return [(ub[0],), g2u(((g1, g2), (open, open))), (ub[-1],)]
    elif b1 and not b2:
        # cleave off exact right endpoint
        return [g2u(((g1, g2), (open, open))), (ub[-1],)]
    elif not b1 and b2:
        # cleave off exact left endpoint
        return [(ub[1],), g2u(((g1, g2), (open, open)))]
    if g1 == float('-inf'):
        if g2 == -maxreal:
            # Cannot split the negative "many" region
            return [ub]
        else:
            return [(negbigu + ubitmask,), (negbigu,), g2u(((-maxreal, g2), (open, open)))]
    elif g2 == float('inf'):
        if g1 == maxreal:
            # Cannot split the positive "many" region
            return [ub]
        else:
            return [g2u(((g1, maxreal), (open, open))), (maxrealu,), (maxrealu + ubitmask)]
    else:
        assert not isinstance(g1, float)
        assert not isinstance(g2, float)
        # See if open interval contains a unum different from either endpoint:
        gm = u2g(x2u(Fraction((g1+g2),2)))
        assert not isinstance(gm[0][0], float)
        assert not isinstance(gm[0][1], float)
        if gm[0][0] > g1:
            return [g2u(((g1, gm[0][0]), (open, open))), (x2u(gm[0][0]),), g2u(((gm[0][0], g2), (open, open)))]
        if gm[0][1] < g2:
            return [g2u(((g1, gm[0][1]), (open, open))), (x2u(gm[0][1]),), g2u(((gm[0][1], g2), (open, open)))]
        # Cannot split; must be the samllest ULP size.
        return [ub]

# Check if an argument is a list of general intervals.
def glistNaNQ(u):
    if isinstance(u, tuple) or isinstance(u, list):
        return reduce(lambda x, y: x or y, [math.isnan(x[0][0]) or math.isnan(x[0][1]) for x in u])

# Check if an argument is a list of general intervals.
def glistQ(u):
    if isinstance(u, tuple) or isinstance(u, list):
        return reduce(lambda x, y: x and y, map(gQ, u))

# Check if an argument is a list of unums or ubounds
def ulistQ(u):
    if isinstance(u, tuple) or isinstance(u, list):
        return reduce(lambda x, y: x and y, map(uQ, u))


# Polynomial helper function that uses powers of x - x0 instead of x.
def polyTg(coeffsg, xg, x0g):
    k = len(coeffsg)
    coeffstg = list(coeffsg)
    if x0g[0][0] == float('-inf') or x0g[0][1] == float('inf'):
        return ((float('-inf'), float('inf')), (closed, closed))
    j = 0
    while j <= k - 1:
        bi = binomial(k - 1, j)
        pg = timesg(coeffsg[-1], ((bi, bi), (closed, closed)))
        i = k - 2
        while i >= j:
            bi = binomial(i, j)
            pg = plusg(timesg(coeffsg[i], ((bi, bi), (closed, closed))), timesg(x0g, pg))
            i -= 1
        coeffstg[j] = pg
        j += 1
    xmg = minusg(xg, x0g)
    pg = coeffstg[k-1]
    i = k - 1
    while i >= 1:
        pg = plusg(timesg(pg, xmg), coeffstg[i-1])
        i -= 1
    return pg

# Polynomial helper function that evaluates a polynomial at the endpoints
# of an inexact unum, and intersects them to tighten the result.
def polyinexactg(coeffsg, xg):
    return intersectg(polyTg(coeffsg, xg, ((xg[0][0], xg[0][0]), (closed, closed))),
            polyTg(coeffsg, xg, ((xg[0][1], xg[0][1]), (closed, closed))))

# Polynomial evaluation of an exact general interval using Horner's rule
def polyexactg(coeffsg, xg):
    assert gQ(xg)
    k = len(coeffsg)
    pg = coeffsg[-1]
    i = k - 2
    while i >= 0:
        pg = plusg(coeffsg[i], timesg(pg, xg))
        assert gQ(pg)
        i -= 1
    return pg

# Bisect an inexact general interval along a coarsest-possible ULP boundary.
def bisect(g):
    #print 'bisect', g
    if gQ(g):
        gL = g[0][0]
        gR = g[0][1]
        if gL < 0 and gR > 0:
            gM = 0
        elif gL == float('-inf') and gR > -maxreal:
            gM = -maxreal
        elif gL < maxreal and gR == float('inf'):
            gM = maxreal
        else:
            #print gL, gR, -maxreal
            m = Fraction(2,1)**floor(log(gR-gL, 2))
            if not math.isinf(m) and Fraction(gL, m).denominator == 1:
                if gR - gL == m:
                    gM = Fraction(gL + gR, 2)
                else:
                    gM = m * floor(gL/m + 1)
            else:
                gM = m * ceil(gL / m)
        # XXX: BUG should this end with g[1][1] or g[1][0] like the upstream?
        return (((gL, gM), (g[1][0], open)), ((gM, gR), (open, g[1][1])))


def tripleEq(x, y):
    if gQ(x) and gQ(y):
        if math.isnan(x[0][0]):
            assert math.isnan(x[0][1])
        if math.isnan(y[0][0]):
            assert math.isnan(y[0][1])
        return x == y or (math.isnan(x[0][0]) and math.isnan(y[0][0]))

# Polynomial evaluation of a general interval without u-layer information loss.
def polyg(coeffsg, xg):
    if glistQ(coeffsg) and gQ(xg):
        k = len(coeffsg)
        if math.isnan(xg[0][0]) or math.isnan(xg[0][1]) or glistNaNQ(coeffsg):
            return ((NaN, NaN), (open, open))
        # Constant case. Just return the first (and only coefficient).
        if k == 1:
            return coeffsg[0]
        # Linear case is a fused multiply-add; no dependency problem.
        if k == 2:
            return fmag(coeffsg[1], xg, coeffsg[1])
        # Exact argument is also easy, since no dependency problem.
        if xg[0][0] == xg[0][1]:
            return polyexactg(coeffsg, xg)
        # Quadratic or higher requires finesse. Intersect the two
        # endpoint-based evaluations.
        trials = (xg,)
        gL = polyexactg(coeffsg, ((xg[0][0], xg[0][0]), (closed, closed)))
        if xg[1][0]:
            gL = (gL[0], (open, open))
        gR = polyexactg(coeffsg, ((xg[0][1], xg[0][1]), (closed, closed)))
        if xg[1][1]:
            gR = (gR[0], (open, open))
        if gL[0][0] < gR[0][0] or (gL[0][0] == gR[0][0] and not gL[1][0]):
            (min, minQ) = transpose(gL)[0]
        else:
            (min, minQ) = transpose(gR)[0]
        if gL[0][1] > gR[0][1] or (gL[0][1] == gR[0][1] and not gL[1][1]):
            (max, maxQ) = transpose(gL)[1]
        else:
            (max, maxQ) = transpose(gR)[1]
        #assert gQ(((min, max), (minQ, maxQ)))
        while len(trials) >= 1:
            # print 'trials', trials
            pg = polyinexactg(coeffsg, trials[0])
            #
            assert gQ(pg)
            #assert gQ(((min, max), (minQ, maxQ)))
            if tripleEq(intersectg(u2g(g2u(pg)), u2g(g2u(((min, max), (minQ, maxQ))))), u2g(g2u(pg))):
                trials = trials[1:]
            else:
                trials = bisect(trials[0]) + trials[1:]
                gM = polyexactg(coeffsg, ((trials[0][0][1], trials[0][0][1]), (closed, closed)))
                #print 'gM', gM, min, max
                if gM[0][0] < min or gM[0][0] == min and not gM[1][0]:
                    (min, minQ) = transpose(gM)[0]
                if gM[0][1] > max or gM[0][1] == max and not gM[1][1]:
                    (max, maxQ) = transpose(gM)[1]
                #print 'prefin', ((min, max), (minQ, maxQ))
                ((min, max), (minQ, maxQ)) = u2g(g2u(((min, max), (minQ, maxQ))))
                assert gQ(((min, max), (minQ, maxQ)))
                #print 'fin', ((min, max), (minQ, maxQ))
        return ((min, max), (minQ, maxQ))

def polyu(coeffsu, u):
    if ulistQ(coeffsu) and uQ(u):
        coeffsg = [u2g(coeff) for coeff in coeffsu]
        return g2u(polyg(coeffsg, u2g(u)))

# Look for alternative unum string that favors the exponent.
def favore(unum):
    if unumQ(unum):
        u = demotef(promotee(unum))
        # BUG: there are two inexQ(u) on this line
        if inexQ(u) or esize(u) == esizemax or fsize(u) == 1 or inexQ(u):
            return u
        while fsize(u) > 1 and exQ(demotef(u)):
            u = demotef(u)
        return u

# Look for alternative unum string that favors the fraction.
def favorf(unum):
    if unumQ(unum):
        u = demotee(padu(unum))
        if inexQ(u):
            return u
        else:
            while fsize(u) > 1 and exQ(demotef(u)):
                u = demotef(u)
            return u

# Find the right neighbor of a unum
def nborhi(u, minpower):
    if unumQ(u):
        ut = x2u(0) if u == ((utagmask + signmask(u)) & u) else u
        s = -1**(sign(ut))
        overflow = False
        if minpower < log(smallsubnormal, 2):
            ulpminu = smallsubnormalu
        elif minpower > log(maxreal, 2):
            overflow = True
            ulpminu = x2u(2**floor(log(maxreal, 2)))
        else:
            ulpminu = x2u(2**minpower)
        ulpmin = u2g(ulpminu)[0][0]
        if u == posinfu or u == sNaNu or u == qNaNu:
            return qNaNu
        elif u == neginfu:
            if overflow:
                if utagmask == 1:
                    return x2u(-2) + ubitmask # Warlpiri enviroment
                else:
                    return x2u(-3) + ubitmask
            else:
                negbigu + ubitmask
                # if -inf, use the (-inf, x) unum with the most x,
                # unless the requested minpower caused overflow.
        elif inexQ(u):
            # If inexact always use the exact upper value
            return x2u(u2g(u)[0][1])
        elif overflow and u == x2u(2) and utagmask == 1:
            return x2u(2) + ubitmask # Warlpiri
        elif overflow and u == x2u(3) and utagmask != 1:
            return x2u(3) + ubitmask
        else:
            # Reduce ULP until it equals ulpmin,
            # or we run out of exponent and fraction bits
            t = u2g(ut)[0][0]
            ut = x2u(t)
            while not Fraction(t,ulpmin).denominator == 1:
                ulpmin //= 2
            while ulphi(ut) < ulpmin and ut != favorf(ut):
                ut = favorf(ut)
            while esize(ut) < esizemax and ulphi(promotee(ut)) >= ulpmin:
                ut = promotee(ut)
            while fsize(ut) < fsizemax and ulphi(promotee(ut)) >= ulpmin:
                ut = promotef(ut)
            return ut + s * ubitmask

# Find the left neighbor of a unum
def nborlo(u, minpower):
    if unumQ(u):
        # Watch out for negative zero, otherwise use nborhi
        if sameuQ((x2u(u),), (x2u(0),)) and minpower < log(smallsubnormal, 2):
            return smallsubnormalu + signbigu - ubitmask
        else:
            return negateu((nborhi(negateu((u,)), minpower),)[0])

def drop(l, which):
    return l[0:which[0]] + l[which[1]+1:]
# The ubinsert insertion routine assumes ubset is a sorted
# set of disjoint ubounds, and has no NaN ubounds.
def ubinsert(ubset, ub):
    if uboundQ(ub):
        k = len(ubset)
        lefttouch = False
        righttouch = False
        newset = ubset
        if k == 0:
            # print "   First element. ", view(ub)
            return (ub,)
        j = k
        print "searching for matches"
        # XXX: this could be written in a much more pythonic style
        while j > 0 and gtuQ(newset[j-1], ub):
            j -= 1
        if j > 0:
            lefttouch = nnequQ((nborhi(newset[j-1][-1], float('-inf')),), ub)
        if j < k:
            righttouch = nnequQ((nborlo(newset[j][0], float('-inf')),), ub)
        if lefttouch and righttouch:
            newset = drop(newset, (j-1, k-1)) + [(newset[j-1][0], newset[j][-1])] + drop(newset, (0, j))
            print "   Joined on both sides. "
        elif lefttouch and not righttouch:
            newset = drop(newset, (j-1, k-1)) + [(newset[j+1][0], ub[-1])] + drop(newset, (0, j-1))
            print "   Joined on left side. "
        elif not lefttouch and righttouch:
            newset = drop(newset, (j, k-1)) + [(ub[-1], newset[j][-1])] + drop(newset, (0, j))
            print "   Joined on right side. "
        else:
            print "   Inserted new ubound, not touching"
            if j + 1 > k:
                newset = newset + (ub,) + drop(newset, (0, j-1))
            else:
                newset = drop(newset, (j,k-1)) + (ub,) + drop(newset, (0, j-1))
        return newset

# The try-everything solver. Has commented-out print statements
# that can be uncommented to see why long-running solvers are taken
# too long.
def solveforub(domain, conditionQ):
    sols = []
    trials = domain
    print "COND-d", domain
    while len(trials) > 0:
        new = []
        for ub in trials:
            b = conditionQ(ub)
            print "COND", b, view(ub), ub
            if b:
                temp = splitub(ub)
                if len(temp) == 1:
                    # unsplittable. Join to existing region or start new one.
                    sols = ubinsert(sols, temp[0])
                    print "joined", view(ub), "to sols: ", "".join([view(sol) for sol in sols])
                else:
                    new = temp + new
        trials = new
        print "The 'new' list:", [view(i) for i in new]
    return sols

ubitsmoved = numbersmoved = 0



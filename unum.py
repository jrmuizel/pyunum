from __future__ import division

import math
open = True
closed = False
print open
from fractions import *
e = 3
f = 4
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
qNaNu = posinfu + ubitmask
sNaNu = neginfu + ubitmask
NaN = float("nan")

maxreal = 2**(2**(esizemax-1)) * (2**fsizemax - 1) / (2**fsizemax)
smallsubnormal = 2**(2 - 2**(esizemax-1) - fsizemax)

def log(x, base):
    return math.log(x, base)

def floor(x):
    return int(math.floor(x))

def unumQ(x):
    return True

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
    return True

def exQ(x):
    return True

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
        return ((-1)**sign(u))*(2**expovalue(u))*(hidden(u) + Fraction(frac(u), 2**fsize(u)))

def bigu(u):
    if unumQ(u):
        return expomask(u) + fracmask(u) + (efsizemask & u) - ulpu * ((u & efsizemask) == efsizemask)

def gQ(x):
    if isinstance(x, tuple):
        if len(x) == 2:
            if isinstance(x[0], tuple):
                if len(x[0]) == 2:
                    if isinstance(x[1], tuple):
                        if len(x[1]) == 2:
                            if isinstance(x[1][0], bool) and isinstance(x[1][1], bool) and floatQ(x[0][0]) and floatQ(x[0][1]):
                                if math.isnan(x[1][1]) or math.isnan(x[0][1]):
                                    return True
                                if (x[0][0] == x[0][1] and not x[1][0] and not x[1][1]) or x[0][0] < x[0][1]:
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
            return 1 + numbits(u)
        else:
            return 1 + numbits(ub[0]) + numbits(ub[1])

def u2g(u):
    if uQ(u):
        if unumQ(u):
            return unum2g(u)
        else:
            return ubound2g(u)

def x2u(x):
    if floatQ(x):
        if math.isnan(x):
            return qNaNu
        elif x == float("inf"):
            return posinfu
        elif x == float("-inf"):
            return neginfu
        elif abs(x) > maxreal:
            return maxrealu + ubitbask + (signbigu if x < 0 else 0)
        elif x == 0:
            return 0
        elif abs(x) < smallsubnormal:
            return utagmask + (signbigu if x < 0 else 0)
        elif abs(x) < u2f(smallnormalu):
            y = abs(x) / smallsubnormal
            y = (signbigu if x < 0 else 0) + efsizemask + (ubitmask if y != math.floor(y) else 0) + (math.floor(y), utagsize)
            # XXX is this right?
            assert(False)
            while ((3 << (utagsize - 1)) & y) == 0:
                y = (y - (efsizemask & y))/2 + (efsizemask & y) - 1
        else:
            n = 0
            y = abs(x) / 2**(scale(x))
            n = 0
            while math.floor(y) != y and n < fsizemax:
                n += 1
                y *= 2
            if y == math.floor(y): # the value is representable
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
                    y1 += signmask(y)
                return y1

def plusg(x, y):
    ((xlo, xhi), (xlob, xhib)) = x
    ((ylo, yhi), (ylob, yhib)) = y
    if math.isnan(xlo) or math.isnan(xhi) or math.isnan(ylo) or math.isnan(yhi):
        return ((NaN, open), (NaN,open))
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

def ubleft(xleft):
    open = xleft[1]
    u = x2u(xleft[0])
    x = xleft[0]
    if x == float('-inf'):
        if open:
            return negopeninfu
        else:
            return neginfu
    if x == 0 and open:
        return negopenzerou
    if u2f(u) == x:
        if open:
            return (u - (ulpu * (x < 0))) | ubitmask
        else:
            return u
    else:
        return u | (open * ubitmask)


def g2u(g):
    if gQ(g):
        ulo = g[0][0]
        uhi = g[0][1]
        blo = g[1][1]
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

ubitsmoved = numbersmoved = 0

five = x2u(34.2)
def print_unum(u):
    print sign(u), expovalue(u), hidden(u), Fraction(frac(u), 2**fsize(u))
print_unum(five)
print frac(five), fsize(five)
print "numbits", numbits(five), "of", maxubits, fsize(five), "of", fsizemax, esize(five), "of", esizemax
print u2f(five)

print plusg(u2g(x2u(34.2)), u2g(x2u(0)))
print u2f(plusu(x2u(34.2), x2u(5)))


"""
def f2g("""

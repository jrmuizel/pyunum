import sys
import unum_config
unum_config.e = 0
unum_config.f = 0
from unum import *
#import unum
"""
for i in range(0, 16):
    assert view(i) == view(g2u(u2g(i)))

for i in range(0,16):
    print bin(i), view(i)

for i in range(0,16):
    for j in range(0, 16):
        print view(i) + " * " + view(j) + " = " + view(timesu(i, j))

for i in range(0,16):
    for j in range(0, 16):
        print view(i) + " + " + view(j) + " = " + view(timesu(i, j))

for i in range(0,16):
    for j in range(0, 16):
        print view(i) + " / " + view(j) + " = " + view(divideu(i, j))
walpri = {}
for i in range(0, 16):
    for j in range(0, 16):
        if uQ((i, j)):
            walpri[view((i,j))] = (i,j)
print len(walpri)
for w in sorted(walpri):
    print w
for i in walpri.values():
    for j in walpri.values():
        print view(i) + " + " + view(j) + " = " + view(plusu(i, j))
sys.exit(1)
"""
unum_config.e = 3
unum_config.f = 4
import unum
reload(unum)
from unum import *
print e
print u2g(5)
print view(x2u(5))
five = x2u(34.2)
print view(five)
def print_unum(u):
    print sign(u), expovalue(u), hidden(u), Fraction(frac(u), 2**fsize(u))
print_unum(five)
print frac(five), fsize(five)
print "numbits", numbits(five), "of", maxubits, fsize(five), "of", fsizemax, esize(five), "of", esizemax
print u2f(five)
print x2u(5)
print u2g(x2u(9005))
print plusg(u2g(x2u(34.2)), u2g(x2u(0)))
print 'times', u2f(x2u(34.2)), u2f(x2u(1))
print 'times', timesg(u2g(x2u(34.2)), u2g(x2u(1)))
print view(plusu(x2u(34.2), x2u(0)))
print view(x2u(34.2)), view(timesu(x2u(34.2), x2u(2)))
print view(x2u(34.2)), view(timesu(divideu(x2u(34.2), x2u(2)), x2u(2)))
print view(squareu(x2u(-5)))
#print view(polyu((x2u(5), x2u(6), x2u(1)), x2u(34.2)))
print 'final'
#print view(polyu((x2u(-5), x2u(6), x2u(1)), g2u(((float('-inf'), float('inf')), (True, True)))))
assert promotee(16) == 32
#t = polyu((x2u(-5), x2u(6), x2u(1)), (4294967038, 896))
#ref = g2u(((-14, float('inf')), (False, True)))
#assert sameuQ(t, ref)
#print view(t), view(ref)
print view((4294967038, 896))
#print neginfu, posinfu
x = solveforub(((neginfu, posinfu),), lambda ub: nnequQ(polyu((x2u(-5), x2u(6), x2u(1)), ub), x2u(0)))
#x = solveforub(((x2u(-10), x2u(10)),), lambda ub: nnequQ(polyu((x2u(-5), x2u(6), x2u(1)), ub), x2u(0)))
print [view(i) for i in x]

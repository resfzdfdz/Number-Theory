import random
import math
from ENTF import *

###   Fast Number Theory Transformation   ###

###   Verify NTT parameters   ###
def verify_NTTP(NTT_Param):
    (a, N, M) = NTT_Param
    for i in range(1, N):
        if (pow(a, i, M) == 1):
            print('NTT Parameter fail!\n')
            return 0
    if (pow(a, N, M) == 1):
        print('NTT Parameter pass!\n')
        return 1
    else:
        print('NTT Parameter fail!\n')
        return 0

###   random polynomial coefficient generation   ###
def vector_generator(n, coef = 2 ** 16):
    r = []
    for i in range(n):
        r.append(random.randint(0, coef - 1))
    return r

###   adding 0 at the end of the vector   ###
###   we use 2n point NTT to get result   ###
###   of polynomial multiplication        ###
def padding(n):
    r = []
    for i in range(n):
        r.append(0)
    return r

###   transform polynimial coefficient to ###
###   integer in order to do integer muls ###
def vector2int(x, coef = 2 ** 16):
    s = 0
    for i in range(len(x)):
        s += x[i] * (coef ** i)
    return s

###   Fermat numbers, when n = 1, 2, 3, 4 ###
###   Fermat number is prime. However,    ###
###   when n >=5, it is likely to be a    ###
###   composite rather than prime         ###
def Fermat(n):
    return (2 ** (2 ** n) ) + 1

##  d = b / a mod p
###  modular division                     ###
def moddiv1(a, b, p):
    u = a % p
    v = p
    m = b % p
    n = 0

    while (u != 1):
        r = v // u
        u, v = v - r * u, u
        m, n = (n - r * m) % p, m % p
    return m

###   Number Theory Transform (not fast)  ###
def NTT(x, NTT_Param):
    assert verify_NTTP(NTT_Param) == 1
    
    (a, N, M) = NTT_Param
    xk = []
    for k in range(N):
        s = 0
        for n in range(N):
##            s += ( x[n] * (a ** (n * k) ) ) % M
            s += ( x[n] * pow(a, n * k, M) )
        xk.append(s % M)
    return xk

###   Inverse Number Theory Transform     ###
def INTT(xk, NTT_Param):
    assert verify_NTTP(NTT_Param) == 1
    
    (a, N, M) = NTT_Param
    x = []
    for n in range(N):
        s = 0
        for k in range(N):
##            s += ( xk[k] * moddiv(a ** (n * k), 1, M ) ) % M
            s += ( xk[k] * moddiv1(pow(a, n * k, M), 1, M ) )
        Ninv = moddiv1(N, 1, M)
        x.append(s * Ninv % M)
    return x

###   Circular Convolution                ###
###   Note that all element must be       ###
###   integers, and the result is modulo  ###
###   by prime M                          ###
def NTT_conv(x, y, NTT_Param):
    (a, N, M) = NTT_Param
    z = []
    for n in range(N):
        s = 0
        for m in range(N):
            s += x[m] * y[(n - m) % N]
        z.append(s % M)
    return z

###   Scalar multiplication               ###
def scalar_mul(x, y, M):
    assert len(x) == len(y), 'Data length error!'
    z = []
    for i in range(len(x)):
        z.append(x[i] * y[i] % M)
    return z

def vecadd(x, y):
    r = []
    for i in range(x):
        r.append(x[i] + t[i])
    return r

def minus_conv(x, y, NTT_Param):
    (a, N, M) = NTT_Param
    psai = modsqroot(a, M)

    v = [(psai ** i) % q for i in range(512)]

    a1 = scalar_mul(x, v, M)
    b1 = scalar_mul(y, v, M)
    ak = NTT(a1, NTT_Param)
    bk = NTT(b1, NTT_Param)    
    ck = scalar_mul(ak, bk, M)
    c1 = INTT(ck, NTT_Param)

    u = [moddiv1(psai ** i, 1, q) for i in range(512)]
    c = scalar_mul(c1, u, M)

    return c

def verify_minus_conv(x, y, NTT_Param):
    (a, N, M) = NTT_Param

    r = x
    s = y
    
    pad = padding(N)
    r.extend(pad)
    s.extend(pad)
    
    z = NTT_conv(r, s, (a, 2* N, M) )

    w = [ (z[i] - z[i + 512]) % M for i in range(512)]

    return w

def verify_NTT_Trans():
    r = 1024
    seg = 16
    n = r // seg    
    N = 2 * n
    
    while (1):
        M = geneprime(2 ** 39, 2 ** 65)
        if ( (M - 1) % N != 0):
            continue
        else:
            k = (M - 1) // N
            break
        
    g = primitive_root(M)
    a = pow(g, k, M)
    NTT_Param = (a, N, M)

    print('NTT_Param = {st}\n'.format(st = NTT_Param) )

    #####  random vector generation   ###
    x = vector_generator(n, 2 ** seg)
    y = vector_generator(n, 2 ** seg)
    ##
    xi = vector2int(x)
    yi = vector2int(y)
    ##
    pad = padding(n)
    ##
    #####  padding   ###
    x.extend(pad)
    y.extend(pad)

    #  Test NTT   ###
    u = NTT(x, NTT_Param)
    v = INTT(u, NTT_Param)
    assert v == x, 'NTT error!'

    ###   NTT to polynomial multiplication   ###
    ref_z = xi * yi

    xk = NTT(x, NTT_Param)
    yk = NTT(y, NTT_Param)

    zk = scalar_mul(xk, yk, M)
    z = INTT(zk, NTT_Param)

    zi = vector2int(z)

    print(hex(ref_z) + '\n')
    print(hex(zi) + '\n')

    if (ref_z == zi):
        print('NTT Transfrom Success!\n')
    else:
        print('NTT Transfrom Fail!\n')    


def redunt2int(x):
    x_int = vector2int(x)
    p = 2 ** 256 - 2 ** 224 - 2 ** 96 + 2 ** 64 - 1
    return x_int % p

if __name__ == '__main__':
##    verify_NTT_Trans()
    h = 0

####   Test Minus Overlap Convolution   ##
##psai = 1321
##wn = 3
##q = 12289
##N = 512
##
##nttp = (wn, N, q)
##
##a = vector_generator(512, 12289)
##b = vector_generator(512, 12289)
##
##v = [(psai ** i) % q for i in range(512)]
##
##a1 = scalar_mul(v, a, q)
##b1 = scalar_mul(v, b, q)
##
##ak = NTT(a1, nttp)
##bk = NTT(b1, nttp)
##
##dk = scalar_mul(ak, bk, q)
##d1 = INTT(dk, nttp)
##
##u = [moddiv1(psai ** i, 1, q) for i in range(512)]
##
##d = scalar_mul(d1, u, q)
##
##
##pad = padding(512)
##a.extend(pad)
##b.extend(pad)
##
##e = NTT_conv(a, b, (3, 1024, 12289) )
##f = []
##
##for i in range(512):
##    f.append ( (e[i] - e[i + 512]) % q )

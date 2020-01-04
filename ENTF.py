import random
import math

#-------------------------------------------#
#    Elementary Number Theory Functions     #
#-------------------------------------------#

########   Greatest Common Divisior      ########
def gcd(a, b):
    if a == 0:
        return b
    else:
        return gcd(b % a, a)

########   Rank2 Matrix Multiplication   ########		
def matmul(m1, m2):
    mm = [[0, 0], [0, 0]]
    mm[0][0] = m1[0][0] * m2[0][0] + m1[0][1] * m2[1][0]
    mm[0][1] = m1[0][0] * m2[0][1] + m1[0][1] * m2[1][1]
    mm[1][0] = m1[1][0] * m2[0][0] + m1[1][1] * m2[1][0]
    mm[1][1] = m1[1][0] * m2[0][1] + m1[1][1] * m2[1][1]
    return mm

########   Extend Euclidean Algorithm    ########
####  Find solutions of ax + by = gcd(a, b)  ####
####  return x, y, gcd(a, b)                 ####
#################################################
def ExEuclid(a, b):
    m=[[1,0],[0,1]]
    while a != 0:
        r = b // a
        a, b = b - r * a, a
        
        ms=[[-r,1],[1,0]]
        m = matmul(m, ms)
    
    x = m[0][1]
    y = m[1][1]
    return x, y, b

########        Modular Inverse          ########
def modinv(a, p):
    x, y, d = ExEuclid(a, p)
    if x < 0:
        return x + p
    else:
        return x
		
########        Modular Devision         ########
####    d = x / y mod m
def moddiv(x, y, m):
    a = y % m
    b = m
    u = x % m
    v = 0
    
    while a != 0:
        r = b // a
        a, b = b - r * a,a
        u, v = (v - r * u) % m, u % m
    
    return v % m
	
########    Modular Square Root          ########
def modsqroot(a, p):
    if (a == 0):
        return 0
    elif (p & 3 == 3):
        return FME(a, (p + 1) >> 2, p)
    elif (p & 7 == 5):
        u = FME(a, (p - 1) >> 2, p)
        if (u == 1):
            return FME(a, (p + 3) >> 3, p)
        else:
            return (FME(a, (p + 3) >> 3, p) * FME(2, (p - 1) >> 2, p)) % p
    else:
        return p8k1(a, p)

def p8k1(a, p):
    for i in range(1, ((p - 1) >> 1) + 1):
        if (i * i % p == a):
            return i

########    Fast Modular Exponent        ########
def FME(a, b, c):
    ans = 1
    base = a
    while (b):
        if (b & 1): ans = ans * base % c
        base = base * base % c
        b = b >> 1
    return ans

########    Probability Prime Test1      ########
def Fermat_test(p, T = 50):
    if ( (p == 1) or (p == 2) or (p == 3)):
        return 'prime'
    else:
        for i in range(2,T+1):
            a = random.randint(2, p - 1)
            s = FME(a, p - 1, p)
            if (s != 1):
                return 'composite'

        return 'prime'
		
########    Probability Prime Test2      ########
def Miller_Labin(p, T = 50):
    if (p == 1) or (p == 2) or (p == 3):
        return 'prime'
    else:
        (s, d) = fac2(p - 1)        
        for i in range(T):
            a = random.randint(2, p - 1)
            base = FME(a, d, p)
            if (base == 1 or base == p - 1):
                continue
            else:
                for i in range(1, s):
                    base = base * base % p
                    if (base == p - 1):
                        break
                if (base == p - 1):
                    continue
                else:
                    return 'composite'
        return 'prime'
		
def fac2(n):
    s = 0
    d = n
    while(d & 1 == 0):
        s += 1
        d >>= 1
    return (s, d)
		
########        Prime Generator          ########
def geneprime(a = 1, b = 2**256):
    while (1):
        p = random.randint(a, b)
        s = Miller_Labin(p, 50)
        if(s == 'prime'):
            return p

########        Prime Table          	 ########
def primetable():
    r = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, \
         67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, \
         139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,    \
         211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277,    \
         281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359,    \
         367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439,    \
         443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521,    \
         523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607,    \
         613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683,    \
         691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773,    \
         787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863,    \
         877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967,    \
         971, 977, 983, 991, 997]
    return r
	
########     Prime Table Generator       ########
######   Generate the first nth prime      ######
def ptn(n):
    if (n <= 168):
        r = primetable()[:n]
    else:
        r = primetable()
        s = 168
        q = 997
        while(s < n):
            q = q + 2
            if (Fermat_test(q) == 'prime'):
                r.append(q)
                s = s + 1
    return r
	
######   Generate the primes less than m   ######
def ptm(m):
    pt = primetable()
    i = 0
    if ( m <= 997 ):
        for per in pt:
            if (per <= m):
                i = i + 1
        r = pt[:i]
    else:
        q = 997
        while (q < m):
            q = q + 2
            r = pt
            if (Fermat_test(q) == 'prime'):
                r.append(q)
    return r
	
########    Chinese Remainder Theory     ########
def CRT(a, m):
    M = mullist(m)
    Mi = []
    ti = []
    for mi in m:
        Mi.append(M // mi)
        ti.append( modinv(Mi[-1], mi) )

    s = 0
    for i in range(len(a)):
        s = s + a[i] * ti[i] * Mi[i]
        
    return s % M

def mullist(r):
    s = 1
    for per in r:
        s = s * per
    return s

	
######## Initial Prime Divisor Factorize ########
def factorize(n):
    sq = math.floor(math.sqrt(n))
    ## slow
    pt = ptm(sq)

    factor = []
    for per in pt:
        if (n % per == 0):
            factor.append(per)
    return factor
	
########     Prime Divisor Factorize     ########
####  We use prime table to trial devision,  ####
####  So it's quite hard to solve number     ####
####  greater than 2 ** 32                   ####
####  A more efficient Algorithm is          ####
####  Pollard-rho Algorithm					 ####
#################################################
def p_resolve(n):
    r = n
    factors = []
    while (Fermat_test(r) == 'composite'):
        factors.append(factorize(r))
        r = r // mullist(factors[-1])
    factors.append([r])

    diff_factors = set(factors[0])
    for per in factors:
        diff_factors |= set(per)
    diff_factors = list(diff_factors)
    diff_num = len(diff_factors)

    diff_power = []
    for i in range(diff_num):
        diff_power.append(0)
    
    for i in range(len(diff_factors)):
        for item in factors:
            if diff_factors[i] in item:
                diff_power[i] += 1
                
    unique_factorize = diff_factors, diff_power
    return unique_factorize
	
########        Pollard-rho             ########
####  Ultrafast prime factorize Algorithm   ####
####  Pollard-rho can resolve about 128 bit ####
####  number in 1 second		    ####
################################################
def pollard_rho(n):
    if Miller_Labin(n, 50) == 'prime':
        return [n]
    else:
        s = []
        c = random.randint(1, n)
        while (1):
            x1 = random.randint(1, n - 1)
            x2 = (x1 * x1 + c) % n
            factor = gcd( (x2 - x1) % n, n)
            if (factor == 1):
                c += 1
                continue
            else:
                break
        s.extend( pollard_rho(factor) )
        u = n // factor
        if (u != 1):
            s.extend( pollard_rho(n // factor) )
        return s

def get_standard_factorize(s):
    diff_factors = sorted(list(set(s)))
    diff_power = []

    for per in diff_factors:
        c = 0
        for item in s:
            if (per == item):
                c += 1
        diff_power.append(c)
		
    return diff_factors, diff_power
	
########        Euler Totient           ########
def phi(n):
    s = pollard_rho(n)
    diff_factors, diff_power = get_standard_factorize(s)
    
    v = n
    for per in diff_factors:
        v = (v // per) * (per - 1)
    return v

########      All Proper Factors        ########
def full_list(r):
    if len(r) == 1:
        s = []
        for i in range(r[0] + 1):
            s.append([i])
        return s
    else:
        v = []
        for i in range(r[-1] + 1):
            u = full_list(r[:-1])
            for per in u:
                per.append(i)
            v.extend(u)
        return v
	
def get_all_proper_factors(unique_factorize):
    diff_factors, diff_power = unique_factorize[0], unique_factorize[1]
    every_power = full_list(diff_power)
    diff_num = len(diff_factors)
    all_factors = []
    for per in every_power:
        s = 1
        for i in range(diff_num):
            s *= diff_factors[i] ** per[i]
        all_factors.append(s)

    all_factors = sorted(all_factors)
    return all_factors
	
########  Find Minimun Primitive Root   ########
def primitive_root(m):
    if m == 2:
        return 1
    elif m == 4:
        return 3
    elif (m & 1):
        euler = phi(m)
        s = pollard_rho(euler)
        unique_factorize = get_standard_factorize(s)
        proper_factors = get_all_proper_factors(unique_factorize)[1:-1]
        
        for i in range(2, m):
            for j in range(len(proper_factors)):
                u = pow(i, proper_factors[j], m)
                if (u == 1):
                    break
            if (j == len(proper_factors) - 1) and (u != 1):
                return i
    else:
        euler = phi(m)
        unique_factorize = p_resolve(euler)
        proper_factors = get_all_proper_factors(unique_factorize)[1:-1]

        for i in range(3, m, 2):
            for j in range(len(proper_factors)):
                u = pow(i, proper_factors[j], m)
                if (u == 1):
                    break
            if (j == len(proper_factors) - 1) and (u != 1):
                return i
	
########        Fermat Numbers          ########
def Fermat(n):
    return 2 ** (2 ** n) + 1
	
########     Quadratic Remainder        ########
########       Legendre Symbol          ########
########        Jacobi Symbol           ########
def LJ1(m):
    return 1

def LJi1(m):
    if (m & 3 == 1):
        return 1
    elif (m & 3 == 3):
        return -1
    else:
        return 0

def LJ2(m):
    if (m & 7 == 1 or m & 7 == 7):
        return 1
    elif (m & 7 == 3 or m & 7 == 5):
        return -1
    else:
        return 0

def quad_reci(m, n):
    return (-1 if (m & 3 == 3 and n & 3 == 3) else 1)

def Jacobi(a, m):
    if (a == 0):
        return 0
    elif (a == 1):
        return LJ1(m)
    elif (a == m - 1):
        return LJi1(m)
    elif (a == 2):
        return LJ2(m)
    else:
        (s, d) = fac2(a)
        h = (LJ2(m) if (s & 1) else 1)
        v = quad_reci(d, m)

        if (d == 1):
            return h * LJ1(m)
        elif (d == m - 1):
            return h * LJi1(m)
        elif (d == 2):
            return h * LJ2(m)
        else:
            return h * v * Jacobi(m % d, d)

def Legendre(x, p):
    return Jacobi(x, p)
	
######## Polynomial Fields Arithmetics  ########
####  This is only implementation on PF     ####
####  Faster Algorithms are not found yet   ####
def GF28_irre():
    return 0b100011011

###   PF_add = PF_sub = (Z2, +) = XOR   ###
def PF_add(x, y):
    return x ^ y

def PF_sub(x, y):
    return x ^ y

def PF_mod(x, p):
    while (len(bin(x)[2:]) > len(bin(p)[2:]) - 1):
        shift_bits = len(bin(x)[2:]) - len(bin(p)[2:])
        x = x ^ (p << shift_bits)
    return x

def PF_modmul(x, y, p):
    shift_y = y
    acc = 0
    n = 0
    while (shift_y):
        if (shift_y & 1):
            acc = PF_add( (x << n), acc)
        n += 1
        shift_y = shift_y >> 1
    return PF_mod(acc, p)

def PF_mul(x, y):
    shift_y = y
    acc = 0
    n = 0
    while (shift_y):
        if (shift_y & 1):
            acc = PF_add( (x << n), acc)
        n += 1
        shift_y = shift_y >> 1
    return acc

def PF_div(y, x, m):
    a = PF_mod(x, m)
    b = m
    u = PF_mod(y, m)
    v = 0
	
    while (a != 0):
        shifts = len(bin(b)[2:]) - len(bin(a)[2:])
        if (shifts >= 0):
            r = shifts
        else:
            r = 0
        a, b = PF_sub(b, a << r), a
        u, v = PF_mod( PF_sub(v, u << r) , m), PF_mod(u, m)
		
    return PF_mod(v, m)

def PF_gcd(x, y):
    if x == 0:
	    return y
    else:
	    return PF_gcd( PF_mod(y, x), x)

def PF_Euclid_div(y, x):
    q = 0
    r = 0
    while (1):
        ly = len(bin(y)[2:])
        lx = len(bin(x)[2:])
        item = ly - lx
        if (item >= 0):
            q += 1 << item
            y = PF_sub(y, x << item)
        else:
            break
    r = y
    return q, r
	
########             End                ########


def gene_d():
    NTT_n = []
    for i in range(Fermat(4), 2 ** 17):
        r = pollard_rho(i)
        s = []
        for per in r:
            s.append(per - 1)

        d = s[0]
        for i in range(len(s) - 1):
            d = gcd(d, s[i+1])
            if (d < 512):
                break

        if (d >= 512) and (d % 512 == 0):
            NTT_n.append(i)
        
    return NTT_n

def nthp_d(lists):
##    lists = gene()
    pts = []
    for per in lists:
        pts.append( primitive_root(per))
        
    phis = []
    for per in lists:
        phis.append( phi(per) )
        
    nthp = []
    for i in range(len(lists)):
        p = lists[i]
        phia = phis[i]
        g = pts[i]
        v = (g ** (phia // 512) ) % p
        nthp.append(v)

    return nthp

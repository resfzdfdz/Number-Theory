import time
import math
import ENTF

#########################################
####       ECSM Implementations      ####
#########################################

###############################################################
####   Affine Coordinate:       y^2 = x^3 + ax + b         ####
####   Projective Coordinate:   zy^2 = axz^2 + bz^3        ####
####   Jacobian Coordinate:     y^2 = x^3 +axz^4 + bz^6    ####
###############################################################

####   Modular Division implemented by Extend Euclidean Algorithm
####   d = b / a mod p
def moddiv(a, b, p):
    u, v, m, n = a % p, p, b % p, 0

    while (u != 1):
        r = v // u
        u, v = v - r * u, u
        m, n = (n - r * m) % p, m
##        print('q =', hex(r))
##        print('u =', hex(u))
##        print('v =', hex(v))
##        print('m =', hex(m))
##        print('n =', hex(n))

    return m

###############################################################
####   Affine Coordinate Elliptic Curve Arithmetic         ####
###############################################################
def PA(P1, P2, p):
    (x1, y1) = P1
    (x2, y2) = P2

    if (P1 == (0, 0) ):
        x3, y3 = P2
    elif (P2 == (0, 0) ):
        x3, y3 = P1
    else:
        k = moddiv(x2 - x1, y2 - y1, p)
        x3 = (k * k - x1 - x2) % p
        y3 = (k * (x1 - x3) - y1) % p

    return x3, y3

def PD(P, p, a):
    (x1, y1) = P

    if (P == (0, 0) ):
        x3, y3 = 0, 0
    else:
        k = moddiv(2 * y1, 3 * x1 * x1 + a, p)  
        x3 = (k * k - x1 - x1) % p
        y3 = (k * (x1 - x3) - y1) % p

    return x3, y3

def ECSM(P, k, p, a):
    index = 255
    Q = (0, 0)
    while (index != -1):
        ref = 2 ** index
        Q = PD(Q, p, a)
        if (k & ref):
            Q = PA(Q, P, p)
        index = index - 1

    return Q
###############################################################
####                Affine Coordinate End                  ####
###############################################################

###############################################################
####                Coordinate Transform                   ####
###############################################################  
def Affine2Projective(P):
    (x, y) = P
    Q = (x, y, 1)
    return Q

def Affine2Jacobian(P):
    (x, y) = P
    Q = (x, y, 1)
    return Q

def Projective2Affine(P, p):
    (X, Y, Z) = P
    x = moddiv(Z, X, p)
    y = moddiv(Z, Y, p)
    return (x, y)

def Jacobian2Affine(P, p):
    (X, Y, Z) = P
    x = moddiv(Z ** 2, X, p)
    y = moddiv(Z ** 3, Y, p)
    return (x, y)
###############################################################
####            Coordinate Transform  End                  ####
###############################################################  

###############################################################
####   Projective Coordinate Elliptic Curve Arithmetic     ####
###############################################################
def PA_P(P, Q, p):
    (x1, y1, z1) = P
    (x2, y2, z2) = Q

    if (P == (0, 1, 0) ):
        x3, y3, z3 = Q
    elif (Q == (0, 1, 0) ):
        x3, y3, z3 = P
    else:
        k1 = x1 * z2 % p
        k2 = x2 * z1 % p
        k3 = (k1 - k2) % p
        k4 = y1 * z2 % p
        k5 = y2 * z1 % p
        k6 = (k4 - k5) % p
        k7 = (k1 + k2) % p
        k8 = z1 * z2 % p
        k9 = k3 * k3 % p
        k10 = k3 * k9 % p
        k11 = (k8 * k6 * k6 - k7 * k9) % p
        
        x3 = k3 * k11 % p
        y3 = (k6 * (k9 * k1 - k11) - k4 * k10) % p
        z3 = k10 * k8 % p

    return (x3, y3, z3)

def PD_P(P, p, a):
    (x1, y1, z1) = P

    if (P == (0, 1, 0)):
        return 0, 1, 0
    else:
        k1 = (3 * x1 * x1 + a * z1 * z1) % p
        k2 = 2 * y1 * z1 % p
        k3 = y1 * y1 % p
        k4 = k3 * x1 * z1 % p
        k5 = k2 * k2 % p
        k6 = (k1 * k1 - 8 * k4) % p
        
        x3 = k2 * k6 % p
        y3 = (k1 * (4 * k4 - k6) - 2 * k5 * k3) % p
        z3 = k2 * k5 % p
        
    return (x3, y3, z3)

def ECSM_P(P, k, p, a):
    index = 255
    Q = (0, 1, 0)
    while (index != -1):
        ref = 2 ** index
        Q = PD_P(Q, p, a)
        if (k & ref):
            Q = PA_P(Q, P, p)
        index = index - 1

    return Q
###############################################################
####            Projective Coordinate End                  ####
###############################################################

###############################################################
####   Jacobian Coordinate Elliptic Curve Arithmetic       ####
###############################################################
def PA_J(P, Q, p):
    (x1, y1, z1) = P
    (x2, y2, z2) = Q
    
    if (P == (1, 1, 0)):
        x3, y3, z3 = Q
    elif (Q == (1, 1, 0) ):
        x3, y3, z3 = P
    else:
        k1 = x1 * z2 * z2 % p
        k2 = x2 * z1 * z1 % p
        k3 = k1 - k2
        k4 = (y1 * z2 ** 3) % p
        k5 = (y2 * z1 ** 3) % p
        k6 = k4 - k5
        k7 = k1 + k2
        k8 = k4 + k5

        x3 = (k6 * k6 - k7 * k3 * k3) % p
        k9 = (k7 * k3 * k3 - 2 * x3) % p
        k10 = (k9 * k6 - k8 * k3 * k3 * k3) % p
        if (k10 & 1):
            y3 = (k10 + p) >> 1
        else:
            y3 = k10 >> 1
        z3 = z1 * z2 * k3 % p

    return x3, y3, z3

def PD_J(P, p, a):
    (x1, y1, z1) = P
    
    if (P == (1, 1, 0)):
        x3, y3, z3 = (1, 1, 0)
    else:
        k1 = (3 * x1 * x1 + a * z1 ** 4) % p
        k2 = (4 * x1 * y1 * y1) % p
        k3 = (8 * y1 ** 4) % p

        x3 = (k1 * k1 - 2 * k2) % p
        y3 = (k1 * (k2 - x3) - k3) % p
        z3 = 2 * y1 * z1 % p

    return x3, y3, z3

def ECSM_J(P, k, p, a):
    index = 255
    Q = (1, 1, 0)
    while (index != -1):
        ref = 2 ** index
        Q = PD_J(Q, p, a)
        if (k & ref):
            Q = PA_J(Q, P, p)
        index = index - 1

    return Q
###############################################################
####                Jacobian Coordinate End                ####
###############################################################


######################################################
####                Curve 25519                   ####
####   Affine Coordinate: y^2 = x^3 + Ax^2 + x    ####
######################################################

##################################################
####    ECSM Implementation of Curve25519     ####
####    and Fp2 Field Operations              ####
##################################################
def Curve25519_Param():
    p = 2 ** 255 - 19
    A = 486662
    return p, A
    
def Fp2_add(x, y):
    p, A = Curve25519_Param()
    
    (a, b) = x
    (c, d) = y

    e = (a + c) % p
    f = (b + d) % p

    return e, f

def Fp2_sub(x, y):
    p, A = Curve25519_Param()

    (a, b) = x
    (c, d) = y

    e = (a - c) % p
    f = (b - d) % p

    return e, f

def Fp2_mul(x, y):
    p, A = Curve25519_Param()

    (a, b) = x
    (c, d) = y

    e = (a * c + 2 * b * d) % p
    f = (a * d + b * c) % p

    return e, f

##   d = y / x mod p in Zp[sqrt(2)] = Fp2
def Fp2_div(x, y):
    p, A = Curve25519_Param()

    (a, b) = x
    (c, d) = y

    num1 = (a * c - 2 * b * d) % p
    num2 = (a * d - b * c) % p
    den = (a * a - 2 * b * b) % p

    e = moddiv(den, num1, p)
    f = moddiv(den, num2, p)

    return e, f

##   Scalar multiplaction in Fp2
def Fp2_smul(k, x):
    p, A = Curve25519_Param()

    (a, b) = x

    e = k * a % p
    f = k * b % p

    return e, f

###############################################################
####                Fp2 Curve25519 Arithmetic              ####
###############################################################
def Fp2_PA(P, Q):
    p, A = Curve25519_Param()

    (x1, y1) = P
    (x2, y2) = Q
    zero = ( (0, 0), (0, 0) )

    if (P == zero):
        return Q
    elif (Q == zero):
        return P
    else:
        t1 = Fp2_sub(y2, y1)
        t2 = Fp2_sub(x2, x1)
        k = Fp2_div(t2, t1)

        t3 = Fp2_mul(k, k)
        t4 = Fp2_sub(t3, (A, 0))
        t5 = Fp2_sub(t4, x1)
        x3 = Fp2_sub(t5, x2)
        
        t6 = Fp2_sub(x1, x3)
        t7 = Fp2_mul(k, t6)
        y3 = Fp2_sub(t7, y1)

    return x3, y3

def Fp2_PD(P):
    p, A = Curve25519_Param()

    (x1, y1) = P
    zero = ( (0, 0), (0, 0) )

    if (P == zero):
        return zero
    else:
        t1 = Fp2_mul(x1, x1)
        t2 = Fp2_smul(3, t1)
        t3 = Fp2_smul(2 * A, x1)
        t4 = Fp2_add(t2, t3)
        num = Fp2_add(t4, (1, 0))
        den = Fp2_smul(2, y1)
        k = Fp2_div(den, num)
        
        t5 = Fp2_mul(k, k)
        t6 = Fp2_sub(t5, (A, 0))
        t7 = Fp2_sub(t6, x1)
        x3 = Fp2_sub(t7, x1)

        t8 = Fp2_sub(x1, x3)
        t9 = Fp2_mul(k, t8)
        y3 = Fp2_sub(t9, y1)

    return x3, y3

def Fp2_ECSM(P, k):
    index = 255
    Q = ((0, 0), (0, 0))
    while (index != -1):
        ref = 2 ** index
        Q = Fp2_PD(Q)
        if (k & ref):
            Q = Fp2_PA(Q, P)
            if (Q[0][1] != 0) or (Q[1][1] != 0):
                print(Q)
        index = index - 1

    return Q

##  X25519 used for Key Sharing
##  K = X25519(a, X25519(b, 9)) = X25519(b, X25519(a, 9))
def ECSM_X25519(a, q):
    p, A = Curve25519_Param()
    alpha = q ** 3 + A * q ** 2 + q
    
    if (q == 9):
        v = 14781619447589544791020593568409986887264606134616475288964881837755586237401
        Q = ( (q, 0), (v, 0) )
    elif ( ENTF.Legendre(alpha, p)  == 1):
        v = ENTF.modsqroot(alpha, p)
        Q = ( (q, 0), (v, 0) )
    else:
        if (alpha & 1):
            r =(alpha + p) >> 1
        else:
            r = alpha >> 1
        v = ENTF.modsqroot(r, p)
        Q = ( (q, 0), (0, v) )

    P = Fp2_ECSM(Q, a)

    return P[0][0] % p
###############################################################
####            Fp2 Curve25519 Arithmetic End              ####
###############################################################

############################################################################
####    Git from https://github.com/nnathan/eccsnacks.git               ####
####    Montgomery Powering Ladder Algorithm for X25519                 ####
####    See Montgomery Powering Ladder in                               ####
####    https://blog.csdn.net/qq_41763108/article/details/88908818      ####
############################################################################
def cswap(swap, x_2, x_3):
    P = 2 ** 255 - 19
    dummy = swap * ((x_2 - x_3) % P)
    x_2 = x_2 - dummy
    x_2 %= P
    x_3 = x_3 + dummy
    x_3 %= P
    return (x_2, x_3)

##  X25519 used for Key Sharing
##  K = X25519(a, X25519(b, 9)) = X25519(b, X25519(a, 9))
def X25519(k, u):
    P = 2 ** 255 - 19
    A24 = 121665
    x_1 = u
    x_2 = 1
    z_2 = 0
    x_3 = u
    z_3 = 1
    swap = 0

    for t in reversed(range(255)):
        k_t = (k >> t) & 1
        swap ^= k_t
        x_2, x_3 = cswap(swap, x_2, x_3)
        z_2, z_3 = cswap(swap, z_2, z_3)
        swap = k_t

        A = x_2 + z_2
        A %= P

        AA = A * A
        AA %= P

        B = x_2 - z_2
        B %= P

        BB = B * B
        BB %= P

        E = AA - BB
        E %= P

        C = x_3 + z_3
        C %= P

        D = x_3 - z_3
        D %= P

        DA = D * A
        DA %= P

        CB = C * B
        CB %= P

        x_3 = ((DA + CB) % P)**2
        x_3 %= P

        z_3 = x_1 * (((DA - CB) % P)**2) % P
        z_3 %= P

        x_2 = AA * BB
        x_2 %= P

        z_2 = E * ((AA + (A24 * E) % P) % P)
        z_2 %= P

    x_2, x_3 = cswap(swap, x_2, x_3)
    z_2, z_3 = cswap(swap, z_2, z_3)

    return (x_2 * pow(z_2, P - 2, P)) % P

################################################

################################################
####    Use Montgomery Ladder to do ECSM    ####
####    in Affine Coornidate y^2=x^3+ax+b   ####
################################################

####    See Montgomery Powering Ladder in
####    https://blog.csdn.net/qq_41763108/article/details/88908818

####    Only X-coordinate Montgomery Powering Ladder Algorithm
####    Finally get X-coordinate of R = [k]P in SM2 Elliptic Curve Arithmetic
def Montgomery_Ladder(G, k, mode = 0):
    p = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF
    a = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC
    b = 0x28e9fa9e9d9f5e344d5a9e4bcf6509a7f39789f515ab8f92ddbcbd414d940e93
    
    Bx = G[0]
    Bz = 1
    Lx = 1
    Lz = 0
    Hx = G[0]
    Hz = 1

    index = 255
    while (index != -1):
        ref = 2 ** index
        kj = k & ref
        if (kj):
            temp_lx = ( (a * Lz * Hz - Lx * Hx) ** 2 - 4 * b * Lz * Hz * (Lx * Hz + Hx * Lz) ) * Bz % p
            temp_lz = ( Bx * (Hx * Lz - Lx * Hz) ** 2 ) % p
            temp_hx = ( (Hx * Hx - a * Hz * Hz) ** 2 - 8 * b * Hx * Hz ** 3 ) % p
            temp_hz = ( Hx * (Hx * Hx + a * Hz * Hz) + b * Hz ** 3) * (4 * Hz) % p
        else:
            temp_lx = ( (Lx * Lx - a * Lz * Lz) ** 2 - 8 * b * Lx * Lz ** 3) % p
            temp_lz = ( Lx * (Lx * Lx + a * Lz * Lz) + b * Lz ** 3) * (4 * Lz) % p
            temp_hx = ( (a * Lz * Hz - Lx * Hx) ** 2 - 4 * b * Lz * Hz * (Lx * Hz + Hx * Lz) ) * Bz % p
            temp_hz = ( Bx * (Hx * Lz - Lx * Hz) ** 2 ) % p
            
        Lx = temp_lx
        Lz = temp_lz
        Hx = temp_hx
        Hz = temp_hz
        index = index - 1

##        if (index < 4):
##            print('Lx =', hex(Lx))
##            print('Lz =', hex(Lz))
##            print('Hx =', hex(Hx))
##            print('Hz =', hex(Hz))


    M1 = (Lx + G[0] * Lz) % p
    M2 = (Lx - G[0] * Lz) % p
    M3 = (Lx * G[0] + a * Lz) % p
    num = ( Hz * ( M1 * M3 + 2 * b * Lz * Lz ) - Hx * M2 * M2 ) % p
    den = (2 * G[1] * Lz ** 2 * Hz) % p

##    print('x3  =', hex(Lx))
##    print('z31 =', hex(Lz))
##    print('y3  =', hex(num))
##    print('z32 =', hex(den))

    xkp = moddiv(Lz, Lx, p)
    ykp = moddiv(den, num, p)

    if (mode == 0):
        return xkp, ykp
    else:
        return Lx, Lz, num, den

#########################################
####       ECSM Test Functions       ####
#########################################
def verify_point(Q):
    (qx, qy) = Q
    k, p, a, b, px, py, n = get_ref()

    left = qy * qy % p
    right = (qx ** 3 + a * qx + b) % p

    if (left == right):
        print('Point on curve')
    else:
        print('Point not on curve')
    
def get_ref():
    k = 0x775513ed70c95a8967b3b780cf3f962c6baccf92d1b7e059110a8391fbccbe21
    p = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF
    a = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC
    b = 0x28e9fa9e9d9f5e344d5a9e4bcf6509a7f39789f515ab8f92ddbcbd414d940e93
    px = 0x32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7
    py = 0xBC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0
    n = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123

    return k, p, a, b, px, py, n

def Curve25519_ref():
    k = 0x77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a
    p = 2 ** 255 - 19
    A = 486662
    u = 9
    v = 14781619447589544791020593568409986887264606134616475288964881837755586237401
    n = 2 ** 252 + 0x14def9dea2f79cd65812631a5cf5d3ed
    return k, p, A, u, v, n

def test_ECSM_ref():
    k, p, a, b, px, py, n = get_ref()
    G = (px, py)
    P = ECSM(G, k, p, a)

    print(hex(P[0]))
    print(hex(P[1]))

def test_ECSM_P():
    k, p, a, b, px, py, n = get_ref()
    G = (px, py)
    GP = Affine2Projective(G)
    PP = ECSM_P(GP, k, p, a)

    P = Projective2Affine(PP, p)

    print(hex(PP[0]))
    print(hex(PP[1]))
    print(hex(PP[2]))

    print(hex(P[0]))
    print(hex(P[1]))

def test_ECSM_J():
    k, p, a, b, px, py, n = get_ref()
    G = (px, py)
    GJ = Affine2Jacobian(G)
    PJ = ECSM_J(GJ, k, p, a)

    P = Jacobian2Affine(PJ, p)

    print(hex(PJ[0]))
    print(hex(PJ[1]))
    print(hex(PJ[2]))

    print(hex(P[0]))
    print(hex(P[1]))    

def test_ECSM_performance():
    k, p, a, b, px, py, n = get_ref()
    G = (px, py)

    N = 100
    start = time.clock()
    for i in range(N):
        P = ECSM(G, k, p, a)
        G = P
    finish = time.clock()

    print('avg time =', (finish - start) / N)

def test_Fp2_ECSM():
    b = 0x5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb
    k, p, A, u, v, n = Curve25519_ref()
    G = ( (u, 0), (v, 0) )
    PB = Fp2_ECSM(G, b)

def test_ECSM_X25519():
    k, p, A, u, v, n = Curve25519_ref()
    a = ENTF.random.randint(1, n)
    N = 50

    start = time.clock()
    for i in range(N):
        b = X25519(a, u)
        a = b
    finish = time.clock()

    print('avg time =', (finish - start) / N)

def test_Montgomary_Ladder():
    k, p, a, b, px, py, n = get_ref()
    G = (px, py)
    for i in range(1000):
        k = ENTF.random.randint(1, 2 ** 256 - 1)
        R0, R1 = ECSM(G, k, p, a)
        R2, R3 = Montgomery_Ladder(G, k)
##        print('R0 =', hex(R0))
##        print('R1 =', hex(R1))
##        print('R2 =', hex(R2))
##        print('R3 =', hex(R3))
##        print('\n')
        if (R0 != R2) or (R1 != R3):
            print('Montgomary Ladder Error!\n')
            print('G = \n')
            print(hex(px))
            print(hex(py))
    

####    Fast Reduction Resolve
def sm2p_resolve(n):
    r = [ [n, 1] ]

    if n < 256:
        return r
    else:
        s = [ [n - 32, 1], [n - 160, 1], [n - 192, -1], [n - 256, 1] ]

        l = []
        for item in s:
            if item[0] >= 256:
                s1 = sm2p_resolve(item[0])
                l.extend(s1)
            else:
                l.append(item)

        t = []                
        for per in l:
            t0 = [sus[0] for sus in t]
            if (per[0] in t0):
                for i in range(len(t)):
                    if (t[i][0] == per[0]):
                        t[i][1] += per[1]
            else:
                t.append(per)

        w = []
        for per in t:
            if (per[1] != 0):
                w.append(per)
        w = sorted(w)
        
        return w

###################################################
####   Elliptic Curve Parameter Generation     ####
###################################################

####   Generate Hilbert class group and class polynomials   ####
def H(D):
    s = math.floor( math.sqrt(D / 3) )

    Ar = []
    MG = []
    for B in range(s + 1):
        facs = ENTF.pollard_rho(D + B * B)
        uq = ENTF.get_standard_factorize(facs)
        afacs = ENTF.get_all_proper_factors(uq)
        Arr = []
        for per in afacs:
            if (per >= 2 * B) and (per <= math.sqrt(D + B * B) ):
                Arr.append(per)

        M = []
        for Ai in Arr:
            C = (D + B * B) // Ai
            if ENTF.gcd(Ai, ENTF.gcd(2 * B, C) ) == 1:
                M.append( (Ai, B, C) )
                if 0 < 2 * B and 2* B < Ai and Ai < C:
                    M.append( (Ai, -B, C) )
        if (M != []):
            MG.append(M)

    hd = 0
    for line in MG:
        for mats in line:
            if (mats != []):
                hd += 1
            
    return MG, hd



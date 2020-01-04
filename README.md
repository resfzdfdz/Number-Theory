# Number Theory Functions
ENTF.py
    --> gcd(a, b): Greatest Common Divisor
    --> moddiv(x, y, m): modular division x/y mod m
    --> modsqroot(a, p): modular square root of a 
    --> Fermat_test(p, T): Fermat Prime Test
    --> Miller_Labin(p, T): Miller Labin Prime Test
    --> geneprime(a, b): Generate a prime in range(a, b)
    --> primetable(): Get primes less than 1000
    --> ptn(n): Get the first nth primes
    --> ptm(m): Get primes less than m
    --> CRT(a, m): Chinese Remainder Theory
    --> pollard_rho(n): Pollard Rho algorithm, prime divisor factorize
    --> get_standard_factorize(s): Standard Factorize
    --> phi(n): Get Euler Totient
    --> get_all_proper_factors(unique_factorize): Get all factors 
    --> primitive_root(m): Get primitive root of modulo m
    --> Fermat(n):  Fermat numbers (2 ** (2 ** n)) + 1
    --> Jacobi(a, m): Jacobi symbol
    --> Legendre(x, p): Legendre symbol
    --> PF_* : GF(2^n) Field Arithmetic

SM2.py
    --> moddiv(a, b, p): modular division b/a mod p
    --> ECSM/PA/PD:   Affine Coornidate Elliptic Curve Arithmetic
    --> ECSM_P/PA_P/PD_P:  Projective Coordinate Elliptic Curve Arithmetic
    --> ECSM_J/PA_J/PD_J:  Jacobian Coordinate Elliptic Curve Arithmetic
    --> Fp2_* : Field Fp2 = Fp(sqrt(2)) Arithmetic and Elliptic Curve Arithmetic
    --> X25519(swap, x_2, x_3): X25519 git from https://github.com/nnathan/eccsnacks.git
    --> Montgomery_Ladder(G, k, mode): Only X-coordinate Montgomery Ladder.
                                       mode = 0: return P = kG = (px, py)
                                       mode = 1: return P = kG = (Lx/Lz, Hx/Hz)
    --> get_ref():  Get parameters of SM2
    --> test_* : Test Functions

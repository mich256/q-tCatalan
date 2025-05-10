from random import randint
from sage.libs.singular.function import singular_function
minbase = singular_function("minbase")

def bivariate_vandermonde(S):
    n = len(S)
    P = PolynomialRing(QQ, ",".join(f"x{i}" for i in [1 .. n]) + "," + ",".join(f"y{i}" for i in [1 .. n]))
    x = P.gens()[:n]
    y = P.gens()[n:]
    M = matrix( [ [ x[i]^a * y[i]^b for i in range(n) ] for a,b in S ])
    return M.det()

def conjectured_basis(n):
    S = []
    for D in DyckWords(n):
        #D.pp()
        a = D.to_area_sequence()
        d = [len([j for j in range(i+1,n) if a[j] == a[i] or a[j] == a[i]-1]) for i in range(n-1)]+[0]
        S.append(list(zip(a,d)))
    return [bivariate_vandermonde(i) for i in S]

def test_basis(n,m, randomized=False):
    Ss = conjectured_basis(n)
    P = Ss[0].parent()
    if randomized:
        p = 4
        while not is_prime(p):
            p = randint(1000, 10000)
        gens = P.gens()
        Q = PolynomialRing(FiniteField(p), gens)
    else:
        Q = P

    I = Q.ideal(Ss)
    return len(minbase(I^m)) == binomial((m+1)*n,n) / (m*n+1)

def test_new(n):
    load('sign-character.sage')
    S = test(n)
    P = S[0].parent()
    I = ideal(S)
    return minbase(I)

def Ddict(n,m):
    load('qtCatalan.sage')
    d = {}
    for D in Dyck_paths(n*m,n):
        d[(D.area(),D.dinv())] = D
    return d

def Dpprint(n,m, randomized=False):
    Ss = conjectured_basis(n)
    P = Ss[0].parent()
    if randomized:
        p = 4
        while not is_prime(p):
            p = randint(1000, 10000)
        gens = P.gens()
        Q = PolynomialRing(FiniteField(p), gens)
    else:
        Q = P

    I = Q.ideal(Ss)
    M = minbase(I^m)
    qtCat = 0
    qtring.<q,t> = QQ['q,t']
    x = Q.gens()[:n]
    y = Q.gens()[n:]
    d = Ddict(n,m)
    for f in M:
        m = f.monomials()[0]
        xdeg = sum( m.degree(x[i]) for i in range(n))
        ydeg = sum( m.degree(y[i]) for i in range(n))
        print(factor(f))
        D = d[(xdeg,ydeg)]
        print(D.dinv_code())
        D.pp()
        qtCat += q^(xdeg) * t^(ydeg)
    return qtCat
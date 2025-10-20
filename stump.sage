from random import randint
from sage.libs.singular.function import singular_function
minbase = singular_function("minbase")

def bivariate_vandermonde(S):
    n = len(S)
    P = PolynomialRing(QQ, ",".join(f"x{i}" for i in [1 .. n]) + "," + ",".join(f"y{i}" for i in [1 .. n]))
    x = P.gens()[:n]
    y = P.gens()[n:]
    M = matrix([[x[i]^a * y[i]^b / (factorial(a) * factorial(b)) for i in range(n)] for a,b in S])
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

def Ddict(n):
    d = {}
    for D in DyckWords(n):
        d[D.area(),D.bounce()] = D
    return d

def Dad(n,m):
    load('qtCatalan.sage')
    d = {}
    for D in Dyck_paths(n*m,n):
        d[(D.area(),D.bounce())] = D
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
    load('qtCatalan.sage')
    tt = qtCatalan(n*m,n)
    R = tt.parent()
    qtCat = 0
    x = Q.gens()[:n]
    y = Q.gens()[n:]
    d = Dad(n,m)
    dd = Ddict(n)
    for f in M:
        m = f.monomials()[0]
        xdeg = [m.degree(x[i]) for i in range(n)]
        ydeg = [m.degree(y[i]) for i in range(n)]
        D = d[(sum(xdeg),sum(ydeg))]
        D.pp()
        print(D.bounce_sequence())
        L = factor(f)
        print(L)
        temp = L[0][0].monomials()
        if len(temp) != 2:
            m1 = temp[0]
            x1 = [m1.degree(x[i]) for i in range(n)]
            y1 = [m1.degree(y[i]) for i in range(n)]
            D = dd[sum(x1),sum(y1)]
            D.pp()
            print(D.bounce_path().touch_points())
        elif temp[0].degree(x[0]) == 1:
            continue
        else:
            continue
        temp = L[-1][0].monomials()
        if len(temp) != 2:
            m2 = temp[0]
            x2 = [m2.degree(x[i]) for i in range(n)]
            y2 = [m2.degree(y[i]) for i in range(n)]
            D = dd[sum(x2),sum(y2)]
            D.pp()
            print(D.bounce_path().touch_points())
        elif temp[0].degree(x[0]) == 1:
            continue
        else:
            continue
    return
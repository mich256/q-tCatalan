def gen_ring(n):
    # Use strings, not SR vars; use a graded order for nicer behavior in Singular
    names = [f"x{i}" for i in range(1, n+1)] + [f"y{i}" for i in range(1, n+1)]
    R = PolynomialRing(QQ, names, order='degrevlex')
    x = R.gens()[:n]
    y = R.gens()[n:]
    return R, x, y

def powersum(x, y, i, j):
    n = len(x)
    return sum(x[k]**i * y[k]**j for k in range(n))

def gen_ideal(R, x, y):
    n = len(x)
    gens = [powersum(x, y, i, j) for i in range(2*n) for j in range(2*n)]
    return R.ideal(gens[1:])

def Poly2Vec(MonList,pol):
    """This function converts polymonial pol into the vector by extracting coeffients in front of monomials from the list MonList."""
    return([pol.monomial_coefficient(z) for z in MonList])

def gb_init(n):
    return gen_ideal(n).normal_basis()

def test_rank(n):
    ListV = [Poly2Vec(gb_init(n), x.reduce(gen_ideal(n))) for x in test(n)]
    return matrix(ListV)

def gen_m(R, x, y, l):
    n = len(l)
    return matrix([[x[i]^(l[j][0])*y[i]^(l[j][1]) for j in range(n)] for i in range(n)])

def gen_det(R, x, y, l):
    return gen_m(R, x, y, l).determinant()

def dinv_code(D):
    n = D.semilength()
    a = D.to_area_sequence()
    return [len([j for j in range(i+1,n) if a[j] == a[i] or a[j] == a[i]-1]) for i in range(n-1)]+[0]

def catdet(R, x, y):
    n = len(x)
    for D in DyckWords(n):
        a = D.to_area_sequence()
        aa = dinv_code(D)
        yield gen_det(R, x, y, list(zip(aa,a)))

def alternate_ideal(R, x, y):
    n = len(x)
    return R.ideal(list(catdet(R, x, y)))

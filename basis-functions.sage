load('sign-character.sage')

def powersum(l1,l2,i,j):
    n = len(l1)
    return sum([l1[k]^i*l2[k]^j for k in range(n)])

def gen_ideal(n):
    R = gen_ring(n)
    xs = R.gens()[:n]
    ys = R.gens()[n:]
    #R.inject_variables(verbose = False)
    return R.ideal([powersum(xs,ys,i,j) for i in range(2*n) for j in range(2*n)][1:])

def Poly2Vec(MonList,pol):
    """This function converts polymonial pol into the vector by extracting coeffients in front of monomials from the list MonList."""
    return([pol.monomial_coefficient(z) for z in MonList])

def test_rank(n):
    R = gen_ring(n)
    JGB = R.ideal(gen_ideal(n).groebner_basis())
    l_quot = JGB.normal_basis()
    ListV = [Poly2Vec(l_quot, x.reduce(JGB)) for x in test(n)]
    return rank(matrix(ListV))
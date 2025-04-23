n = 3

# pp is random prime for defining base field
pp = random_prime(1000)

variables = ['x'+str(i) for i in range(1,n+1)] + ['y'+str(i) for i in range(1,n+1)]

load('sign-character.sage')

#R = gen_ring(n)
R = gen_ring(n,GF(pp))
#R.inject_variables(verbose = False)
 
def Pij(i,j):
    """Here we define the power sum of x_k^iy_k^j"""
    return(sum([R.gens()[k]^i*R.gens()[k+n]^j for k in range(n)]))

# J is the defining ideal for the double coinvariants
J = ideal([Pij(i,j) for i in range(2*n) for j in range(2*n)][1:])
JGB = ideal(J.groebner_basis())

# Groebner basis for DH
l_quot = JGB.normal_basis()

def Poly2Vec(MonList,pol):
    """This function converts polymonial pol into the vector by extracting coeffients in front of monomials from the list MonList."""
    return([pol.monomial_coefficient(z) for z in MonList])


S = SymmetricGroup(n)

def DescMon(g):
    """This function computes the descent monomials in x variables"""
    z = R.gens()
    p = Permutation(g)
    return R.prod(prod(z[p[j]-1] for j in range(i+1)) for i in range(len(p) - 1) if p[i] > p[i + 1])


def Desc(g):
  """ The list of descents of the g from symmetric group S"""  
  ll = []
  for i in range(1,n):
      if g(i) > g(i+1):
          ll = ll+[i]
  return(ll)

def w_tau(g):
    """The list of w_i(g) from the first version of our paper"""
    ll = []
    DesPrime =[0]+ Desc(g)+[n]
    for m in range(0,len(DesPrime)-2):
        for i in range(DesPrime[m]+1,DesPrime[m+1]+1):
            cc = 0
            for jj in range(DesPrime[m+1]+1,DesPrime[m+2]+1):
                if g(i) > g(jj):
                    cc = cc + 1
            ll = ll + [DesPrime[m+1]-i+cc]
    for i in range(DesPrime[-2]+1,n+1):
       ll = ll + [n-i+1]
    return(ll)


# Generating function for our monomial basis
MonBasisGenF = sum(DescMon(Permutation(g))* prod(sum(R.gen(g(i+1)+n-1)^k for k in range(w_tau(g)[i])) for i in range(n)) for g in S)

# The list mononial in our basis
ListMonBasis = MonBasisGenF.monomials()


# The matrix of coordinates of our monomial basis in Groebber basis
ListV = [Poly2Vec(l_quot, mm.reduce(JGB)) for mm in ListMonBasis]
MM = matrix(ListV)


# Here we check that our theorem is valid. That is we have a basis of DH
print('The rank of subspace of spanned by our monomial basis is '+str(rank(MM)))
print('The dimension of the space of the double coinvariants is '+ str(JGB.vector_space_dimension()))



def Act(g,mm):
    """Action of g from the symmetric group S on the monomial mm"""
    degs = list(mm.exponents()[0])
    return(prod(R.gen(g(i+1)-1)^(degs[i]) for i in range(n))*prod(R.gen(g(i+1)+n-1)^(degs[i+n]) for i in range(n)))

def AntiSym(mm):
    """Antisymmetrization of the monomial mm"""
    return(sum(Act(g,mm)*sign(g) for g in S))
    
CatGens = [ mm for mm in ListMonBasis if AntiSym(mm) != 0]

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3



def ClearRepeats(ll):
    """It keeps one monomial per S orbit"""
    aux = ll
    aux1 = []
    flag = True
    while flag:
        flag = False
        for x in aux:
            if len(intersection([Act(g,x) for g in S],aux))>1:
                flag = True
                aux.remove(x)
    return(aux)

    
CatGens1 = ClearRepeats(CatGens)


# This a check that we get a correct number of elements in the Catalan basis
print('Number of elements in Catalan subset of monomials is '+ str(len(CatGens1)))

print('The Catalan number is '+str(catalan_number(n)))

# Check of linear independence of the AntiSyms of our monomials
ListCatalan = [Poly2Vec(l_quot, AntiSym(mm).reduce(JGB)) for mm in CatGens1]
MM = matrix(ListCatalan)

print('The rank of the subspace spanned by the antisyms of our monomials is '+str(rank(MM)))

M = ideal([AntiSym(mm) for mm in CatGens])
for gb in ideal(ideal(M^2).groebner_basis()).gens():
    temp = gb.lm().exponents()[0]
    list(zip(temp[:3],temp[3:]))



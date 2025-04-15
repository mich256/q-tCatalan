load('sign-character.sage')

def cox_gen_ring(W):
    n = W.rank()
    return gen_ring(n)

def act(g,mm):
    """Action of g from the symmetric group S on the monomial mm"""
    degs = list(mm.exponents()[0])
    R = mm.parent()
    n = R.ngens()//2
    return(prod(R.gen(g(i+1)-1)^(degs[i]) for i in range(n))*prod(R.gen(g(i+1)+n-1)^(degs[i+n]) for i in range(n)))

def symmetrizer(W,mm):
    return(sum(act(g,mm) for g in W))

def antisymmetrizer(W,mm):
    s = W.sign_representation()
    """Antisymmetrization of the monomial mm"""
    return(sum(act(g,mm)*s.sign_function(g) for g in W))
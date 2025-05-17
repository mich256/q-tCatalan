n = 3

x = [var('x_%d'%i) for i in range(1,n+1)]
y = [var('y_%d'%i) for i in range(1,n+1)]
R = PolynomialRing(QQ, x+y, order = 'lex')

vn = det(matrix.vandermonde(x))

def comp(a,b):
	return lambda x: a(b(x))

def E(p):
	return lambda f: sum(y[j] * f.derivative(x[j], p) for j in range(n))

def psm(mu):
	temp = vn
	for i in mu:
		temp = E(i)(temp)
	return temp

Sym = SymmetricFunctions(QQ)
s = Sym.schur()
p = Sym.powersum()

def sch(l):
	expanded = p(s(l))
	parts = expanded.support()
	return sum(expanded.coefficient(par) * psm(par) for par in parts)

def diagram(D):
	P = [i - D.to_area_sequence()[i] for i in range(D.semilength())]
	return Partition(reversed(P))

def allen(n):
	return [sch(diagram(D)) for D in DyckWords(n)]
n = 4

x = [var('x_%d'%i) for i in range(1,n+1)]
y = [var('y_%d'%i) for i in range(1,n+1)]
R = PolynomialRing(QQ, x+y, order = 'lex')

vn = prod(x[i]-x[j] for i in range(n) for j in range(i))

def gen_m(l):
	return matrix([[x[i]^(l[j][0])*y[i]^(l[j][1]) for j in range(n)] for i in range(n)])

def gen_det(l):
	return gen_m(l).determinant()

def gen_ideal(n):
	return R.ideal([(sum([x[k]^i*y[k]^j for k in range(n)])) for i in range(2*n) for j in range(2*n)][1:])

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
e = Sym.elementary()
h = Sym.homogeneous()

def sch(l):
	expanded = p(s(l))
	print(l,expanded)
	parts = expanded.support()
	return sum(expanded.coefficient(par) * psm(par) for par in parts)

def diagram(D):
	P = [i - D.to_area_sequence()[i] for i in range(D.semilength())]
	return Partition(reversed(P))

def dinv_code(D):
	n = D.semilength()
	a = D.to_area_sequence()
	return [len([j for j in range(i+1,n) if a[j] == a[i] or a[j] == a[i]-1]) for i in range(n-1)]+[0]

def allen(n):
	return [sch(diagram(D)) for D in DyckWords(n)]

def prodfact(tu):
	return prod(factorial(i) for i in tu)

def mon_to_diff(m, f):
	temp = f
	for v in m.variables():
		temp = f.derivative(v, m.degree(v))
	return temp

def poly_to_diff(p, f):
	return sum(coeff * mon_to_diff(mon, f) for coeff, mon in p)
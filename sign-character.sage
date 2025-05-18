#load('Basis.sage')
load('qtCatalan.sage')
import itertools

n = 3

x = [var('x_%d'%i) for i in range(1,n+1)]
y = [var('y_%d'%i) for i in range(1,n+1)]
PolynomialRing(QQ, x+y, order = 'lex')

def gen_m(l):
	return matrix([[x[i]^(l[j][0])*y[i]^(l[j][1]) for j in range(n)] for i in range(n)])

def gen_det(l):
	return gen_m(l).determinant()

def gen_ideal(n):
	R = gen_ring(n)
	x = R.gens()[:n]
	y = R.gens()[n:]
	return ideal([(sum([x[k]^i*y[k]^j for k in range(n)])) for i in range(2*n) for j in range(2*n)][1:])

def dinv_code(D):
	n = D.semilength()
	a = D.to_area_sequence()
	return [len([j for j in range(i+1,n) if a[j] == a[i] or a[j] == a[i]-1]) for i in range(n-1)]+[0]

def normalcat(n):
	t = []
	for D in DyckWords(n):
		a = D.to_area_sequence()
		aa = dinv_code(D)
		# t.append(list(zip(a,aa)))
		temp = gen_det(list(zip(a,aa)))
		if temp == 0:
			raise Exception('zero determinant')
		else:
			t.append(temp)
	return t

def pf(D):
	s = D.to_area_sequence()
	n = len(s)
	# Create a list of (index, label) pairs
	indexed = list(enumerate(s))
	# Sort by label, maintaining original order (stable sort)
	indexed.sort(key=lambda x: x[1])
	# Assign numbers 1 to n in this order
	perm = [0] * n
	for i, (idx, _) in enumerate(indexed):
		perm[idx] = i + 1
	return ParkingFunction(labelling=perm,area_sequence=s)

def pfmaj(pf):
	w = pf.to_labelling_permutation()
	a = pf.to_area_sequence()
	return [a[w(i)-1] for i in range(1,len(a)+1)]

def pfdinv(pf):
	w = pf.to_labelling_permutation()
	a = pf.to_area_sequence()
	n = len(a)
	return [len([j for j in range(i+1,n) if (a[j] == a[i] and w(j+1) > w(i+1)) or (a[j] == a[i] - 1 and w(j+1) < w(i+1))]) for i in range(n-1)]+[0]

def test(n):
	t = set()
	for D in DyckWords(n):
		pfd = pf(D)
		t2 = list(zip(pfdinv(pfd),D.to_area_sequence()))
		t.add(gen_det(t2))
	return t
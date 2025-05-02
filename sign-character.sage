#load('Basis.sage')
load('qtCatalan.sage')
import itertools

def gen_ring(n, F = QQ):
	xs = [var('x_%d'%i) for i in range(1,n+1)]
	ys = [var('y_%d'%i) for i in range(1,n+1)]
	return PolynomialRing(F, xs+ys, order = 'lex')

def gen_m(l):
	n = len(l)
	R = gen_ring(n)
	xy = R.gens()
	return matrix([[xy[i]^(l[j][0])*xy[n+i]^(l[j][1]) for j in range(n)] for i in range(n)])

def gen_det(l):
	return gen_m(l).determinant()

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

def pfmaj(pf):
	w = pf.to_labelling_permutation()
	a = pf.to_area_sequence()
	return [a[w(i)-1] for i in range(1,len(a)+1)]

def pfdinv(pf):
	w = pf.to_labelling_permutation()
	a = pf.to_area_sequence()
	n = len(a)
	return [len([j for j in range(i+1,n) if (a[j] == a[i] and w(j+1) > w(i+1)) or (a[j] == a[i] - 1 and w(j+1) < w(i))]) for i in range(n-1)]+[0]

def test(n):
	for pf in ParkingFunctions(n):
		temp = gen_det(list(zip(pfmaj(pf),pfdinv(pf))))
		if temp != 0:
			print(pfmaj(pf),pfdinv(pf))
			pf.pretty_print()
	return

def test1(n):
	for D in DyckWords(n):
		pf = ParkingFunction(labelling = list(range(1,n+1)), area_sequence=D.to_area_sequence())
		t1 = list(zip(D.to_area_sequence(),dinv_code(D)))
		t2 = list(zip(pfmaj(pf),pfdinv(pf)))
		temp = gen_det(t2)
		print(t1,t2)
		if temp == 0:
			raise Exception('zero')
	return
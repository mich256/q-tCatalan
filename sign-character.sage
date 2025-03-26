def gen_ring(n):
	xs = [var('x_%d'%i) for i in range(1,n+1)]
	ys = [var('y_%d'%i) for i in range(1,n+1)]
	return PolynomialRing(QQ, xs+ys)

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

def test(n):
	t = []
	for D in DyckWords(n):
		a = D.to_area_sequence()
		aa = dinv_code(D)
		t.append(gen_det(list(zip(a,aa))))
	return t

def gen_m_from_parts(p):
	return gen_m(Partition(p).cells())

def gen_det_from_parts(p):
	return gen_det(Partition(p).cells())
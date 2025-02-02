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

def gen_m_from_parts(p):
	return gen_m(Partition(p).cells())

def gen_det_from_parts(p):
	return gen_det(Partition(p).cells())


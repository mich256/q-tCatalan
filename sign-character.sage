#load('Basis.sage')
load('qtCatalan.sage')

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
		#t.append(list(zip(a,aa)))
		temp = gen_det(list(zip(a,aa)))
		if temp == 0:
			raise Exception('zero determinant')
		else:
			t.append(temp)
	return t

def test2(n):
	t = []
	for D in DyckWords(n):
		a = D.to_area_sequence()
		tp = [0]*n
		for i in range(n):
			if i in D.bounce_path().touch_points():
				tp[n-i-1] = i
		temp = gen_det(list(zip(a,tp)))
		if temp == 0:
			#raise Exception('zero')
			print(a,tp)
			D.bounce_path().pp()
		else:
			print(a,tp)
			t.append(temp)
	#return t

def rationalcat(n,m):
	t = []
	for D in Dyck_paths(m,n):
		a = D.area_sequence()
		d = D.dinv_code()
		temp = gen_det(list(zip(a,d)))
		if temp == 0:
			#raise Exception('zero determinant')
			print('*')
			continue
		else:
			t.append(temp)
	return t
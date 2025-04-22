load('sign-character.sage')

def schedule(w):
	n = w.size()
	r = w.runs()
	v = w.inverse()
	r.append([0])
	s = [0]*n
	for j in range(len(r)-1):
		for i in r[j]:
			s[v(i)-1] = len([k for k in r[j] if k > i]) + len([k for k in r[j+1] if k < i])-1
	return s

def dinv_code(D):
	n = D.semilength()
	a = D.to_area_sequence()
	return [len([j for j in range(i+1,n) if a[j] == a[i] or a[j] == a[i]-1]) for i in range(n-1)]+[0]

def maj(w):
	des = w.descents()
	v = w.inverse()
	n = len(w)
	return [len([j for j in des if j >= v(i)]) for i in range(1,n+1)]

def CO_basis(n):
	R = gen_ring(n)
	v = R.gens()
	xs = v[:n]
	ys = v[n:]
	return sum(prod(ys[i]^(maj(w)[i]) * sum(xs[w[i]-1]^k for k in range(schedule(w)[i]+1)) for i in range(n)) for w in Permutations(n)).monomials()
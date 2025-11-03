load('sign-character.sage')
load('bounce.sage')

def schedule(w):
	n = w.size()
	r = w.runs()
	r.append([0])
	v = w.inverse()
	s = [0]*n
	for j in range(len(r)-1):
		for i in r[j]:
			s[v(i)-1] = len([k for k in r[j] if k > i]) + len([k for k in r[j+1] if k < i])
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
	return sum(prod(ys[i]^(maj(w)[i]) * sum(xs[w[i]-1]^k for k in range(schedule(w)[i])) for i in range(n)) for w in Permutations(n)).monomials()

def dw_latex(dw, aa = False, dd = False, bb = False):
	n = dw.semilength()
	res = '\\begin{tikzpicture}[scale=0.5]\n'
	res += '\\draw[dotted] (0,0) grid (%d,%d);\n' % (n,n)
	res += '\\draw[thick] (0,0)'
	stats = ''
	coord = [0,0]
	for i in range(len(dw)):
		if dw[i] == 1:
			coord[1] += 1
		if dw[i] == 0:
			coord[0] += 1
		res += '--(%d,%d)'% tuple(coord)
	if aa:
		stats += '\\draw node at (1,-.5) {area: %s};\n' % str(dw.to_area_sequence())
	if dd:
		stats += '\\draw node at (1,-1.5) {dinv: %s};\n' % str(dinv_code(dw))
	if bb:
		stats += '\\draw node at (1,-2) {bounce: %s};\n' % str(bounce_seq(dw))
	print(res + ';\n' + stats + '\\end{tikzpicture}')

def pf_dinv_code(pf):
	n = len(pf)
	d = [0]*n
	for i,j in pf.dinversion_pairs():
		d[i] += 1
	return d

def pf_latex(pf, aa = False, dd = False):
	w = pf.to_labelling_permutation()
	dw = pf.to_dyck_word()
	n = len(pf)
	res = '\\begin{tikzpicture}[scale=0.5]\n'
	res += '\\draw[dotted] (0,0) grid (%d,%d);\n' % (n,n)
	res += '\\draw[thick] (0,0)'
	label = ''
	stats = ''
	coord = [0,0]
	for i in range(2*n):
		if dw[i] == 1:
			coord[1] += 1
			label += '\\draw node at (%f,%f) {%d};\n' % (coord[0]+0.5,coord[1]-0.5, w(coord[1]))
		if dw[i] == 0:
			coord[0] += 1
		res += '--(%d,%d)'% tuple(coord)
	if aa:
		stats += '\\draw node at (1,-.5) {%s};\n' % str(pf.to_area_sequence())
	if dd:
		stats += '\\draw node at (1,-1.5) {%s};\n' % str(pf_dinv_code(pf))
	print(res + ';\n' + label + stats + '\\end{tikzpicture}')

def bounce_seq(d):
	temp = []
	l = d.bounce_path().touch_points()
	for j in range(len(l)):
		if j == 0:
			temp += l[j]*[0]
		else:
			temp += (l[j]-l[j-1])*[j]
	return temp

def bounce_code(d):
	l = d.bounce_path().touch_points()
	n = l[-1]
	return [n-l[j] for j in range(len(l)-1)[::1]]

def pf_latex(pf, aa = False, dd = False):
	w = pf.to_labelling_permutation()
	dw = pf.to_dyck_word()
	n = len(pf)
	res = '\\begin{tikzpicture}[scale=0.5]\n'
	res += '\\draw[dotted] (0,0) grid (%d,%d);\n' % (n,n)
	res += '\\draw[thick] (0,0)'
	label = '\\draw node at (0.5,0.5) {%d};\n' % w[0]
	stats = ''
	coord = [0,0]
	for i in range(len(dw)):
		if dw[i] == 1:
			coord[1] += 1
			label += '\\draw node at (%f,%f) {%d};\n' % (coord[0]+0.5,coord[1]-0.5, w(coord[1]))
		if dw[i] == 0:
			coord[0] += 1
		res += '--(%d,%d)'% tuple(coord)
	if aa:
		stats += '\\draw node at (1,-.5) {area: %d};\n' % pf.area()
	if dd:
		stats += '\\draw node at (1,-1.5) {dinv: %d};\n' % pf.dinv()
	return res + ';\n' + label + stats + '\\end{tikzpicture}\n'
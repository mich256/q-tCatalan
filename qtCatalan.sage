from itertools import chain

def get_intervals(binary_list):
	if not binary_list:
		return []

	intervals = []     # This will store the run lengths.
	current_value = binary_list[0]
	if current_value != 0 and current_value != 1:
		raise Exception('input a binary list')
	count = 1
	# Start from the second element of the list.
	for item in binary_list[1:]:
		if item != 0 and item != 1:
			raise Exception('input a binary list')
		if item == current_value:
			count += 1
		else:
			intervals.append(count)
			current_value = item
			count = 1

	# Append the final run length.
	intervals.append(count)
	return intervals

def sc(m,p):
	if 1 <= p <= m:
		return m+1-p
	if -m <= p <= 0:
		return m+p
	else:
		return 0

class RationalDyckPath:
	def __init__(self, l: list):
		self.DyckWord = l
		self.vertical = sum(l)
		self.horizontal = len(l) - sum(l)
		self.slope = self.vertical/self.horizontal
		temp = get_intervals(self.DyckWord)
		v = temp[::2]
		self.cumulative_v = [sum(v[:i+1]) for i in range(len(v))]
		h = temp[1:][::2]
		self.cumulative_h = [sum(h[:i+1]) for i in range(len(h))]
		
	def area_sequence(self):
		v = self.cumulative_v
		h = self.cumulative_h
		t = [max(0,floor(y/self.slope)) for y in range(v[0])]
		for i in range(1,len(v)):
			t += [max(0,floor(y/self.slope) - h[i-1]) for y in range(v[i-1],v[i])]
		return t

	def diagram(self):
		P = [floor(i / self.slope) - self.area_sequence()[i] for i in range(self.vertical)]
		P = Partition(reversed(P))
		#P.pp()
		return P

	def dinv_boxes(self):
		La = self.diagram()
		t = []
		tt = [[' ']*2*i for i in La] + ['']
		for (i,j) in La.cells():
			a = La.arm_length(i,j)
			l = La.leg_length(i,j)
			if l/(a+1) <= self.slope and a/(l+1) < 1/self.slope:
				t.append((i,j))
				tt[i][2*j] = '.'
		return t, tt

	def dinv_code(self):
		a = self.area_sequence()
		m = 1/self.slope
		n = self.vertical
		return [sum(sc(m,a[i] - a[j]) for j in range(i+1,n)) for i in range(n-1)] + [0]

	def dinv(self):
		return len(self.dinv_boxes()[0])

	def reverse(self):
		return RationalDyckPath([1-i for i in self.DyckWord[::-1]])

	def pp(self) -> None:
		n = len(self.DyckWord)
		if n == 0:
			return
		temp = get_intervals(self.DyckWord)
		tt = len(temp)
		a = self.area_sequence()
		k = self.vertical-1
		st = self.dinv_boxes()[1]
		#x = sum([2*temp[i]-1 for i in range(1,tt-2)[::2]]) + tt//2
		if temp[0] == self.vertical:
			s = ' ' + (2*self.horizontal-1)* '_' + '\n'
		else:
			s = ''.join(st[0]) + ' ' + (2*temp[-1]-1)* '_' + '\n'
		for i in list(reversed(range(tt-1)))[::2]:
			if i == 0:
				for j in range(temp[0]):
					s += '|' + a[k]* ' x' + '\n'
					k -= 1
				print(s)
				return
			for j in range(temp[i]-1):
				s += ''.join(st[self.vertical-k])+ '|' + a[k]* ' x' + '\n'
				k -= 1
			s += ''.join(st[self.vertical-k]) + ' ' + (2*temp[i-1]-1)* '_' + '|' + a[k]*' x' + '\n'
			k -= 1

	def area(self):
		return sum(self.area_sequence())

	def rank(self):
		r = [0]
		counter = 0
		h = self.horizontal
		v = self.vertical
		g = gcd(h,v)
		for i in self.DyckWord:
			if i == 1:
				counter += h
			else:
				counter -= v
			r.append(counter)
		r.pop()
		return r

	def zeta(self):
		B = self.rank()
		sorted_indices = sorted(range(len(B)), key=lambda i: (B[i], -i))
		rearranged_A = [self.DyckWord[i] for i in sorted_indices]
		return RationalDyckPath(rearranged_A)

	def bounce_v(self):
		v = []
		h = [0]
		horizontal_sum = 0
		vertical_sum = 0
		index = 0
		while horizontal_sum < self.horizontal and vertical_sum < self.vertical:
			v.append(self.cumulative_v[index] - vertical_sum)
			vertical_sum = self.cumulative_v[index]
			if len(v) <= 1/self.slope:
				h.append(h[-1]+v[-1])
			else:
				h.append(h[-1]-v[-1/self.slope-1]+v[-1])
			horizontal_sum += h[-1]
			while horizontal_sum >= self.cumulative_h[index]:
				index += 1
		return v

	def bounce_sequence(self):
		v = self.bounce_v()
		return [j for j in range(len(v)) for i in range(v[j])]

	def bounce_code(self):
		v = self.bounce_v()
		return Composition(v[1:][::-1]).partial_sums(final = True)

	def bounce(self):
		return sum(self.bounce_sequence())

	def latex(self):
		res = '\\begin{tikzpicture}[scale=0.5]\n'
		res += '\\draw[dotted] (0,0) grid (%d,%d);\n' % (self.horizontal,self.vertical)
		res += '\\draw[thick] (0,0)'
		i,j = 0,0
		for k in self.DyckWord:
			if k == 1:
				j += 1
			else:
				i += 1
			res += '--(%d,%d)'%(i,j)
		print(res + ';\n' + '\\end{tikzpicture}')

def Dyck_paths(h: int, v: int):
	t = []
	slope = v/h
	for x in Subsets(h+v,h):
		tt = [0 if y in x else 1 for y in range(1,v+h+1)]
		switch = True
		horizontal_sum = 0
		vertical_sum = 0
		for j in range(v+h):
			horizontal_sum += 1-tt[j]
			vertical_sum += tt[j]
			if horizontal_sum == 0:
				continue
			elif vertical_sum/horizontal_sum < v/h:
				switch = False
				break
		if switch:
			t.append(RationalDyckPath(tt))
	return t

def qtCatalan(a: int, b: int):
	R.<q,t> = QQ['q,t']
	Cqt = 0
	for D in Dyck_paths(a, b):
		Cqt += q^(D.area()) * t^(D.dinv())
	return Cqt

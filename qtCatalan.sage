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

class RationalDyckPath:
	def __init__(self, l: list):
		self.DyckWord = l
		self.vertical = sum(l)
		self.horizontal = len(l) - sum(l)
		# if gcd(self.vertical,self.horizontal) != 1:
		# 	raise Exception('coprime')
		self.slope = self.vertical/self.horizontal
		# self.vertices = [(0,0)]
		self.dinv_boxes = []
		temp = get_intervals(self.DyckWord)
		v = temp[::2]
		self.cumulative_v = [sum(v[:i+1]) for i in range(len(v))]
		h = temp[1:][::2]
		self.cumulative_h = [sum(h[:i+1]) for i in range(len(h))]
		self.dinv_box_strs = list(chain.from_iterable([[' ']*(1+sum([2*h[j] for j in range(len(h)-i-1)]))] + [[' ']*sum([2*h[j] for j in range(len(h)-i-1)])]*(max(v[len(h)-i-1]-1,0)) for i in range(len(h)-1))) + [' ']

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

	def dinv_code(self):
		La = self.diagram()
		t = [0]*self.vertical
		h = self.cumulative_h
		for (i,j) in La.cells():
			a = La.arm_length(i,j)
			l = La.leg_length(i,j)
			if l/(a+1) < self.slope and ((a == 0 and l+1 >= self.slope) or (a != 0 and (l+1)/a >= self.slope)):
				self.dinv_boxes.append((i,j))
				self.dinv_box_strs[i][2*j] = '.'
				t[floor(self.vertical-i-1-self.slope*j)-1] += 1
		return t

	def pp(self) -> None:
		self.dinv_code()
		n = len(self.DyckWord)
		if n == 0:
			return
		temp = get_intervals(self.DyckWord)
		tt = len(temp)
		a = self.area_sequence()
		k = self.vertical-1
		#x = sum([2*temp[i]-1 for i in range(1,tt-2)[::2]]) + tt//2
		if temp[0] == self.vertical:
			s = ' ' + (2*self.horizontal-1)* '_' + '\n'
		else:
			s = ''.join(self.dinv_box_strs[0]) + (2*temp[-1]-1)* '_' + '\n'
		for i in list(reversed(range(tt-1)))[::2]:
			if i == 0:
				for j in range(temp[0]):
					s += '|' + a[k]* ' x' + '\n'
					k -= 1
				print(s)
				return
			for j in range(temp[i]-1):
				s += ''.join(self.dinv_box_strs[self.vertical-k])+ '|' + a[k]* ' x' + '\n'
				k -= 1
			s += ''.join(self.dinv_box_strs[self.vertical-k]) + (2*temp[i-1]-1)* '_' + '|' + a[k]*' x' + '\n'
			k -= 1

	def pretty_print(self):
		self.pp()

	def area(self):
		return sum(self.area_sequence())

	def dinv(self):
		return sum(self.dinv_code())

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
	if a == 0 or b == 0:
		raise Exception('cannot be zeros')
	if gcd(a,b) != 1:
		raise Exception('coprime')
	R.<q,t> = PolynomialRing(ZZ, 'q,t')
	Cqt = 0
	for i in Dyck_paths(a, b):
		p = RationalDyckPath(i)
		Cqt += q^(p.area()) * t^(p.dinv())
	return Cqt
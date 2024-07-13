class RationalDyckPath:
	def __init__(self, l: list):
		self.DyckWord = l
		self.vertical = sum(l)
		self.horizontal = len(l) - sum(l)
		self.vertices = [(0,0)]
		h = 0
		v = 0
		for i in l:
			if i == 0:
				h += 1
			elif i == 1:
				v += 1
			else:
				raise Exception('input a binary list')
			if h != 0 and v/h < self.vertical/self.horizontal:
				raise Exception('stay above the diagonal')
			self.vertices.append((h,v))

	def is_interior(self, x: tuple):
		l = self.DyckWord
		a = self.vertical
		b = self.horizontal
		if x[0] >= b or x[0] == 0:
			return False
		if x[1] >= a or x[1] == 0:
			return False
		if x[0] != 0 and x[1]/x[0] <= a/b:
			return False
		n = len(self.DyckWord)
		v = 0
		h = 0
		i = 0
		while i < n and v < x[1]:
			v += l[i]
			h += 1-l[i]
			i += 1
		while i < n and h < x[0]:
			if l[i] == 1:
				return True
			h += 1
			i += 1
		return False

	def area(self):
		counter = 0
		for x in range(self.horizontal):
			for y in range(self.vertical):
				if self.is_interior((x,y)):
					#print(x,y)
					counter += 1
		return counter

	def slope_pair(self, x: tuple, y: tuple):
		a = self.vertical
		b = self.horizontal
		if x[0] >= y[0]:
			return False
		if y[0] == x[0] + 1:
			if (y[1] - x[1]) / (y[0] - x[0]) <= a/b:
				return True
			return False
		s1 = (y[1] - x[1]) / (y[0] - x[0])
		s2 = (y[1] - x[1] + 1) / (y[0] - x[0] - 1)
		if s1 <= a/b and a/b <= s2:
			return True
		return False

	def dinv(self):
		i = 0
		l = self.DyckWord
		n = len(l)
		counter = 0
		while i < n:
			if l[i] == 0:
				x = self.vertices[i]
				j = i+1
				while j < n:
					if l[j] == 1:
						y = self.vertices[j]
						if self.slope_pair(x,y):
							#print(x,y)
							counter += 1
					j += 1
			i += 1
		return counter

def paths(a, b, end: tuple):
	if gcd(a,b) != 1:
		raise Exception('coprime')
	if end[0] < 0 or end[1] < 1:
		return []
	if end[0] == 0:
		return [[1]*end[1]]
	if end[1]/end[0] < a/b:
		return []
	if end == (0,1):
		return [[1]]
	if end == (b,a):
		foo = paths(a,b,(b-1,a))
		for i in range(len(foo)):
			foo[i].append(0)
		return foo
	foo = paths(a,b,(end[0]-1,end[1]))
	bar = paths(a,b,(end[0],end[1]-1))
	for i in foo:
		i.append(0)
	for j in bar:
		j.append(1)
	return foo + bar

def Dyck_paths(a: int, b: int):
	if gcd(a,b) != 1:
		raise Exception('coprime')
	return paths(a,b,(b,a))

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
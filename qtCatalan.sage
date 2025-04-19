def get_intervals(binary_list):
    """
    Given a list of 0's and 1's that always starts with 1, this function returns a list
    whose entries correspond to the lengths of consecutive intervals (runs) of 1's and 0's.

    For example:
    >>> get_intervals([1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1])
    [2, 3, 3, 2, 1]
    """
    if not binary_list:
        return []

    intervals = []     # This will store the run lengths.
    current_value = binary_list[0]
    count = 1

    # Start from the second element of the list.
    for item in binary_list[1:]:
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
		self.slope = self.vertical/self.horizontal
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
			if h != 0 and v/h < self.slope:
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
		if x[0] != 0 and x[1]/x[0] <= self.slope:
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

	def area_sequence(self):
		temp = get_intervals(self.DyckWord)
		v = temp[::2]
		v = [sum(v[:i+1]) for i in range(len(v))]
		h = temp[1:][::2]
		h = [sum(h[:i+1]) for i in range(len(h))]
		t = [max(0,floor(y/self.slope)) for y in range(v[0])]
		for i in range(1,len(v)):
			t += [max(0,floor(y/self.slope) - h[i-1]) for y in range(v[i-1],v[i])]
		return t

	def pp(self) -> None:
		n = len(self.DyckWord)
		if n == 0:
			return
		temp = get_intervals(self.DyckWord)
		tt = len(temp)
		a = self.area_sequence()
		k = self.vertical-1
		x = sum([2*temp[i]-1 for i in range(1,tt-2)[::2]]) + tt//2
		s = x *' ' + (2*temp[tt-1]-1)* '_' + '\n'
		for i in list(reversed(range(tt-1)))[::2]:
			x -= 1
			if temp[i] > 1:
				for j in range(temp[i]-1):
					s += x*' '+ '|' + a[k]* ' x' + '\n'
					k -= 1
			x -= 2*temp[i-1]-1
			if i != 0:
				s += x*' ' + (2*temp[i-1]-1)* '_' + '|' + a[k]*' x' + '\n'
			else:
				s += '|'
			k -= 1
		print(s)

	def pretty_print(self):
		self.pp()

	def area(self):
		return sum(self.area_sequence())

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
	#if gcd(a,b) != 1:
	#	raise Exception('coprime')
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
	#if gcd(a,b) != 1:
	#	raise Exception('coprime')
	return [RationalDyckPath(x) for x in paths(a,b,(b,a))]

def qtCatalan(a: int, b: int):
	if a == 0 or b == 0:
		raise Exception('cannot be zeros')
	#if gcd(a,b) != 1:
	#	raise Exception('coprime')
	R.<q,t> = PolynomialRing(ZZ, 'q,t')
	Cqt = 0
	for i in Dyck_paths(a, b):
		p = RationalDyckPath(i)
		Cqt += q^(p.area()) * t^(p.dinv())
	return Cqt
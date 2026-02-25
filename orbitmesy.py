from sage.all import *
from numpy import mean

def A(n: int):
	# A^n = {[i, j]: 1 <= i <= j <= n}
	E = [(i, j) for i in range(1, n + 1) for j in range(i, n + 1)]
	def f(x, y):
		(xi, xj) = x
		(yi, yj) = y
		return yi <= xi and xj <= yj
	return Poset((E, f))

def Exc(pi: Permutation) -> set:
	# Exc(pi) = {[i, pi(i) - 1]: pi(i) > i}
	result = set()
	for i in range(pi.size()):
		if pi[i] > i + 1:
			result.add((i + 1, pi[i] - 1))
	return result

def Exc_inverse(antichain: set, poset) -> Permutation:
	n = len(poset.minimal_elements()) + 1
	pi = [0] * n
	# 1. Map the Excedances: For every interval [i,j]∈A, set π(i)=j+1.
	for (i, j) in antichain:
		pi[i - 1] = j + 1
	# 2. Map the Non-Excedances: Let U={1,…,n}∖{i:[i,j]∈A} be the set of unmapped indices...
	U = set(range(1, n + 1)) - {i for (i, _) in antichain}
	# and V={1,…,n}∖{j+1:[i,j]∈A} be the set of unused values.
	V = set(range(1, n + 1)) - {j + 1 for (_, j) in antichain}
	# Map the elements of U to the elements of V in increasing order.
	for (i, j) in zip(sorted(U), sorted(V)):
		pi[i - 1] = j
	return Permutation(pi)

def rowmotion_for_321(pi: Permutation) -> Permutation:
	# permutation -> antichain (via Exc) -> ideal (all elements <= ones in antichain) -> rowmotion on ideals
	assert pi.avoids([3, 2, 1]), f"{pi} is not 321-avoiding"
	poset = A(pi.size() - 1)
	ideal = set(poset.order_ideal(Exc(pi)))
	antichain = poset.subposet(poset.rowmotion(ideal)).maximal_elements()
	return Exc_inverse(antichain, poset)

class FDSWithOrbitmesy(finite_dynamical_systems.FiniteDynamicalSystem):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self._orbits = self.cycles()
	def orbitmesic_orbits(self, eta):
		averages = [mean([eta(x) for x in O]) for O in self._orbits]
		global_avg = mean(averages)
		return [O for (i, O) in enumerate(self._orbits) if averages[i] == global_avg]

def dyck_path_from_132(pi: Permutation) -> path_tableaux.DyckPath:
	assert pi.avoids([1, 3, 2]), f"{pi} is not 132-avoiding"
	result = [0]
	for i in range(pi.size()):
		right_greater = 0
		for j in pi[i + 1 : ]:
			if j > pi[i]:
				right_greater += 1
		while result[-1] <= right_greater:
			result.append(result[-1] + 1)
		result.append(result[-1] - 1)
	return path_tableaux.DyckPath(result)

def dyck_path_to_321(d: path_tableaux.DyckPath) -> Permutation:
	peaks = {}
	for i in range(1, len(d) - 1):
		if d[i - 1] + 1 == d[i] == d[i + 1] + 1:
			peaks[i] = d[i]
	pi = [0 for _ in range((len(d) - 1) // 2)]
	for (x, y) in peaks.items():
		pi[(x - y) // 2] = (x + y) // 2
	unused_x = sorted([x for x in range(len(pi)) if pi[x] == 0])
	unused_y = sorted([y for y in range(1, len(pi) + 1) if y not in pi])
	for (x, y) in zip(unused_x, unused_y):
		pi[x] = y
	assert Permutation(pi).avoids([3, 2, 1]), f"{pi} is not 321-avoiding (something is wrong)"
	return Permutation(pi)

def dyck_path_from_321(pi: Permutation) -> path_tableaux.DyckPath:
	assert pi.avoids([3, 2, 1]), f"{pi} is not 321-avoiding"
	excedances = {}
	for i in range(pi.size()):
		if pi[i] >= i + 1:
			excedances[i] = pi[i] # get upper left corner on square
	word = []
	(x, y) = (0, 0)
	for (ex, ey) in excedances.items():
		while x < ex:
			word.append(0)
			x += 1
		while y < ey:
			word.append(1)
			y += 1
		word.append(0)
		x += 1
	while x < pi.size():
		word.append(0)
		x += 1
	return path_tableaux.DyckPath(DyckWord(word))

def dyck_path_to_132(d: path_tableaux.DyckPath) -> Permutation:
	right_greater = []
	y = 0
	for step in d.to_DyckWord():
		if step == 1:
			y += 1
		else:
			y -= 1
			right_greater.append(y)
	pi = []
	unused = list(range(1, (len(d) + 1) // 2))
	for x in right_greater:
		pi.append(unused.pop(-x - 1))
	assert Permutation(pi).avoids([1, 3, 2]), f"{pi} is not 132-avoiding (something is wrong)"
	return Permutation(pi)

def rowmotion_for_132(pi: Permutation) -> Permutation:
	# 132 -> dyck path -> 321 -> 321 via rowmotion -> dyck path -> 132
	return dyck_path_to_132(dyck_path_from_321(rowmotion_for_321(dyck_path_to_321(dyck_path_from_132(pi)))))

n = 1
while True:
	fds = FDSWithOrbitmesy(Permutations(n, avoiding=[1, 3, 2]), rowmotion_for_132)
	print(f"{n}: {fds.orbitmesic_orbits(lambda x: x.number_of_fixed_points())}")
	n += 1
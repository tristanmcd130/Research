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

def rowmotion_for_321_avoiding(pi: Permutation) -> Permutation:
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

n = 4
fds = FDSWithOrbitmesy(Permutations(n, avoiding=[3, 2, 1]), rowmotion_for_321_avoiding)
print(f"Orbitmesic orbits of S_{n}(321) wrt number of descents: {fds.orbitmesic_orbits(lambda x: x.number_of_descents())}")
print(f"Orbitmesic orbits of S_{n}(321) wrt number of fixed points: {fds.orbitmesic_orbits(lambda x: x.number_of_fixed_points())}")
print(f"Orbitmesic orbits of S_{n}(321) wrt number of inversions: {fds.orbitmesic_orbits(lambda x: x.number_of_inversions())}")

# try implementing phi (binary tree -> noncrossing permutation) from https://arxiv.org/pdf/2401.14318
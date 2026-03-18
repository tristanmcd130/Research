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

# 321 avoiding -> tableaux via rsk -> another tableaux via K_B. after that keep doing K_1 or K_2 until you can do K_B again. find a way to go directly from the permutation represented by the very 1st tableaux to the last one
# rsk accepts permutations directly

def K1(pi: Permutation) -> set[Permutation]:
	# If x < y < z and pi = x1...yxz..xn, then K1(pi) = x1..yzx...xn
	results = set()
	for i in range(len(pi) - 2):
		y = pi[i]
		x = pi[i + 1]
		z = pi[i + 2]
		if x < y < z:
			pi_copy = list(pi)
			pi_copy[i + 1], pi_copy[i + 2] = pi_copy[i + 2], pi_copy[i + 1]
			results.add(Permutation(pi_copy))
	return results

def K2(pi: Permutation) -> set[Permutation]:
	# If x < y < z and pi = x1...xzy..xn, then K1(pi) = x1..zxy...xn
	results = set()
	for i in range(len(pi) - 2):
		x = pi[i]
		z = pi[i + 1]
		y = pi[i + 2]
		if x < y < z:
			pi_copy = list(pi)
			pi_copy[i], pi_copy[i + 1] = pi_copy[i + 1], pi_copy[i]
			results.add(Permutation(pi_copy))
	return results

def KB(pi: Permutation) -> set[Permutation]:
	# KB(pi) = K1(pi) ∩ K2(pi)
	return K1(pi) & K2(pi)

def increasing_subseqs(elements: list, length: int) -> list[list]:
	if len(elements) < length:
		return []
	if length == 0:
		return [[]]
	results = []
	for e in elements:
		for s in increasing_subseqs([x for x in elements if x > e], length - 1):
			results.append([e] + s)
	return results

def row_index(rows: list[list], x) -> int:
	for i in range(len(rows)):
		if x in rows[i]:
			return i
	raise ValueError(f"{x} not in {rows}")

def abc(Q: Tableau) -> bool:
	"""
	A tableau has an abc pattern if there are 3 numbers x, y, and z such that:
	1. x < y < z
	2. x is at the bottom of its column
	3. z is to the right of x (but not necessarily in the same row)
	4. y is in the same column as z or to the right of it
	"""
	Q_T = [[row[i] for row in Q if i < len(row)] for i in range(max(Q.shape()))] # "transpose" of Q to help see if x is at the bottom of a column
	for x, y, z, in increasing_subseqs(range(1, Q.size() + 1), 3):
		x_row = row_index(Q_T, x)
		y_row = row_index(Q_T, y)
		z_row = row_index(Q_T, z)
		if Q_T[x_row][-1] == x and z_row > x_row and y_row >= z_row:
			return True
	return False

def abcd(Q: Tableau) -> bool:
	"""
	A tableau has an abcd pattern if there are 4 numbers x, y, z, and w such that:
	1. x < y < z < w
	2. x is not at the bottom of its column
	3. w is right below x
	4. z is somewhere to the right of x
	5. y is in the same column as z or to the right of it
	"""
	Q_T = [[row[i] for row in Q if i < len(row)] for i in range(max(Q.shape()))]
	for x, y, z, w in increasing_subseqs(range(1, Q.size() + 1), 4):
		x_row = row_index(Q_T, x)
		y_row = row_index(Q_T, y)
		z_row = row_index(Q_T, z)
		if Q_T[x_row][-1] != x and Q_T[x_row][Q_T[x_row].index(x) + 1] == w and z_row > x_row and y_row >= z_row:
			return True
	return False

def plot_connections(n):
	S_n = list(Permutations(n))
	K1_edges = []
	K2_edges = []
	KB_edges = []
	for pi in S_n:
		kb = KB(pi)
		for x in K1(pi) - kb:
			K1_edges.append((x, pi, "K1"))
		for x in K2(pi) - kb:
			K2_edges.append((x, pi, "K2"))
		for x in kb:
			KB_edges.append((x, pi, "KB"))

	g = DiGraph([S_n, K1_edges + K2_edges + KB_edges], multiedges=True)
	g.plot(
		vertex_colors={
			"yellow": [pi for pi in S_n if not abc(RSK(pi)[1]) and not abcd(RSK(pi)[1])],
		},
		edge_colors={"red": K1_edges, "blue": K2_edges, "green": KB_edges},
		layout="graphviz",
		figsize=[48, 48]
	).save("graph.png")

plot_connections(5)

# SD isnt in sage, its defined in the GHOSS25motzkin paper

"""
RSK is a bijection between S_n and 2 tableaux P and Q.
Thm: 2 perms have the same P tableau iff they're related by some Knuth moves.
In that image I made, each component is called a Knuth graph (all connected nodes there share tableau P with shape lambda). In that graph, each node corresponds with a standard tableau (increasing rows and cols) with shape lambda.
Superstandard tableau: read it top to bottom, left to right as 123456789...
1 3 5
2 4
Lemma: If w goes to v via K1 or KB, swapping numbers in positions i and i + 1, then Q(w) -> Q(v) involves swapping i and i + 1 in the tableaux.
Another lemma: A superstandard tableau's graph only involves K1 and KB arrows (K2 is possible, but never on its own; you can always also do K1 and thus KB).
1st row of P will always be the same as the 1st row of the soliton decomposition.
These are equivalent:
- P(w) = SD(w)
- shape of P(w) = shape of SD(w)
- SD(w) is standard
In the dominance partial order(?), SD(w) <= P(w).
If shape of P(w) is:
x x x x
x
x
x
then P(w) is good, otherwise it's bad.

Try expanding your code for all of S_5 and not just 321-avoiding (have to add abcd condition to good tableau check), and then swap the directions of all edges.
"""
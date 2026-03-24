from sage.all import *

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

def superstandard(shape: tuple[int]) -> Tableau:
	result = []
	for n in shape:
		result.append([0] * n)
	n = 1
	for i in range(max(shape)):
		try:
			j = 0
			while True:
				result[j][i] = n
				n += 1
				j += 1
		except IndexError:
			pass
	return Tableau(result)

# somehow make it so it only generates permutations with a given p tableau, or a superstandard p tableau for a given shape
# allow a specific permutation as an input, then graph only the component with that in it
# print lists of good and bad?
def plot_connections(n, *, P=None, shape=None):
	assert P is None or shape is None, f"P and shape cannot both be specified"
	assert sum(shape) == n, f"Sum of shape must be n"
	S_n = list(filter(lambda pi: (P is None or RSK(pi)[0] == P) and (shape is None or RSK(pi)[0] == superstandard(shape)), Permutations(n)))
	K1_edges = []
	K2_edges = []
	KB_edges = []
	for pi in S_n:
		kb = KB(pi)
		for x in K1(pi) - kb:
			K1_edges.append((pi, x, "K1"))
		for x in K2(pi) - kb:
			K2_edges.append((pi, x, "K2"))
		for x in kb:
			KB_edges.append((pi, x, "KB"))

	g = DiGraph([S_n, K1_edges + K2_edges + KB_edges], multiedges=True)
	g.plot(
		vertex_colors={
			"yellow": [pi for pi in S_n if not abc(RSK(pi)[1]) and not abcd(RSK(pi)[1])],
		},
		edge_colors={"red": K1_edges, "blue": K2_edges, "green": KB_edges},
		layout="graphviz"
	).save("graph.png")

plot_connections(7, shape=(4, 2, 1))

"""
SD isnt in sage, its defined in the GHOSS25motzkin paper
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
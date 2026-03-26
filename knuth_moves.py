from sage.all import *
from matplotlib import colormaps as cm
import matplotlib.colors as mcolors
from collections import defaultdict

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

def abc(Q: Tableau) -> int:
	"""
	A tableau has an abc pattern if there are 3 numbers x, y, and z such that:
	1. x < y < z
	2. x is at the bottom of its column
	3. z is to the right of x (but not necessarily in the same row)
	4. y is in the same column as z or to the right of it
	"""
	count = 0
	Q_T = [[row[i] for row in Q if i < len(row)] for i in range(max(Q.shape()))] # "transpose" of Q to help see if x is at the bottom of a column
	for x, y, z, in increasing_subseqs(range(1, Q.size() + 1), 3):
		x_row = row_index(Q_T, x)
		y_row = row_index(Q_T, y)
		z_row = row_index(Q_T, z)
		if Q_T[x_row][-1] == x and z_row > x_row and y_row >= z_row:
			count += 1
	return count

def abcd(Q: Tableau) -> int:
	"""
	A tableau has an abcd pattern if there are 4 numbers x, y, z, and w such that:
	1. x < y < z < w
	2. x is not at the bottom of its column
	3. w is right below x
	4. z is somewhere to the right of x
	5. y is in the same column as z or to the right of it
	"""
	count = 0
	Q_T = [[row[i] for row in Q if i < len(row)] for i in range(max(Q.shape()))]
	for x, y, z, w in increasing_subseqs(range(1, Q.size() + 1), 4):
		x_row = row_index(Q_T, x)
		y_row = row_index(Q_T, y)
		z_row = row_index(Q_T, z)
		if Q_T[x_row][-1] != x and Q_T[x_row][Q_T[x_row].index(x) + 1] == w and z_row > x_row and y_row >= z_row:
			count += 1
	return count

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

def plot_connections(*,
					 n: int | None = None,
					 P: Permutation | None = None,
					 shape: tuple[int, ...] | None = None,
					 target_pi: Permutation | None = None,
					 print_good: bool = False,
					 figsize: int | None = 12,
					 filename: str | None = "graph.png") -> None:
	"""
	Provide exactly 1 of the following parameters (which must be provided as keyword arguments):
	- n: an integer. All permutations in S_n will be graphed.
	- shape: a tuple of integers. All permutations in S_len(shape) whose P tableau has this shape will be graphed.
	- target_pi: a Permutation. Only this permutation, and all the other ones in S_n related to it via Knuth moves, will be graphed.

	Optionally provide any number of these parameters:
	- print_good: a boolean. If True, it will print all the good and bad permutations in the graph, how many good and bad permutations were graphed, and the "badness" of each bad permutation (number of abc patterns + abcd patterns).
	- figsize: an integer. Making this bigger will increase the size of the graph.
	- filename: a string. The name of the file to which the graph will be saved. You can use extensions other than .png, such as .pdf.
	"""
	assert [n, P, shape, target_pi].count(None) == 3, f"Exactly 1 of n, P, shape, or target_pi must be provided"
	if n is None:
		if P is not None:
			n = len(P)
		elif shape is not None:
			n = sum(shape)
		elif target_pi is not None:
			n = len(target_pi)
	S_n = list(filter(lambda pi: (P is None or RSK(pi)[0] == P) and (shape is None or RSK(pi)[0] == superstandard(shape)), Permutations(n)))
	assert target_pi is None or target_pi in S_n, f"{target_pi} is not in {S_n}"

	pi_graph = set() if target_pi is None else {target_pi}
	K1_edges = []
	K2_edges = []
	KB_edges = []
	for pi1 in S_n:
		kb = KB(pi1)
		for pi2 in K1(pi1) - kb:
			if target_pi is None or pi1 in pi_graph or pi2 in pi_graph:
				K1_edges.append((pi1, pi2, "K1"))
				pi_graph |= {pi1, pi2}
		for pi2 in K2(pi1) - kb:
			if target_pi is None or pi1 in pi_graph or pi2 in pi_graph:
				K2_edges.append((pi1, pi2, "K2"))
				pi_graph |= {pi1, pi2}
		for pi2 in kb:
			if target_pi is None or pi1 in pi_graph or pi2 in pi_graph:
				KB_edges.append((pi1, pi2, "KB"))
				pi_graph |= {pi1, pi2}

	S_n = list(pi_graph)
	badness = {pi: abc(RSK(pi)[1]) + abcd(RSK(pi)[1]) for pi in S_n}
	cmap = cm.get_cmap("viridis")
	norm = mcolors.Normalize(vmin=0, vmax=max(badness.values()))
	vertex_colors = defaultdict(list)
	for pi in S_n:
		vertex_colors[mcolors.to_hex(cmap(norm(badness[pi])))].append(pi)
	g = DiGraph([S_n, K1_edges + K2_edges + KB_edges], multiedges=True)
	g.plot(
		vertex_colors=vertex_colors,
		edge_colors={"red": K1_edges, "blue": K2_edges, "green": KB_edges},
		figsize=(figsize, figsize),
		layout="graphviz"
	).save(filename)
	if print_good:
		good = [k for (k, v) in badness.items() if v == 0]
		print(f"Good permutations ({len(good)}):")
		for pi in good:
			print(pi)
		bad = list(set(S_n) - set(good))
		print(f"\nBad permutations ({len(bad)}):")
		for pi in bad:
			print(f"{pi} (badness: {badness[pi]})")

for [x, y] in Partitions(7, length=2):
	plot_connections(shape=(x, y), print_good=True, figsize=6, filename=f"{x}_{y}.png")

"""
3/24:
Number of descents (pi(i) > pi(i + 1)) can tell you if its bad: 
Docker env? Make it somehow compatible with Jupyter
Count how many times abc and abcd patterns appear and see how that affects how far away a bad is from a cluster of goods in a graph

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
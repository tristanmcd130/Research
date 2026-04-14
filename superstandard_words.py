from sage.all import *
from knuth_moves import *
from collections import Counter
import matplotlib.pyplot as plt

def superstandard_word(Q: Tableau) -> Permutation:
	return Permutation(RSK_inverse(superstandard(Q.shape()), Q)[1])

if __name__ == "__main__":
	counter = Counter()
	for Q in StandardTableaux(10):
		w = superstandard_word(Q)
		abc_patterns = abc(Q)
		abcd_patterns = abcd(Q)
		abcd_patterns_in_1st_column = len([p[0] in tableau_transpose(Q)[0] for p in abc_patterns | abcd_patterns])
		height_diff_from_sd = abs(Q.height() - (w.number_of_descents() + 1))
		counter[(abcd_patterns_in_1st_column, height_diff_from_sd)] += 1
		# if len(abc_patterns | abcd_patterns) > 0:
		# 	print(f"{Q} has {abcd_patterns_in_1st_column} abc(d) patterns involving the 1st column, and a height difference of {height_diff_from_sd} from SD(w)")
		# 	if abcd_patterns_in_1st_column != height_diff_from_sd:
		# 		print(f"abc patterns: {abc_patterns}")
		# 		print(f"abcd patterns: {abcd_patterns}")
	max_x = max([x for (x, _) in counter.keys()])
	max_y = max([y for (_, y) in counter.keys()])
	data = [[0 for _ in range(max_x + 1)] for _ in range(max_y + 1)]
	fig, ax = plt.subplots()
	for ((x, y), v) in counter.items():
		data[y][x] = v
		ax.text(x, y, str(v), ha="center", va="center", color="w")
	im = ax.imshow(data)
	ax.set_xlabel("Number of abc(d) patterns in 1st column")
	ax.set_ylabel("Height difference from SD(w)")
	fig.tight_layout()
	plt.show()

"""
4/14:
(dominance partial order: if (p_1, p_2... p_m) and (q_1, q_2... q_n) are shapes of a tableau, then p <= q iff for all k >= 1, sum from i = 1 to k of p_i <= sum from i = 1 to k of q_i)
proving descents(w) + 1 == height(Q) => 1st column of Q has no abc(d) patterns:
assume contrapositive, 1st column has abc(d) patterns => descents(w) + 1 != height(Q)
how does number of abcd patterns influence number of descents?
maybe difference in height between SD(w) and Q corresponds to number of abc(d) patterns involving 1st row?

4/9:
alternate way of seeing if a tableau is bad: in motzkin paper, you do something to turn columns of tableaus into partitions, then if theres a descent anywhere, its bad
how are 2 bad tableaux with height equal to SD(w), where 1 becomes "bad" (has a descent) in u_2 vs in u_1, different?
for the height of a tableau Q to be the same as the height of SD(w) (Q.height() == w.number_of_descents() + 1), maybe it needs to not have an abc(d) pattern where x (and w) are in the 1st column? check other direction of implication too
"""
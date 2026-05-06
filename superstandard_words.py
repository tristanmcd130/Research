from sage.all import *
from knuth_moves import *
from bbs import *
from collections import Counter
import matplotlib.pyplot as plt

def superstandard_word(Q: Tableau) -> Permutation:
	return Permutation(RSK_inverse(superstandard(Q.shape()), Q)[1])

def SD_tableau(Q: Tableau) -> Tableau:
	return SD(superstandard_word(Q))

def cut(Q: Tableau, n: int) -> Tableau:
	return Tableau([row[n - 1 : ] for row in Q if len(row) - n + 1 > 0])

if __name__ == "__main__":
	counter = Counter()
	for Q in StandardTableaux(10):
		w = superstandard_word(Q)
		abc_patterns = abc(Q)
		abcd_patterns = abcd(Q)
		bad_patterns = len([p for p in abc_patterns | abcd_patterns if p[0] in tableau_transpose(Q)[0]])
		height_diff = abs(SD_tableau(Q).height() + cut(SD_tableau(Q), 2).height() - Q.height() - cut(Q, 2).height())
		counter[(bad_patterns, height_diff)] += 1
	max_x = max([x for (x, _) in counter.keys()])
	max_y = max([y for (_, y) in counter.keys()])
	data = [[0 for _ in range(max_x + 1)] for _ in range(max_y + 1)]
	fig, ax = plt.subplots()
	for ((x, y), v) in counter.items():
		data[y][x] = v
		ax.text(x, y, str(v), ha="center", va="center", color="w")
	im = ax.imshow(data)
	ax.set_xlabel("Number of bad patterns")
	ax.set_ylabel("Height difference between 1st 2 columns of SD(w) and Q")
	ax.set_xticks(range(max_x + 1))
	ax.set_yticks(range(max_y + 1))
	fig.tight_layout()
	plt.show()

"""
4/30:
lemma: k is a descent in u2 iff k is a descent in standardize(cut(Q, 2))
conjecture: height(2nd col of SD(w)) = height(2nd col of Q) iff there are no abc(d) patterns with a in 2nd col of Q
make a function in bbs.py that takes a tableau, converts it to superstandard word, then passes it through SD

4/16:
write bbs code
try to find "equivalence classes" between abc(d) patterns such that the number of classes present in a tableau = difference in height between Q and SD(w)

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
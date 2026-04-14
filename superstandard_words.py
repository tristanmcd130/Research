from sage.all import *
from knuth_moves import *

def superstandard_word(Q: Tableau) -> Permutation:
	return Permutation(RSK_inverse(superstandard(Q.shape()), Q)[1])

if __name__ == "__main__":
	for n in range(4, 11):
		sd_height = 0
		no_1st_column_abcd = 0
		sd_height_implies_no_1st_column_abcd = 0
		no_1st_column_abcd_implies_sd_height = 0
		bad = 0
		for Q in StandardTableaux(n):
			w = superstandard_word(Q)
			abc_patterns = abc(Q)
			abcd_patterns = abcd(Q)
			if len(abc_patterns | abcd_patterns) > 0:
				# print(f"{Q} is bad")
				bad += 1
				same_height_as_sd = Q.height() == w.number_of_descents() + 1
				no_abcd_patterns_in_1st_column = not any(map(lambda x: x[0] in tableau_transpose(Q)[0], abc_patterns | abcd_patterns))
				if same_height_as_sd:
					sd_height += 1
				if no_abcd_patterns_in_1st_column:
					no_1st_column_abcd += 1
				if not same_height_as_sd or no_abcd_patterns_in_1st_column:
					# print(f"Same height as SD(w) => no abcd patterns with x in 1st column")
					sd_height_implies_no_1st_column_abcd += 1
				if not no_abcd_patterns_in_1st_column or same_height_as_sd:
					# print(f"No abcd patterns with x in 1st column => same height as SD(w)")
					no_1st_column_abcd_implies_sd_height += 1
				print(f"{Q} has {len(abcd_patterns)} abcd patterns and {w.number_of_descents()} descents")
		print(f"n = {n}:")
		print(f"# of tableau with same height as SD(w): {sd_height}")
		print(f"# of tableau with no abcd patterns with x in 1st column: {no_1st_column_abcd}")
		print(f"# of tableau where same height as SD(w) => no abcd patterns with x in 1st column: {sd_height_implies_no_1st_column_abcd}")
		print(f"# of tableau where no abcd patterns with x in 1st column => same height as SD(w): {no_1st_column_abcd_implies_sd_height}")
		print(f"# of bad tableau: {bad}\n")

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
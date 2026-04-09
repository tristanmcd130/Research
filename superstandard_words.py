from sage.all import *
from knuth_moves import *

def superstandard_word(Q: Tableau) -> Permutation:
	return Permutation(RSK_inverse(superstandard(Q.shape()), Q)[1])

if __name__ == "__main__":
	# no_abcd = set()
	# has_abcd = set()
	# same_height = set()
	# diff_height = set()
	for Q in StandardTableaux(10):
		w = superstandard_word(Q)
		# if abcd(Q) == 0:
		# 	no_abcd.add(w)
		# else:
		# 	has_abcd.add(w)
		# if Q.height() == len(w.descents()) + 1:
		# 	same_height.add(w)
		# else:
		# 	diff_height.add(w)
		if Q.height() == w.number_of_descents() + 1:
			print(f"{Q} has height {'equal to' if Q.height() == w.number_of_descents() + 1 else 'different from'} the height of SD(w), {abcd(Q)} abcd patterns, and its superstandard word is {w}")
	# for x in [no_abcd, has_abcd]:
	# 	for y in [same_height, diff_height]:
	# 		print(len(x & y))

"""
4/9:
alternate way of seeing if a tableau is bad: in motzkin paper, you do something to turn columns of tableaus into partitions, then if theres a descent anywhere, its bad
how are 2 bad tableaux with height equal to SD(w), where 1 becomes "bad" (has a descent) in u_2 vs in u_1, different?
for the height of a tableau Q to be the same as the height of SD(w) (Q.height() == w.number_of_descents() + 1), maybe it needs to not have an abc(d) pattern where x (and w) are in the 1st column? check other direction of implication too
"""
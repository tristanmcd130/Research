"""
Consider all Q tableau where w = RSK^-1(P, Q) and P is superstandard. For n >= 7, how many of those tableau have a 1st column with height equal to the height of the 1st column of SD(w)?
2 facts: height of the 1st column of Q = length of the longest decreasing subsequence of the permutation, and height of the 1st column of SD(w) = number of descents in w + 1
For Q tableau where height(SD(w)) = height(Q), does it always avoid abcd patterns?
How to generate these Q tableau: for each Q tableau of size n (Tableaux(n)), do w = RSK^-1(superstandard(Q.shape()), Q), then see if len(w.descents()) == height(Q)
"""

from sage.all import *
from knuth_moves import *

def equal_heights(n: int) -> int:
	abcd_count = 0
	count = 0
	total = 0
	for Q in StandardTableaux(n):
		if len(abc(Q) | abcd(Q)) > 0:
			w = Permutation(RSK_inverse(superstandard(tuple(Q.shape())), Q)[1])
			if Q.height() == len(w.descents()) + 1:
				# print(f"{Q} has height equal to height(SD(w))", end="")
				count += 1
				if len(abcd(Q)) > 0:
					abcd_count += 1
					print(f"{Q} has abcd patterns")
				# print()
			total += 1
	print(f"{total} bad tableaux of size {n}, {count} of which have height equal to height(SD(w)), and {abcd_count} of those have abcd patterns")
	return count

if __name__ == "__main__":
	equal_heights(11)
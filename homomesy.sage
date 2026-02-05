from sage.all import *

def A(n: int):
	# A^n = {[i, j]: 1 <= i <= j <= n}
	E = [(i, j) for i in range(1, n + 1) for j in range(i, n + 1)]
	def f(x, y):
		(xi, xj) = x
		(yi, yj) = y
		return yi <= xi and xj <= yj
	return Poset((E, f))

# xxx-avoiding permutation bijections exist as DyckWord.to_xxx_avoiding_permutation()

def E_p(pi: Permutation) -> list[list[int]]:
	n = pi.size()
	array = [[0 for _ in range(n)] for _ in range(n)]
	for i in range(1, n):
		for j in range(n):
			if pi[j] == i - 1:
				

# permutation -> array
# array -> dyck path (peaks at excedances)
# dyck path -> antichain
# antichain -> antichain (rowmotion)

# try doing rowmotion directly on dyck paths
# 312/132-avoiding have a bijection with a binary search tree, try doing rowmotion on them
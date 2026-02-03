from sage.all import *

def A(n: int):
	# A^n = {[i, j]: 1 <= i <= j <= n}
	E = [(i, j) for i in range(1, n + 1) for j in range(i, n + 1)]
	def f(x, y):
		(xi, xj) = x
		(yi, yj) = y
		return yi <= xi and xj <= yj
	return Poset((E, f))

def permutation_to_array(pi: Permutation):
	

# permutation -> array
# array -> dyck path (peaks at excedances)
# dyck path -> antichain
# antichain -> antichain (rowmotion)
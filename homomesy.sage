from sage.all import *

def A(n: int):
	# A^n = {[i, j]: 1 <= i <= j <= n}
	E = [(i, j) for i in range(1, n + 1) for j in range(i, n + 1)]
	def f(x, y):
		(xi, xj) = x
		(yi, yj) = y
		return yi <= xi and xj <= yj
	return Poset((E, f))

def Exc(pi: Permutation):
	result = set()
	for i in range(pi.size()):
		print(f"pi({i + 1}) = {pi[i]}")
		if pi[i] > i + 1:
			result.add((i + 1, pi[i] - 1))
	return result
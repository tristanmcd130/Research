from sage.all import *

def A(n):
	# A^n = {[i, j]: 1 <= i <= j <= n}
	E = set()
	R = set()
	def interval_cover(i, j):
		if i < j:
			E.add((i, j))
			interval_cover(i, j - 1)
			R.add(((i, j - 1), (i, j)))
			interval_cover(i + 1, j)
			R.add(((i + 1, j), (i, j)))
		elif i == j:
			E.add((i, j))
	interval_cover(1, n)
	return Poset((E, R))
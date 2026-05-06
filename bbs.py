from sage.all import *

def SD(pi: Permutation) -> Tableau:
	# print(f"SD BEGIN: {pi}")
	boxes = list(pi)
	for _ in range(len(pi) - 2): # technically just a conjecture
		for n in range(1, len(pi) + 1):
			i = [k for (k, v) in enumerate(boxes) if v == n][0]
			new_i = i
			while boxes[new_i] != 0:
				new_i += 1
				if new_i >= len(boxes):
					boxes.append(0)
			boxes[new_i] = boxes[i]
			boxes[i] = 0
			while boxes[0] == 0:
				boxes.pop(0)
			# print(boxes)
	solitons = []
	while len(boxes) > 0:
		solitons.append([])
		while len(boxes) > 0 and boxes[0] != 0:
			solitons[-1].append(boxes.pop(0))
		while len(boxes) > 0 and boxes[0] == 0:
			boxes.pop(0)
	# print(f"SD END: {sorted(solitons, key=len, reverse=True)}")
	return Tableau(sorted(solitons, key=len, reverse=True)) # SD isn't guaranteed to be standard, but we need to at least make sure it's rows get shorter

# The version of SD that takes a Tableau is in superstandard_words.py.
# It has to be there because otherwise these files would mutually import each other, which isn't allowed.

if __name__ == "__main__":
	print(SD(Permutation([4, 5, 2, 3, 6, 1])))
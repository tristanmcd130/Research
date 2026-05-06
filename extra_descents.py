from sage.all import *
from superstandard_words import *

# We proved that no abc(d) in 1st column <=> height(Q) = height(SD(w))
# The following is always true: height(SD(w)) = descents(w) + 1 = descents(Q) + 1
# When height(Q) != height(SD(w)), we posit the difference is equal to the number of "equivalence classes" between 1st column patterns
# What this shows is that 2 patterns are "equivalent" if c = b + 1:

if __name__ == "__main__":
	counter = Counter()
	for Q in StandardTableaux(10):
		patterns = abc(Q) | abcd(Q)
		first_column_patterns = [p for p in patterns if p[0] in tableau_transpose(Q)[0]]
		c_succeeds_b = [p for p in first_column_patterns if p[2] == p[1] + 1]
		print(f"{Q}\nabc(d) patterns in 1st column: {first_column_patterns}\nHeight difference with SD(w): {superstandard_word(Q).number_of_descents() + 1 - Q.height()}\n")
		# assert len(c_succeeds_b) == abs(Q.height() - superstandard_word(Q).number_of_descents() - 1), "Conjecture is wrong"
		counter[(len(c_succeeds_b), abs(Q.height() - superstandard_word(Q).number_of_descents() - 1))] += 1
	max_x = max([x for (x, _) in counter.keys()])
	max_y = max([y for (_, y) in counter.keys()])
	data = [[0 for _ in range(max_x + 1)] for _ in range(max_y + 1)]
	fig, ax = plt.subplots()
	for ((x, y), v) in counter.items():
		data[y][x] = v
		ax.text(x, y, str(v), ha="center", va="center", color="w")
	im = ax.imshow(data)
	ax.set_xlabel("Number of abc(d) patterns in 1st column where c = b + 1")
	ax.set_ylabel("Height difference from SD(w)")
	ax.set_xticks(range(max_x + 1))
	ax.set_yticks(range(max_y + 1))
	fig.tight_layout()
	plt.show()
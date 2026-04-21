# Proof of the Conjecture

## 1. Definitions and Theorem
Let $`Q`$ be a standard Young tableau.
*   **$`h(Q)`$**: The height of the 1st column of $`Q`$.
*   **$`c(Q)`$**: The number of descents in $`Q`$. An index $`i`$ is a descent if $`i+1`$ appears in a strictly lower row than $`i`$ (i.e., $`R(i+1) > R(i)`$).

We analyze two specific forbidden patterns:
*   **`abc` pattern**: Three values $`x < y < z`$ such that $`x`$ is at the bottom of its column, and $`C(x) < C(z) \le C(y)`$.
*   **`abcd` pattern**: Four values $`x < y < z < w`$ such that $`w`$ is directly below $`x`$ in the same column ($`C(x) = C(w)`$), and $`C(x) < C(z) \le C(y)`$.

**Theorem**: $`c(Q) = h(Q) - 1`$ if and only if $`Q`$ contains no `abc` or `abcd` patterns where the smallest element $`x`$ is in the 1st column.

## 2. The Core Lemma of Standard Tableaux
To prove this, we rely on a fundamental property of standard tableaux:

**Core Lemma**: For any standard Young tableau and any entry $`y`$, let $`z = y+1`$. Exactly one of the following is true:
1.  **Ascent**: $`z`$ is placed strictly to the right and weakly above $`y`$ ($`C(z) > C(y)`$ and $`R(z) \le R(y)`$).
2.  **Descent**: $`z`$ is placed weakly to the left and strictly below $`y`$ ($`C(z) \le C(y)`$ and $`R(z) > R(y)`$).

**Proof of Core Lemma**:
Since $`y`$ and $`z`$ are consecutive integers, they correspond to adding two cells to a valid Young diagram one after the other.
Can $`z`$ be placed strictly below and strictly to the right of $`y`$? That is, can we have both $`R(z) > R(y)`$ and $`C(z) > C(y)`$?
For the cell of $`z`$ at $`(R(z), C(z))`$ to be a valid addition to the Young diagram, the cells immediately above it at $`(R(z)-1, C(z))`$ and immediately to its left at $`(R(z), C(z)-1)`$ must already be in the diagram.
Since $`y`$ is added immediately before $`z`$, these two cells must be present in the diagram *before* $`y`$ is added, or one of them must be $`y`$ itself.
*   If one of them is $`y`$, then either $`R(y) = R(z)-1`$ and $`C(y) = C(z)`$ (contradicting $`C(z) > C(y)`$) or $`R(y) = R(z)`$ and $`C(y) = C(z)-1`$ (contradicting $`R(z) > R(y)`$).
*   Therefore, both cells must be present before $`y`$ is added. By the shape properties of a Young diagram, this implies the entire rectangle of cells $`(r,c)`$ with $`r \le R(z)-1`$ and $`c \le C(z)-1`$ must be present before $`y`$ is added.
Since we assumed $`R(y) \le R(z)-1`$ and $`C(y) \le C(z)-1`$, the cell $`(R(y), C(y))`$ would already be occupied, meaning $`y`$ cannot be placed there. This is a contradiction.
Thus, we cannot have both $`R(z) > R(y)`$ and $`C(z) > C(y)`$. Since $`y`$ and $`z`$ cannot occupy the same cell, and $`R(z) \le R(y)`$ with $`C(z) \le C(y)`$ is impossible (the cell would already be occupied), we are left with exactly the two cases described: an ascent or a descent. $`\blacksquare`$

## 3. Mandatory Descents
Let the elements in the 1st column of $`Q`$ be $`v_0 < v_1 < \dots < v_{h-1}`$. 
Because $`v_k`$ is in the 1st column, $`C(v_k) = 0`$. 

Consider the element $`v_k - 1`$ for any $`k \in \{1, \dots, h-1\}`$. Because $`v_k > v_{k-1} \ge 1`$, this element exists.
By the Core Lemma, the transition from $`v_k - 1`$ to $`v_k`$ is either an ascent or a descent.
If it were an ascent, we would have $`C(v_k) > C(v_k - 1)`$. But $`C(v_k) = 0`$, and column indices cannot be negative. This is a contradiction.
Therefore, the transition from $`v_k - 1`$ to $`v_k`$ MUST be a descent.

There are exactly $`h(Q) - 1`$ such mandatory descents (one entering each cell of the first column, except the very first cell $`v_0 = 1`$). Therefore, $`c(Q) \ge h(Q) - 1`$ is always true.
Having $`c(Q) = h(Q) - 1`$ means that these mandatory descents are the *only* descents in the entire tableau. Any other descent is an "extra" descent.

## 4. Extra Descents Implies Forbidden Patterns
Suppose $`c(Q) > h(Q) - 1`$. Then there is an extra descent $`y \to z`$ (where $`z = y+1`$) that is not one of the mandatory descents entering the 1st column.
Since $`y`$ is a descent, the Core Lemma dictates $`C(z) \le C(y)`$. 
Because $`y \to z`$ is not a mandatory descent, $`z`$ cannot be in the 1st column, so $`C(z) > 0`$. Thus, $`0 < C(z) \le C(y)`$.
We map this extra descent to our forbidden patterns:

**Case A: The extra descent occurs after the 1st column is finished**
This means $`y > v_{h-1}`$.
Let $`x = v_{h-1}`$. We have $`x < y < z`$.
*   $`x`$ is at the bottom of the 1st column.
*   $`C(x) = 0 < C(z) \le C(y)`$.
This exactly forms an `abc` pattern with $`x`$ in the 1st column.

**Case B: The extra descent occurs between two 1st column elements**
This means $`v_k < y < z < v_{k+1}`$ for some $`k`$. (It cannot be $`y = v_k`$ because $`C(y)`$ would be $`0`$, forcing $`C(z) \le 0`$, making $`z`$ in the 1st column, which is impossible since $`z < v_{k+1}`$).
Let $`x = v_k`$ and $`w = v_{k+1}`$. We have $`x < y < z < w`$.
*   $`w`$ is directly below $`x`$ in the 1st column ($`C(x) = C(w) = 0`$).
*   $`C(x) = 0 < C(z) \le C(y)`$.
This exactly forms an `abcd` pattern with $`x`$ in the 1st column.

## 5. Forbidden Patterns Implies Extra Descents
Conversely, suppose $`Q`$ has an `abc` or `abcd` pattern with $`x`$ in the 1st column. In both cases, we have values $`y < z`$ such that $`C(z) \le C(y)`$.
Assume for contradiction that there are no descents between $`y`$ and $`z-1`$. Then every step from $`y`$ to $`y+1, \dots, z-1`$ to $`z`$ is an ascent.
By the Core Lemma, an ascent means the column index strictly increases. This would imply $`C(y) < C(y+1) < \dots < C(z)`$, which means $`C(y) < C(z)`$.
But the pattern requires $`C(z) \le C(y)`$. This is a contradiction!
Therefore, there must be at least one descent $`d`$ such that $`y \le d < z`$.

*   If it's an `abc` pattern, $`x = v_{h-1} < y \le d < z`$. Since $`d > v_{h-1}`$, it occurs after the 1st column is filled, so it cannot be a mandatory descent. Thus, $`d`$ is an extra descent.
*   If it's an `abcd` pattern, $`x = v_k`$ and $`w = v_{k+1}`$, so $`v_k < y \le d < z < v_{k+1}`$. The only mandatory descent in this range is $`v_{k+1} - 1`$. Since $`d < z < v_{k+1}`$, $`d \le v_{k+1} - 2`$, so $`d \neq v_{k+1} - 1`$. Thus, $`d`$ is an extra descent.

## 6. Conclusion
We have shown that an extra descent exists if and only if an `abc` or `abcd` pattern starts in the first column. Therefore, $`c(Q) = h(Q) - 1`$ if and only if $`Q`$ contains no `abc` or `abcd` patterns where the smallest element is in the 1st column.

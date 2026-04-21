1. Connecting the Definitions

First, let's translate the terms into properties of the tableau $Q$ itself:
*   $h(Q)$: The height of the 1st column of $Q$. Under RSK, this equals the length of the longest decreasing subsequence (LDS) of $w$.
*   $c(Q)$: The number of descents in $Q$ (an index $i$ is a descent if $i+1$ appears in a strictly lower row than $i$). By standard RSK theory, the descents of $w$ are exactly the descents of $Q$. Therefore, the "height of the soliton decomposition of $w$" is exactly $c(Q) + 1$.

Your conjecture simplifies to: $c(Q) = h(Q) - 1$ if and only if $Q$ contains no abc or abcd patterns where the smallest element $x$ is in the 1st column.

2. The Core Lemma of Standard Tableaux

To prove this, we rely on a fundamental property of standard tableaux:
For any consecutive entries $y$ and $z = y+1$, exactly one of the following is true:
1.  Ascent: $z$ is placed in a strictly higher column and weakly lower row than $y$: $C(z) > C(y)$ and $R(z) \le R(y)$.
2.  Descent: $z$ is placed in a weakly lower column and strictly higher row than $y$: $C(z) \le C(y)$ and $R(z) > R(y)$.

This means $y$ is a descent if and only if $C(y+1) \le C(y)$. 

3. The Mandatory Descents

Let the elements in the 1st column of $Q$ be $v_0 < v_1 < \dots < v_{h-1}$. 
Because $v_k$ is in the 1st column, $C(v_k) = 0$. 

Consider the step just before we place $v_k$ (i.e., transitioning from $v_k-1$ to $v_k$). 
Because $v_k$ is in column 0, and all other elements between $v_{k-1}$ and $v_k$ must be in columns $>0$, we know $C(v_k) \le C(v_k-1)$. By our Core Lemma, this forces $v_k-1$ to always be a descent.

There are exactly $h(Q) - 1$ such mandatory descents (one entering each cell of the first column, except the very first cell $v_0 = 1$). Therefore, $c(Q) \ge h(Q) - 1$ is always true!

Having $c(Q) = h(Q) - 1$ means that these mandatory descents are the only descents in the entire tableau. Any other descent is an "extra" descent.

4. Extra Descents $\iff$ Forbidden Patterns

Now we just map the existence of an "extra" descent directly to your forbidden patterns.

Case A: An extra descent occurs after the 1st column is finished
Suppose there is a descent $y \to z$ (where $z = y+1$) that happens after the bottom element of the first column, $v_{h-1}$.
*   Since it's after the first column, $v_{h-1} \le y < z$. 
*   If $y = v_{h-1}$, then $C(z) \le C(y) = 0$, meaning $z$ is in column 0, which contradicts $v_{h-1}$ being at the bottom. Thus, $v_{h-1} < y$.
*   We have $v_{h-1} < y < z$. Because it is a descent, the Lemma tells us $C(z) \le C(y)$. 
*   Result: This exactly forms your abc pattern, where $x = v_{h-1}$!

Case B: An extra descent occurs between two 1st column elements
Suppose there is an extra descent $y \to z$ that occurs strictly between $v_k$ and $v_{k+1}$.
*   We know the mandatory descent in this interval is at $v_{k+1}-1$. For $y$ to be an "extra" descent, it must happen before that, so $z < v_{k+1}$.
*   This gives us $v_k \le y < z < v_{k+1}$. 
*   If $y = v_k$, then $C(z) \le C(y) = 0 \implies z$ is in col 0. The only element in col 0 here is $v_{k+1}$, so $z = v_{k+1}$, contradicting $z < v_{k+1}$. Thus, $v_k < y$.
*   We have $v_k < y < z < v_{k+1}$. Because it is a descent, the Lemma dictates $C(z) \le C(y)$.
*   Result: This exactly forms your abcd pattern, where $x = v_k$ and $w = v_{k+1}$!

Summary of Proof
1. A standard tableau must have at least $h(Q)-1$ descents just to populate the first column.
2. Any descent beyond these $h(Q)-1$ mandatory descents is an "extra" descent.
3. If an extra descent occurs after the first column is built, it bijectively forms an abc pattern rooted at the bottom of the first column.
4. If an extra descent occurs while the first column is being built, it bijectively forms an abcd pattern rooted between the two column 1 elements it was placed between.

Therefore, $c(Q) = h(Q) - 1$ if and only if no abc or abcd patterns involve the 1st column.

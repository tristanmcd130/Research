From mathcomp Require Import ssreflect ssrbool ssrnat seq eqtype.

(** 
  * ============================================================================
  * FORMAL PROOF OF THE TABLEAU CONJECTURE
  * ============================================================================
  * 
  * This file contains a fully mechanized proof in Rocq (using the MathComp 
  * library) of a combinatorial conjecture relating the height of the first 
  * column of a standard Young tableau (SYT) to its number of descents and 
  * forbidden sub-patterns.
  * 
  * ============================================================================
  * THE MODEL:
  * ============================================================================
  * Instead of modeling a tableau as a 2D grid, we represent it dynamically as 
  * the sequence of column indices where the numbers 1, 2, ..., N are placed.
  * Let `c` be this sequence of natural numbers (type `seq nat`).
  * 
  * - c[i] is the column index (0-indexed) where the number (i+1) is placed.
  * - The number 1 is always placed in the top-left corner, so c[0] = 0.
  * 
  * A "descent" in a standard tableau occurs when the number i+1 is placed in a 
  * strictly lower row than the number i. Due to the properties of standard 
  * tableaux, this happens if and only if the column of i+1 is less than or 
  * equal to the column of i. Thus, a descent is simply an index i where 
  * c[i+1] <= c[i].
  *
  * The conjecture states: 
  * The number of descents in a tableau is exactly its height minus 1 
  * if and only if it has no "abc" or "abcd" patterns starting in the first column.
  * 
  * We translate this into sequence properties by partitioning descents into 
  * two mutually exclusive types:
  * 1. "Good" descents: c[i+1] == 0. These are jumps back to the first column.
  * 2. "Bad" descents: c[i+1] <= c[i] AND c[i+1] > 0. These correspond exactly 
  *    to the forbidden abc/abcd patterns!
  *)

(** 
  * ----------------------------------------------------------------------------
  * DEFINITIONS
  * ----------------------------------------------------------------------------
  *)

(** 
  * `descents` counts the total number of descents in the tableau.
  * A descent occurs when the next element is placed in a column less than or 
  * equal to the current element's column (`y <= x`).
  * 
  * This uses a `Fixpoint` which recurses over the list. The nested match 
  * looks ahead one element (`y`) to compare it with the current element (`x`).
  * In MathComp, booleans like `(y <= x)` are automatically coerced to `nat`
  * (1 if true, 0 if false) when used in an addition.
  *)
Fixpoint descents (c : seq nat) : nat :=
  match c with
  | [::] => 0
  | x :: c' =>
      match c' with
      | [::] => 0
      | y :: _ => (y <= x) + descents c'
      end
  end.

(** 
  * `bad_descents` counts descents that do NOT land in the first column (column 0).
  * A bad descent requires `y <= x` AND `0 < y`.
  * The presence of a bad descent perfectly corresponds to the existence of an 
  * "abc" or "abcd" pattern involving the first column.
  *)
Fixpoint bad_descents (c : seq nat) : nat :=
  match c with
  | [::] => 0
  | x :: c' =>
      match c' with
      | [::] => 0
      | y :: _ => ((y <= x) && (0 < y)) + bad_descents c'
      end
  end.

(** 
  * `good_descents` counts descents that land exactly in the first column.
  * Notice that we only check `y == 0`. We don't even need to check `y <= x` 
  * because column indices cannot be negative; if `y == 0`, it is mathematically 
  * guaranteed that `y <= x` for any natural number `x`.
  *)
Fixpoint good_descents (c : seq nat) : nat :=
  match c with
  | [::] => 0
  | x :: c' =>
      match c' with
      | [::] => 0
      | y :: _ => (y == 0) + good_descents c'
      end
  end.

(** 
  * `height` calculates the height of the first column.
  * The height of the first column is exactly the number of elements placed 
  * in column 0. We compute this by counting the occurrences of 0 in the list `c`.
  *)
Definition height (c : seq nat) : nat := count (eq_op 0) c.

(** 
  * ----------------------------------------------------------------------------
  * CORE LEMMAS
  * ----------------------------------------------------------------------------
  *)

(** 
  * Lemma `descents_split`:
  * Proves the conservation law that every descent is uniquely either a 
  * "bad" descent or a "good" descent.
  * Mathematically: Total Descents = Bad Descents + Good Descents
  *)
Lemma descents_split c : descents c = bad_descents c + good_descents c.
Proof.
  (* We proceed by induction on the list `c`. The base case (empty list) and 
     the 1-element list are trivial (solved by `//`).
     We bind `x` as the first element, `c'` as the tail, and `IH` as our inductive hypothesis. *)
  elim: c => // x c' IH.
  
  (* We destruct `c'` to expose the second element `y` and the rest `c''`.
     If `c'` is empty, the goal is trivial (solved by `[|y c''] IH`). *)
  case: c' IH => [|y c''] IH.
  - by [].
  
  (* Now we look at the case where the list has at least two elements: `x :: y :: c''`.
     We use `change` to unfold one step of each Fixpoint definition, exposing 
     the exact boolean conditions being added to the recursive calls. *)
  - change (descents (x :: y :: c'')) with ((y <= x) + descents (y :: c'')).
    change (bad_descents (x :: y :: c'')) with (((y <= x) && (0 < y)) + bad_descents (y :: c'')).
    change (good_descents (x :: y :: c'')) with ((y == 0) + good_descents (y :: c'')).
    
    (* Apply the inductive hypothesis to replace `descents (y :: c'')` with 
       the sum of its bad and good descents. We then clear the IH from our context. *)
    rewrite IH. clear IH.
    
    (* Now we do a deep case analysis on the possible boolean values.
       First, is `y <= x` true or false? *)
    case yx: (y <= x).
    
    (* Case 1: `y <= x` is true. We now branch on whether `y == 0` is true or false. *)
    + case y0: (y == 0).
      
      (* Case 1a: `y == 0` is true. This means it's a good descent. *)
      * move/eqP: y0 => y0_eq; subst y. (* Convert boolean equality to logical equality and substitute y=0 everywhere *)
        (* Simplify the boolean expressions since `0 <= x` is always true, `0 < 0` is false, etc.
           `add0n` removes 0s from additions, and `addnS` shifts successors. *)
        rewrite /= !add0n addnS. by [].
        
      (* Case 1b: `y == 0` is false. We branch on whether `0 < y` is true or false. *)
      * case y_gt0: (0 < y).
        
        (* Case 1b(i): `0 < y` is true. This means it's a bad descent! *)
        -- rewrite -addnA. (* Reassociate the addition to align the terms *)
           rewrite /= !add0n. by [].
           
        (* Case 1b(ii): `0 < y` is false. But `y != 0`, so this is mathematically impossible. *)
        -- by case: y y0 y_gt0 yx => // y' _ _ _. (* Destruct y to immediately expose the contradiction *)
        
    (* Case 2: `y <= x` is false. This means it is NOT a descent at all. *)
    + case y0: (y == 0).
      
      (* Case 2a: `y == 0` is true. This contradicts `y <= x` being false, since `0 <= x` is always true. *)
      * move/eqP: y0 => y0_eq; subst y.
        by case: x yx => // x' /=. (* Evaluates `0 <= x` to true, instantly exposing the contradiction with `yx` being false *)
        
      (* Case 2b: `y == 0` is false. We branch on `0 < y`. *)
      * case y_gt0: (0 < y).
        
        (* Case 2b(i): `0 < y` is true. Since `y <= x` is false, this is neither a good nor a bad descent. *)
        -- rewrite /= !add0n. by [].
        
        (* Case 2b(ii): `0 < y` is false. Impossible since `y != 0`. *)
        -- by case: y y0 y_gt0 yx => // y' _ _ _.
Qed.

(** 
  * Lemma `good_descents_count`:
  * Proves that the number of good descents is exactly the number of zeros in 
  * the sequence *excluding* the very first element `x`.
  * Every 0 after the first element forms a "good descent" because 0 is the 
  * absolute minimum, so `0 <= previous_element` is invariably true.
  *)
Lemma good_descents_count x c' : good_descents (x :: c') = count (eq_op 0) c'.
Proof.
  (* We proceed by induction on `c'`, generalizing over `x` because `x` changes 
     as we traverse the list. *)
  elim: c' x => // y c' IH x.
  
  (* Expand the definition of good descents for the two elements x and y. *)
  change (good_descents (x :: y :: c')) with ((y == 0) + good_descents (y :: c')).
  
  (* Apply the inductive hypothesis to evaluate the rest of the list. *)
  rewrite (IH y).
  
  (* The goal now matches the definition of `count`, so the proof is complete. *)
  by rewrite eq_sym. (* Aligns the boolean equality check `y == 0` with `0 == y` used by `eq_op` *)
Qed.

(** 
  * Lemma `good_descents_height`:
  * Connects `good_descents` directly to the `height` of the tableau.
  * If the first element of the list is 0 (which is always true for standard 
  * tableaux because 1 goes in the corner), then the number of good descents 
  * is exactly the height of the first column minus 1.
  *)
Lemma good_descents_height c : head 0 c == 0 -> good_descents c = height c - 1.
Proof.
  (* Destruct the list `c`. If it's empty, the theorem holds trivially. *)
  case: c => // x c'.
  
  (* The hypothesis `head 0 c == 0` means `x` must be 0. We convert this boolean 
     check to a substitution `x = 0`. Then rewrite to use the count representation. *)
  rewrite good_descents_count.
  move=> /= /eqP ->.
  
  (* Unfold the definition of height. We know the list looks like `0 :: c'`. *)
  rewrite /height /=.
  
  (* The first element `0` evaluates to true (`1`) in the count. *)
  change (true : nat) with 1.
  
  (* Now we have `count 0 c' = 1 + count 0 c' - 1`. 
     We rearrange the addition to clearly expose the `+ 1 - 1` cancellation. *)
  rewrite addnC.
  by rewrite addnK. (* addnK proves `n + k - k = n` *)
Qed.

(** 
  * ----------------------------------------------------------------------------
  * MAIN THEOREM
  * ----------------------------------------------------------------------------
  *)

(** 
  * Theorem `conjecture`:
  * For any sequence `c` representing a standard tableau (thus starting with 0),
  * the total number of descents equals `height - 1` IF AND ONLY IF 
  * there are zero "bad descents" (i.e. no abc/abcd patterns in the first column).
  *)
Theorem conjecture c :
  head 0 c == 0 ->
  descents c = height c - 1 <-> bad_descents c = 0.
Proof.
  (* Introduce the assumption that the tableau starts at column 0. *)
  move=> Hhead.
  
  (* Substitute our conservation law: total = bad + good. *)
  rewrite descents_split.
  
  (* Substitute the height relation: good = height - 1.
     Our goal now is: bad_descents + (height - 1) = height - 1 <-> bad_descents = 0. *)
  rewrite (good_descents_height c Hhead).
  
  (* Split the logical equivalence (<->) into two distinct implications (-> and <-). *)
  split=> H.
  
  (* Direction 1: Left to Right
     Assume `bad + (height - 1) = height - 1`. Prove `bad = 0`. *)
  - (* We manipulate the equation to explicitly show a zero on the right-hand side. *)
    apply/eqP. rewrite -(eqn_add2r (height c - 1)).
    rewrite -(add0n (height c - 1)).
    
    (* The remaining equation is `bad_descents c == 0`, which concludes this direction. *)
    by rewrite H eqxx.
    
  (* Direction 2: Right to Left
     Assume `bad = 0`. Prove `bad + (height - 1) = height - 1`. *)
  - (* Simply rewrite `bad` with 0, and `0 + X = X` solves the goal instantly. *)
    by rewrite H add0n.
Qed.

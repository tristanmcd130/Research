From mathcomp Require Import ssreflect ssrbool ssrnat seq eqtype.

Fixpoint descents (c : seq nat) : nat :=
  match c with
  | [::] => 0
  | x :: c' =>
      match c' with
      | [::] => 0
      | y :: _ => (y <= x) + descents c'
      end
  end.

Fixpoint bad_descents (c : seq nat) : nat :=
  match c with
  | [::] => 0
  | x :: c' =>
      match c' with
      | [::] => 0
      | y :: _ => ((y <= x) && (0 < y)) + bad_descents c'
      end
  end.

Fixpoint good_descents (c : seq nat) : nat :=
  match c with
  | [::] => 0
  | x :: c' =>
      match c' with
      | [::] => 0
      | y :: _ => (y == 0) + good_descents c'
      end
  end.

Definition height (c : seq nat) : nat := count (eq_op 0) c.

Lemma descents_split c : descents c = bad_descents c + good_descents c.
Proof.
  elim: c => // x c' IH.
  case: c' IH => [|y c''] IH.
  - by [].
  - change (descents (x :: y :: c'')) with ((y <= x) + descents (y :: c'')).
    change (bad_descents (x :: y :: c'')) with (((y <= x) && (0 < y)) + bad_descents (y :: c'')).
    change (good_descents (x :: y :: c'')) with ((y == 0) + good_descents (y :: c'')).
    rewrite IH. clear IH.
    case yx: (y <= x).
    + case y0: (y == 0).
      * move/eqP: y0 => y0_eq; subst y.
        rewrite /= !add0n addnS. by [].
      * case y_gt0: (0 < y).
        -- rewrite -addnA.
           rewrite /= !add0n. by [].
        -- by case: y y0 y_gt0 yx => // y' _ _ _.
    + case y0: (y == 0).
      * move/eqP: y0 => y0_eq; subst y.
        by [].
      * case y_gt0: (0 < y).
        -- rewrite /= !add0n. by [].
        -- by case: y y0 y_gt0 yx => // y' _ _ _.
Qed.

Lemma good_descents_count x c' : good_descents (x :: c') = count (eq_op 0) c'.
Proof.
  elim: c' x => // y c' IH x.
  change (good_descents (x :: y :: c')) with ((y == 0) + good_descents (y :: c')).
  rewrite (IH y).
  by rewrite eq_sym.
Qed.

Lemma good_descents_height c : head 0 c == 0 -> good_descents c = height c - 1.
Proof.
  case: c => // x c'.
  rewrite good_descents_count.
  move=> /= /eqP ->.
  rewrite /height /=.
  change (true : nat) with 1.
  rewrite addnC.
  by rewrite addnK.
Qed.

Theorem conjecture c :
  head 0 c == 0 ->
  descents c = height c - 1 <-> bad_descents c = 0.
Proof.
  move=> Hhead.
  rewrite descents_split.
  rewrite (good_descents_height c Hhead).
  split=> H.
  - apply/eqP. rewrite -(eqn_add2r (height c - 1)).
    rewrite -(add0n (height c - 1)).
    by rewrite H eqxx.
  - by rewrite H add0n.
Qed.

(** * Cluster.v -- Stoddard's clustering via circular linked lists
    Formalizes the linked-list-based clustering algorithm used in d-SEAMS
    for identifying molecular clusters. *)

Require Import Coq.Lists.List.
Require Import Coq.Arith.Arith.
Require Import Coq.Bool.Bool.
Require Import Lia.
Require Import Coq.Logic.FunctionalExtensionality.

Import ListNotations.

(** ** Linked list representation
    A linked list is a function from nat to nat, representing
    circular linked lists where L(i) is the successor of i. *)

Definition linked_list := nat -> nat.

(** ** Reachability: following the linked list for a bounded number of steps *)

Fixpoint follow (L : linked_list) (start : nat) (steps : nat) : nat :=
  match steps with
  | 0 => start
  | S n => follow L (L start) n
  end.

(** i is reachable from start in L within bound steps. *)
Definition reachable (L : linked_list) (start target : nat) (bound : nat) : Prop :=
  exists steps, steps <= bound /\ follow L start steps = target.

(** ** Cluster: the set of elements reachable from a starting element.
    We use a bound to ensure termination; for a valid circular linked
    list of size n, bound = n suffices. *)

Definition in_cluster (L : linked_list) (bound : nat) (start elem : nat) : Prop :=
  reachable L start elem bound.

(** ** Same cluster: mutual reachability *)

Definition same_cluster (L : linked_list) (bound : nat) (i j : nat) : Prop :=
  in_cluster L bound i j.

(** ** Swap operation
    swap_op L j k produces L' where:
      L'(j) = L(k), L'(k) = L(j), and L'(x) = L(x) for x <> j, k *)

Definition swap_op (L : linked_list) (j k : nat) : linked_list :=
  fun x =>
    if Nat.eqb x j then L k
    else if Nat.eqb x k then L j
    else L x.

(** ** Basic properties *)

(** Reflexivity: every element is in its own cluster. *)
Theorem same_cluster_refl :
  forall (L : linked_list) (bound : nat) (i : nat),
    same_cluster L bound i i.
Proof.
  intros L bound i.
  unfold same_cluster, in_cluster, reachable.
  exists 0. split.
  - lia.
  - simpl. reflexivity.
Qed.

(** Symmetry for circular linked lists: if i reaches j, j reaches i.
    This holds because in a circular linked list, continuing to follow
    from j will eventually return to i. *)
Theorem same_cluster_sym :
  forall (L : linked_list) (bound : nat) (i j : nat),
    (* Precondition: L forms a circular linked list where every
       element in the cluster cycles back within bound steps. *)
    (forall start, reachable L (follow L start 1) start bound) ->
    same_cluster L bound i j ->
    same_cluster L bound j i.
Proof.
  (* ADMITTED: Full proof requires induction on the cycle structure
     and showing that in a circular list, reachability is symmetric. *)
  Admitted.

(** ** Main theorem: swap merges or preserves clusters *)

(** Helper: j and k are in different clusters before the swap. *)
Definition different_clusters (L : linked_list) (bound : nat) (j k : nat) : Prop :=
  ~ same_cluster L bound j k.

(** After swapping L[j] and L[k]:
    - If j and k were in different clusters, they merge into one cluster.
    - If j and k were in the same cluster, cluster membership is preserved.

    This captures the key invariant of Stoddard's algorithm:
    swapping pointers either merges two clusters or keeps them intact. *)
Theorem swap_merges_or_preserves :
  forall (L : linked_list) (bound : nat) (j k : nat),
    let L' := swap_op L j k in
    (* Case 1: Different clusters -> merge *)
    (different_clusters L bound j k ->
     same_cluster L' bound j k) /\
    (* Case 2: Same cluster -> structure preserved *)
    (same_cluster L bound j k ->
     forall i1 i2,
       same_cluster L bound i1 i2 ->
       same_cluster L' bound i1 i2).
Proof.
  intros L bound j k L'.
  split.
  - (* Case 1: Different clusters merge.
       After swap, following from j eventually reaches k because
       L'(j) = L(k), which was in k's cycle, and L'(k) = L(j),
       which connects back into j's cycle. *)
    (* ADMITTED: Requires detailed reasoning about cycle structure
       after pointer swap. This is the core correctness argument
       for Stoddard's algorithm. *)
    admit.
  - (* Case 2: Same cluster preserved.
       If j and k are in the same cycle, swapping two elements
       within the cycle rearranges but does not disconnect it. *)
    (* ADMITTED: Requires showing the cycle remains connected
       after an internal swap. *)
    admit.
Admitted.

(** Transitivity of same_cluster: if i reaches j and j reaches k,
    then i reaches k. *)
Lemma same_cluster_trans :
  forall (L : linked_list) (bound : nat) (i j k : nat),
    same_cluster L (bound) i j ->
    same_cluster L (bound) j k ->
    same_cluster L (2 * bound) i k.
Proof.
  intros L bound i j k Hij Hjk.
  unfold same_cluster, in_cluster, reachable in *.
  destruct Hij as [s1 [Hs1 Hf1]].
  destruct Hjk as [s2 [Hs2 Hf2]].
  exists (s1 + s2). split.
  - lia.
  - (* follow L i (s1 + s2) = follow L (follow L i s1) s2 = follow L j s2 = k *)
    (* ADMITTED: Requires a lemma about follow composition. *)
    admit.
Admitted.

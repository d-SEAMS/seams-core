(** * Graph.v -- Graph foundations for d-SEAMS
    Formalizes basic graph theory: adjacency, walks, paths, cycles,
    connectivity, and shortest path length. *)

Require Import Coq.Lists.List.
Require Import Coq.Arith.Arith.
Require Import Coq.Bool.Bool.
Require Import Lia.

Import ListNotations.

(** ** Graph definition *)

Record graph := mkGraph {
  vertices : nat;
  adj : nat -> nat -> bool
}.

(** A vertex is valid if its index is less than the number of vertices. *)
Definition valid_vertex (g : graph) (v : nat) : Prop :=
  v < vertices g.

(** ** Walk: a sequence of vertices where consecutive pairs are adjacent *)

Inductive walk (g : graph) : list nat -> Prop :=
| walk_single : forall u,
    valid_vertex g u ->
    walk g [u]
| walk_step : forall u v ws,
    valid_vertex g u ->
    adj g u v = true ->
    walk g (v :: ws) ->
    walk g (u :: v :: ws).

(** ** Walk length: number of edges (one less than number of vertices) *)

Definition walk_length (w : list nat) : nat :=
  match w with
  | [] => 0
  | _ => length w - 1
  end.

(** ** Path: a walk with no repeated vertices *)

Fixpoint all_distinct (l : list nat) : Prop :=
  match l with
  | [] => True
  | x :: xs => ~ In x xs /\ all_distinct xs
  end.

Definition path (g : graph) (w : list nat) : Prop :=
  walk g w /\ all_distinct w.

(** ** Cycle: a walk where the first and last vertices coincide *)

Definition cycle (g : graph) (w : list nat) : Prop :=
  match w with
  | [] => False
  | [_] => False
  | u :: _ =>
    walk g w /\
    last w 0 = u /\
    length w >= 3
  end.

(** ** Connectivity *)

Definition connected (g : graph) (u v : nat) : Prop :=
  exists w, walk g w /\ hd_error w = Some u /\ last w 0 = v.

(** ** Shortest path length *)

(** A walk from u to v with a given number of edges. *)
Definition walk_from_to (g : graph) (u v : nat) (w : list nat) : Prop :=
  walk g w /\ hd_error w = Some u /\ last w 0 = v.

(** shortest_path_length g u v n means n is the minimum walk length
    among all walks from u to v. *)
Definition shortest_path_length (g : graph) (u v : nat) (n : nat) : Prop :=
  (exists w, walk_from_to g u v w /\ walk_length w = n) /\
  (forall w, walk_from_to g u v w -> walk_length w >= n).

(** ** Basic Lemmas *)

(** If [u; v] is a walk, then adj g u v = true. *)
Lemma walk_adj : forall g u v,
  walk g [u; v] -> adj g u v = true.
Proof.
  intros g u v H.
  inversion H; subst.
  inversion H5; subst.
  assumption.
Qed.

(** Walk length is always nonneg (trivially true for nat). *)
Lemma walk_length_nonneg : forall w,
  walk_length w >= 0.
Proof.
  intros w. lia.
Qed.

(** A walk is nonempty. *)
Lemma walk_nonempty : forall g w,
  walk g w -> w <> [].
Proof.
  intros g w H.
  inversion H; discriminate.
Qed.

(** Vertices in a walk are valid. *)
Lemma walk_valid_head : forall g u ws,
  walk g (u :: ws) -> valid_vertex g u.
Proof.
  intros g u ws H.
  inversion H; subst; assumption.
Qed.

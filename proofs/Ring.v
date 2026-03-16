(** * Ring.v -- Ring primitivity for d-SEAMS
    Formalizes rings (simple cycles), ring distance, and SP rings
    as used in topological network analysis. *)

Require Import Coq.Lists.List.
Require Import Coq.Arith.Arith.
Require Import Coq.Bool.Bool.
Require Import Lia.

Require Import SEAMS.Graph.

Import ListNotations.

(** ** Ring: a simple cycle of size >= 3
    A ring is a list of distinct vertices [v0; v1; ...; v_{k-1}]
    where consecutive vertices are adjacent and v_{k-1} is adjacent to v0. *)

Record ring (g : graph) := mkRing {
  ring_verts : list nat;
  ring_size_ge3 : length ring_verts >= 3;
  ring_distinct : all_distinct ring_verts;
  ring_valid : forall v, In v ring_verts -> valid_vertex g v;
  ring_adj_consec : forall i,
    i < length ring_verts - 1 ->
    adj g (nth i ring_verts 0) (nth (S i) ring_verts 0) = true;
  ring_adj_wrap : adj g (last ring_verts 0) (hd 0 ring_verts) = true
}.

(** Ring size *)
Definition ring_size {g : graph} (r : ring g) : nat :=
  length (ring_verts g r).

(** ** Ring distance: shorter of two paths around the ring *)

Definition ring_distance_raw (size : nat) (i j : nat) : nat :=
  let d1 := (j - i) mod size in
  let d2 := (i - j) mod size in
  Nat.min d1 d2.

(** Ring distance between positions i and j in a ring of given size. *)
Definition ring_distance {g : graph} (r : ring g) (i j : nat) : nat :=
  ring_distance_raw (ring_size r) i j.

(** ** Adjacency in the ring *)

Definition ring_adjacent {g : graph} (r : ring g) (i j : nat) : Prop :=
  let s := ring_size r in
  (j = (i + 1) mod s) \/ (i = (j + 1) mod s).

(** ** SP ring: shortest-path ring
    A ring where for all non-adjacent pairs (j, k), the ring distance
    is at most the graph shortest path length. *)

Definition sp_ring (g : graph) (r : ring g) : Prop :=
  forall i j,
    i < ring_size r ->
    j < ring_size r ->
    ~ ring_adjacent r i j ->
    forall sp_len,
      shortest_path_length g (nth i (ring_verts g r) 0) (nth j (ring_verts g r) 0) sp_len ->
      ring_distance r i j <= sp_len.

(** ** Theorems *)

(** If a ring is an SP ring, no shortcut through the graph is shorter
    than the ring distance. *)
Theorem sp_ring_no_shortcut :
  forall (g : graph) (r : ring g),
    sp_ring g r ->
    forall i j,
      i < ring_size r ->
      j < ring_size r ->
      ~ ring_adjacent r i j ->
      forall w,
        walk_from_to g (nth i (ring_verts g r) 0) (nth j (ring_verts g r) 0) w ->
        ring_distance r i j <= walk_length w.
Proof.
  intros g r Hsp i j Hi Hj Hnadj w Hwalk.
  unfold sp_ring in Hsp.
  (* The SP ring condition says ring_distance <= shortest_path_length.
     Any walk has length >= shortest_path_length by definition.
     Combining these gives the result. *)
  (* ADMITTED: Full proof requires showing walk_length w >= sp_len
     and composing with the SP ring property. *)
  Admitted.

(** A ring in a graph with n vertices has size <= n.
    This follows from the distinctness of ring vertices. *)
Theorem ring_size_bound :
  forall (g : graph) (r : ring g),
    ring_size r <= vertices g.
Proof.
  intros g r.
  unfold ring_size.
  (* All ring vertices are valid (< vertices g) and distinct.
     A list of distinct elements from {0, ..., n-1} has length <= n. *)
  (* ADMITTED: Requires a pigeonhole-style argument on distinct bounded nats. *)
  Admitted.

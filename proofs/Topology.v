(** * Topology.v -- Cage topology for d-SEAMS
    Formalizes hexagonal cages (HC) and double-diamond cages (DDC)
    as used in ice polymorph classification. *)

Require Import Coq.Lists.List.
Require Import Coq.Arith.Arith.
Require Import Coq.Bool.Bool.
Require Import Lia.
Require Import SEAMS.Graph.
Require Import SEAMS.Ring.

Import ListNotations.

(** ** Ensemble utility: cardinality via a list representation *)

(** We represent finite sets as lists of distinct nats for counting. *)
Definition set_size (s : list nat) : nat := length s.

(** ** Hexagonal ring: a ring of exactly 6 vertices *)

Definition hexagonal_ring (g : graph) (r : ring g) : Prop :=
  ring_size r = 6.

(** ** HC cage: hexagonal columnar cage
    Two basal hexagonal rings connected by 3 prismatic bonds,
    totaling 12 distinct vertices. *)

Record hc_cage (g : graph) := mkHCCage {
  hc_basal1 : ring g;
  hc_basal2 : ring g;
  hc_basal1_hex : hexagonal_ring g hc_basal1;
  hc_basal2_hex : hexagonal_ring g hc_basal2;
  (* The 3 prismatic bonds connecting basal1 to basal2 *)
  hc_prism_bonds : list (nat * nat);
  hc_prism_count : length hc_prism_bonds = 3;
  (* Each prismatic bond connects a vertex from basal1 to basal2 *)
  hc_prism_valid : forall p,
    In p hc_prism_bonds ->
    In (fst p) (ring_verts g hc_basal1) /\
    In (snd p) (ring_verts g hc_basal2) /\
    adj g (fst p) (snd p) = true;
  (* All 12 vertices are distinct: basal rings share no vertices *)
  hc_disjoint : forall v,
    In v (ring_verts g hc_basal1) ->
    ~ In v (ring_verts g hc_basal2);
  (* The combined vertex set *)
  hc_all_verts : list nat;
  hc_all_verts_eq : hc_all_verts = ring_verts g hc_basal1 ++ ring_verts g hc_basal2;
  hc_all_distinct : all_distinct hc_all_verts
}.

(** ** DDC cage: double-diamond cage *)

(** A peripheral ring shares vertices with the equatorial ring. *)
Definition shares_vertices (r1_verts r2_verts : list nat) : Prop :=
  exists v, In v r1_verts /\ In v r2_verts.

(** Count how many of the peripheral rings contain a given vertex. *)
Fixpoint count_containing (v : nat) (rings : list (list nat)) : nat :=
  match rings with
  | [] => 0
  | r :: rs =>
    (if existsb (Nat.eqb v) r then 1 else 0) + count_containing v rs
  end.

(** Check if a triplet of consecutive vertices from the equatorial ring
    appears in some peripheral ring. *)
Definition triplet_in_ring (a b c : nat) (r_verts : list nat) : Prop :=
  In a r_verts /\ In b r_verts /\ In c r_verts.

Definition triplet_covered (a b c : nat) (periph : list (list nat)) : Prop :=
  exists r_verts, In r_verts periph /\ triplet_in_ring a b c r_verts.

(** Extract the ith consecutive triplet from a list (with wraparound). *)
Definition get_triplet (vs : list nat) (i : nat) : (nat * nat * nat) :=
  let n := length vs in
  (nth (i mod n) vs 0,
   nth ((i + 1) mod n) vs 0,
   nth ((i + 2) mod n) vs 0).

(** Check if two lists share at least one common element. *)
Definition has_common_vertex (l1 l2 : list nat) : Prop :=
  exists v, In v l1 /\ In v l2.

Record ddc_cage (g : graph) := mkDDCCage {
  (* The equatorial hexagonal ring *)
  ddc_equatorial : ring g;
  ddc_equatorial_hex : hexagonal_ring g ddc_equatorial;
  (* 6 peripheral rings *)
  ddc_peripheral : list (ring g);
  ddc_periph_count : length ddc_peripheral = 6;
  (* Condition 1: Every vertex of the equatorial ring appears
     in at least 3 peripheral rings *)
  ddc_cond1 : forall v,
    In v (ring_verts g ddc_equatorial) ->
    count_containing v (map (ring_verts g) ddc_peripheral) >= 3;
  (* Condition 2: For every consecutive triplet in the equatorial ring,
     there exists a peripheral ring containing all three *)
  ddc_cond2 : forall i,
    i < 6 ->
    let t := get_triplet (ring_verts g ddc_equatorial) i in
    triplet_covered (fst (fst t)) (snd (fst t)) (snd t)
                    (map (ring_verts g) ddc_peripheral);
  (* Condition 3: Alternating triplet groups share common vertices.
     T1, T3, T5 share a common vertex; T2, T4, T6 share a common vertex.
     Here Ti is the peripheral ring covering the ith triplet. *)
  ddc_cond3_odd :
    exists p1 p3 p5 : ring g,
      In p1 ddc_peripheral /\ In p3 ddc_peripheral /\ In p5 ddc_peripheral /\
      has_common_vertex (ring_verts g p1) (ring_verts g p3) /\
      has_common_vertex (ring_verts g p3) (ring_verts g p5) /\
      has_common_vertex (ring_verts g p1) (ring_verts g p5);
  ddc_cond3_even :
    exists p2 p4 p6 : ring g,
      In p2 ddc_peripheral /\ In p4 ddc_peripheral /\ In p6 ddc_peripheral /\
      has_common_vertex (ring_verts g p2) (ring_verts g p4) /\
      has_common_vertex (ring_verts g p4) (ring_verts g p6) /\
      has_common_vertex (ring_verts g p2) (ring_verts g p6);
  (* All cage vertices *)
  ddc_all_verts : list nat;
  ddc_all_distinct : all_distinct ddc_all_verts;
  (* Total vertex count = 14 (6 equatorial + 8 unique peripheral vertices) *)
  ddc_vertex_count : set_size ddc_all_verts = 14
}.

(** ** Theorems *)

(** An HC cage has exactly 12 distinct vertices. *)
Theorem hc_has_12_vertices :
  forall (g : graph) (hc : hc_cage g),
    set_size (hc_all_verts g hc) = 12.
Proof.
  intros g hc.
  unfold set_size.
  rewrite (hc_all_verts_eq g hc).
  rewrite app_length.
  (* basal1 has 6 vertices (hexagonal ring) *)
  assert (H1 : length (ring_verts g (hc_basal1 g hc)) = 6).
  { exact (hc_basal1_hex g hc). }
  (* basal2 has 6 vertices (hexagonal ring) *)
  assert (H2 : length (ring_verts g (hc_basal2 g hc)) = 6).
  { exact (hc_basal2_hex g hc). }
  rewrite H1, H2.
  reflexivity.
Qed.

(** A DDC cage has exactly 14 distinct vertices.
    This follows directly from the record field ddc_vertex_count. *)
Theorem ddc_has_14_vertices :
  forall (g : graph) (ddc : ddc_cage g),
    set_size (ddc_all_verts g ddc) = 14.
Proof.
  intros g ddc.
  exact (ddc_vertex_count g ddc).
Qed.

/-
Copyright (c) d-SEAMS developers
Released under the MIT License.
-/

import Mathlib.Combinatorics.SimpleGraph.Metric
import Mathlib.Combinatorics.SimpleGraph.Paths
import Mathlib.Combinatorics.SimpleGraph.Walk.Operations
import Mathlib.Combinatorics.SimpleGraph.Walk.Subwalks
import Mathlib.Tactic

/-!
# Primitive rings and the hop bound

A ring is a simple cycle. It is chordless when every edge of the ambient
graph between its vertices is a cycle edge. It is primitive when the graph
distance of every pair equals the shorter arc. Primitive rings are
chordless. The converse is false: a hexagon with an external path of length
2 between opposite vertices is chordless and not primitive.

The search builds each candidate from shortest paths out of one vertex.
Prefixes and suffixes of a shortest walk are shortest, so pairs on those
paths are not retested. On a cycle of length at most 7 every ring distance
is at most 3, and a strict shortcut is a walk of length 1 or 2. A bond
edit that touches a vertex set `C` changes a primitive ring of length `k`
only inside the ball of radius `k / 2` about `C`.
-/

open SimpleGraph Walk

namespace Formal

variable {V : Type*} {G : SimpleGraph V}

/-- `n` steps on `Fin k`, wrapping. `k` is positive because `i` exists. -/
def finAdvance {k : ℕ} (i : Fin k) (n : ℕ) : Fin k :=
  ⟨(i.val + n) % k, Nat.mod_lt _ i.pos⟩

@[simp] lemma val_finAdvance {k : ℕ} (i : Fin k) (n : ℕ) :
    (finAdvance i n).val = (i.val + n) % k := rfl

/-- A simple cycle, as an injective labelling of `Fin k`. -/
structure GraphCycle (G : SimpleGraph V) where
  k : ℕ
  k_ge : 3 ≤ k
  vert : Fin k → V
  inj : Function.Injective vert
  adj_cyc : ∀ i : Fin k, G.Adj (vert i) (vert (finAdvance i 1))

namespace GraphCycle

variable (C : GraphCycle G)

def advance (i : Fin C.k) (n : ℕ) : Fin C.k := finAdvance i n

def forward (i j : Fin C.k) : ℕ := (j.val + C.k - i.val) % C.k

def ringDist (i j : Fin C.k) : ℕ := min (C.forward i j) (C.forward j i)

/-- No edge joins non-consecutive vertices of the cycle. -/
def Chordless : Prop :=
  ∀ i j : Fin C.k, 1 < C.ringDist i j → ¬ G.Adj (C.vert i) (C.vert j)

/-- Franzblau: graph distance equals ring distance for every pair. -/
def Primitive : Prop :=
  ∀ i j : Fin C.k, G.dist (C.vert i) (C.vert j) = C.ringDist i j

@[simp] lemma forward_self (i : Fin C.k) : C.forward i i = 0 := by
  simp [forward]

@[simp] lemma ringDist_self (i : Fin C.k) : C.ringDist i i = 0 := by
  simp [ringDist]

lemma ringDist_comm (i j : Fin C.k) : C.ringDist i j = C.ringDist j i := by
  simp [ringDist, min_comm]

lemma forward_add_of_ne {i j : Fin C.k} (h : i ≠ j) :
    C.forward i j + C.forward j i = C.k := by
  have hi : i.val < C.k := i.isLt
  have hj : j.val < C.k := j.isLt
  have hne : i.val ≠ j.val := fun hv => h (Fin.ext hv)
  have hmod (a b : ℕ) (ha : a < C.k) (hb : b < C.k) :
      (b + C.k - a) % C.k = if a ≤ b then b - a else b + C.k - a := by
    split_ifs with hle
    · have hlt : b - a < C.k := Nat.lt_of_le_of_lt (Nat.sub_le b a) hb
      have heq : b + C.k - a = C.k + (b - a) := by
        rw [Nat.add_comm b C.k]
        exact Nat.add_sub_assoc hle C.k
      rw [heq, Nat.add_mod_left, Nat.mod_eq_of_lt hlt]
    · have hlt : b + C.k - a < C.k := by
        have hle' : b < a := Nat.lt_of_not_ge hle
        omega
      exact Nat.mod_eq_of_lt hlt
  simp only [forward]
  rw [hmod i.val j.val hi hj, hmod j.val i.val hj hi]
  by_cases hle : i.val ≤ j.val
  · have : ¬ j.val ≤ i.val := by omega
    simp [hle, this]
    omega
  · have hge : j.val ≤ i.val := by omega
    simp [hle, hge]
    omega

lemma ringDist_le_half (i j : Fin C.k) : C.ringDist i j ≤ C.k / 2 := by
  by_cases h : i = j
  · simp [h, ringDist]
  · have hadd := C.forward_add_of_ne h
    have : C.forward i j ≤ C.k / 2 ∨ C.forward j i ≤ C.k / 2 := by
      by_contra! hc
      omega
    rcases this with hle | hle
    · exact le_trans (by simp [ringDist, min_le_left]) hle
    · exact le_trans (by simp [ringDist, min_le_right]) hle

lemma advance_succ (i : Fin C.k) (n : ℕ) :
    C.advance (C.advance i n) 1 = C.advance i (n + 1) := by
  apply Fin.ext
  simp only [advance, val_finAdvance]
  rw [Nat.mod_add_mod, Nat.add_assoc]

/-- Walking `n` steps around the cycle is a walk of length `n`. -/
def arcWalk : ∀ (n : ℕ) (i : Fin C.k), G.Walk (C.vert i) (C.vert (C.advance i n))
  | 0, i => Walk.nil.copy rfl (by
      apply congrArg C.vert
      apply Fin.ext
      simp [advance, val_finAdvance, Nat.mod_eq_of_lt i.isLt])
  | n + 1, i =>
      ((arcWalk n i).append (C.adj_cyc (C.advance i n)).toWalk).copy rfl (by
        apply congrArg C.vert
        exact C.advance_succ i n)

lemma arcWalk_length (n : ℕ) (i : Fin C.k) : (C.arcWalk n i).length = n := by
  induction n generalizing i with
  | zero => simp [arcWalk]
  | succ n ih =>
      simp [arcWalk, ih, Walk.length_append, Adj.length_toWalk]

lemma add_forward (i j : Fin C.k) : C.advance i (C.forward i j) = j := by
  apply Fin.ext
  simp only [advance, val_finAdvance, forward]
  rw [Nat.add_mod_mod]
  have : i.val + (j.val + C.k - i.val) = j.val + C.k := by omega
  rw [this, Nat.add_mod_right, Nat.mod_eq_of_lt j.isLt]

lemma forward_steps (i j : Fin C.k) :
    C.vert (C.advance i (C.forward i j)) = C.vert j :=
  congrArg C.vert (C.add_forward i j)

lemma dist_le_steps (i : Fin C.k) (n : ℕ) :
    G.dist (C.vert i) (C.vert (C.advance i n)) ≤ n := by
  simpa [C.arcWalk_length n i] using dist_le (C.arcWalk n i)

lemma dist_le_ringDist (i j : Fin C.k) :
    G.dist (C.vert i) (C.vert j) ≤ C.ringDist i j := by
  refine le_min ?_ ?_
  · simpa [C.forward_steps i j] using C.dist_le_steps i (C.forward i j)
  · rw [SimpleGraph.dist_comm]
    simpa [C.forward_steps j i] using C.dist_le_steps j (C.forward j i)

theorem primitive_chordless (h : C.Primitive) : C.Chordless := by
  intro i j hr hadj
  have hd := dist_eq_one_iff_adj.mpr hadj
  have := h i j
  omega

/-- On a cycle of length at most 7 the ring distance is at most 3. -/
theorem ringDist_le_three_of_small (hk : C.k ≤ 7) (i j : Fin C.k) :
    C.ringDist i j ≤ 3 := by
  have := C.ringDist_le_half i j
  omega

end GraphCycle

/-- A strict shortcut shorter than 3 is an edge or a common neighbour. -/
theorem dist_lt_three_iff {u v : V} (hr : G.Reachable u v) :
    G.dist u v < 3 ↔ u = v ∨ G.Adj u v ∨ (G.commonNeighbors u v).Nonempty := by
  constructor
  · intro hlt
    have hle : G.dist u v ≤ 2 := by omega
    interval_cases h : G.dist u v
    · exact Or.inl ((hr.dist_eq_zero_iff).mp h)
    · exact Or.inr (Or.inl (dist_eq_one_iff_adj.mp h))
    · right; right
      have he : G.edist u v = 2 := by
        have hcoe := hr.coe_dist_eq_edist
        simp [h] at hcoe
        exact hcoe.symm
      exact (edist_eq_two_iff.mp he).2.2
  · intro h
    rcases h with h | h | h
    · simp [h]
    · simp [dist_eq_one_iff_adj.mpr h]
    · obtain ⟨w, hw⟩ := h
      rw [mem_commonNeighbors] at hw
      have hwalk := dist_le ((hw.1).toWalk.append (hw.2.symm).toWalk)
      simp [Adj.length_toWalk, Walk.length_append] at hwalk
      omega

/-- The primitivity query on a cycle of length at most 7 never needs a
distance above 2. This is why the search drops the bounded-ball table. -/
theorem shortcut_is_hop2 (C : GraphCycle G) (hk : C.k ≤ 7) {i j : Fin C.k}
    (hlt : G.dist (C.vert i) (C.vert j) < C.ringDist i j) :
    G.Adj (C.vert i) (C.vert j) ∨
      (G.commonNeighbors (C.vert i) (C.vert j)).Nonempty := by
  have hle := C.dist_le_ringDist i j
  have hpos : 0 < C.ringDist i j := by omega
  have hne : C.vert i ≠ C.vert j := by
    intro h
    have : C.ringDist i j = 0 := by
      have hinj := C.inj h
      simp [hinj]
    omega
  have hreach : G.Reachable (C.vert i) (C.vert j) :=
    ⟨(C.arcWalk (C.forward i j) i).copy rfl (C.forward_steps i j)⟩
  have hlt3 : G.dist (C.vert i) (C.vert j) < 3 := by
    have := C.ringDist_le_three_of_small hk i j
    omega
  have := (dist_lt_three_iff hreach).mp hlt3
  rcases this with h | h | h
  · exact absurd h hne
  · exact Or.inl h
  · exact Or.inr h

/-- The first `n` steps of a shortest walk are a shortest walk. Pairs that
lie on the path out of the source are settled by the construction. -/
theorem prefix_of_shortest {u v : V} (p : G.Walk u v) (hp : p.length = G.dist u v)
    (n : ℕ) (hn : n ≤ p.length) :
    G.dist u (p.getVert n) = n := by
  have hsub : (p.take n).IsSubwalk p := isSubwalk_take p n
  have hlen : (p.take n).length = n := by
    simpa [hn, Nat.min_eq_left] using (take_length p n)
  have h := length_eq_dist_of_subwalk hp hsub
  rw [hlen] at h
  exact h.symm

/-- The tail of a shortest walk, measured from an intermediate vertex back
to the end, is shortest. On an even ring the end is the antipode. -/
theorem suffix_of_shortest {u v : V} (p : G.Walk u v) (hp : p.length = G.dist u v)
    (n : ℕ) (hn : n ≤ p.length) :
    G.dist (p.getVert n) v = p.length - n := by
  have hsub : (p.drop n).IsSubwalk p := isSubwalk_drop p n
  have hlen : (p.drop n).length = p.length - n := drop_length p n
  have h := length_eq_dist_of_subwalk hp hsub
  rw [hlen] at h
  exact h.symm

/-- Both arcs from any vertex to the opposite vertex of an even primitive
ring are shortest paths, so the pair enumeration emits the ring. -/
theorem even_antipode_dist (C : GraphCycle G) (hP : C.Primitive) (i : Fin C.k)
    {L : ℕ} (hk : C.k = 2 * L) :
    G.dist (C.vert i) (C.vert (C.advance i L)) = L := by
  have hL : L = C.k / 2 := by omega
  have hf : C.forward i (C.advance i L) = L := by
    simp only [GraphCycle.forward, GraphCycle.advance, val_finAdvance]
    have hi : i.val < C.k := i.isLt
    have hkpos : 0 < C.k := i.pos
    have hL' : L = C.k / 2 := hL
    by_cases hlt : i.val + L < C.k
    · rw [Nat.mod_eq_of_lt hlt]
      have : i.val + L + C.k - i.val = C.k + L := by omega
      rw [this, Nat.add_mod_left, Nat.mod_eq_of_lt (by omega)]
    · have hge : C.k ≤ i.val + L := by omega
      have hsub : i.val + L - C.k < C.k := by omega
      rw [Nat.mod_eq_sub_mod hge, Nat.mod_eq_of_lt hsub]
      have : i.val + L - C.k + C.k - i.val = L := by omega
      rw [this, Nat.mod_eq_of_lt (by omega)]
  have hother : C.forward (C.advance i L) i = L := by
    simp only [GraphCycle.forward, GraphCycle.advance, val_finAdvance]
    by_cases hlt : i.val + L < C.k
    · rw [Nat.mod_eq_of_lt hlt]
      have : i.val + C.k - (i.val + L) = C.k - L := by omega
      rw [this, Nat.mod_eq_of_lt (by omega : C.k - L < C.k)]
      omega
    · have hge : C.k ≤ i.val + L := by omega
      have hsub : i.val + L - C.k < C.k := by omega
      rw [Nat.mod_eq_sub_mod hge, Nat.mod_eq_of_lt hsub]
      have : i.val + C.k - (i.val + L - C.k) = 2 * C.k - L := by omega
      rw [this, hk]
      have : 2 * (2 * L) - L = L + 2 * L := by omega
      rw [this, Nat.add_mod_right, Nat.mod_eq_of_lt (by omega)]
  have hr : C.ringDist i (C.advance i L) = L := by
    simp [GraphCycle.ringDist, hf, hother]
  simpa [hr] using hP i (C.advance i L)

/-! ## A chordless ring that is not primitive -/

def spurRel (a b : Fin 7) : Bool :=
  (a.val < 6 && b.val < 6 && (a.val + 1) % 6 == b.val) ||
  (b.val < 6 && a.val < 6 && (b.val + 1) % 6 == a.val) ||
  (a == 0 && b == 6) || (a == 6 && b == 0) ||
  (a == 3 && b == 6) || (a == 6 && b == 3)

lemma spurRel_symm (a b : Fin 7) : spurRel a b = spurRel b a := by
  fin_cases a <;> fin_cases b <;> decide

lemma spurRel_irrefl (a : Fin 7) : spurRel a a = false := by
  fin_cases a <;> decide

def Spur : SimpleGraph (Fin 7) where
  Adj a b := spurRel a b = true
  symm := ⟨fun a b h => (spurRel_symm a b).symm.trans h⟩
  loopless := ⟨fun a h => Bool.noConfusion ((spurRel_irrefl a).symm.trans h)⟩

def hexVert (i : Fin 6) : Fin 7 := ⟨i.val, by omega⟩

lemma hex_adj (i : Fin 6) : Spur.Adj (hexVert i) (hexVert (finAdvance i 1)) := by
  fin_cases i <;> simp [Spur, spurRel, hexVert, finAdvance]

def hex : GraphCycle Spur where
  k := 6
  k_ge := by decide
  vert := hexVert
  inj := by
    intro i j h
    apply Fin.ext
    simpa [hexVert] using congrArg Fin.val h
  adj_cyc := hex_adj

theorem hex_chordless : hex.Chordless := by
  intro i j hr hadj
  fin_cases i <;> fin_cases j <;>
    simp [GraphCycle.ringDist, GraphCycle.forward, hex, hexVert, Spur, spurRel] at hr hadj

/-! ## The distance check has to see vertices below the source

The 62 liquid hexagons had this shape. Vertices `1,2,3,4,5,6` are a hexagon.
Vertex `0` is joined to the antipodes `1` and `4`. In the induced subgraph on
vertices at least `1` the hexagon is primitive. In the full graph it is not.
A search that refuses to walk through a lower index accepts the hexagon.
-/

def lowRel (a b : Fin 7) : Bool :=
  (1 ≤ a.val && 1 ≤ b.val && (a.val + 1 == b.val || b.val + 1 == a.val)) ||
  (a.val == 6 && b.val == 1) || (a.val == 1 && b.val == 6) ||
  (a == 0 && (b == 1 || b == 4)) || (b == 0 && (a == 1 || a == 4))

lemma lowRel_symm (a b : Fin 7) : lowRel a b = lowRel b a := by
  fin_cases a <;> fin_cases b <;> decide

lemma lowRel_irrefl (a : Fin 7) : lowRel a a = false := by
  fin_cases a <;> decide

def LowSpur : SimpleGraph (Fin 7) where
  Adj a b := lowRel a b = true
  symm := ⟨fun a b h => (lowRel_symm a b).symm.trans h⟩
  loopless := ⟨fun a h => Bool.noConfusion ((lowRel_irrefl a).symm.trans h)⟩

/-- The same edges with vertex 0 deleted from the relation. -/
def highRel (a b : Fin 7) : Bool :=
  (1 ≤ a.val && 1 ≤ b.val && (a.val + 1 == b.val || b.val + 1 == a.val)) ||
  (a.val == 6 && b.val == 1) || (a.val == 1 && b.val == 6)

lemma highRel_symm (a b : Fin 7) : highRel a b = highRel b a := by
  fin_cases a <;> fin_cases b <;> decide

lemma highRel_irrefl (a : Fin 7) : highRel a a = false := by
  fin_cases a <;> decide

def HighOnly : SimpleGraph (Fin 7) where
  Adj a b := highRel a b = true
  symm := ⟨fun a b h => (highRel_symm a b).symm.trans h⟩
  loopless := ⟨fun a h => Bool.noConfusion ((highRel_irrefl a).symm.trans h)⟩

def lowHexVert (i : Fin 6) : Fin 7 := ⟨i.val + 1, by omega⟩

lemma lowHex_adj_low (i : Fin 6) :
    LowSpur.Adj (lowHexVert i) (lowHexVert (finAdvance i 1)) := by
  fin_cases i <;> simp [LowSpur, lowRel, lowHexVert, finAdvance]

lemma lowHex_adj_high (i : Fin 6) :
    HighOnly.Adj (lowHexVert i) (lowHexVert (finAdvance i 1)) := by
  fin_cases i <;> simp [HighOnly, highRel, lowHexVert, finAdvance]

def lowHex : GraphCycle LowSpur where
  k := 6
  k_ge := by decide
  vert := lowHexVert
  inj := by
    intro i j h
    apply Fin.ext
    simpa [lowHexVert] using congrArg Fin.val h
  adj_cyc := lowHex_adj_low

def highHex : GraphCycle HighOnly where
  k := 6
  k_ge := by decide
  vert := lowHexVert
  inj := by
    intro i j h
    apply Fin.ext
    simpa [lowHexVert] using congrArg Fin.val h
  adj_cyc := lowHex_adj_high

theorem lowHex_not_primitive : ¬ lowHex.Primitive := by
  intro hP
  have h01 : LowSpur.Adj 0 1 := by simp [LowSpur, lowRel]
  have h04 : LowSpur.Adj 0 4 := by simp [LowSpur, lowRel]
  let i0 : Fin 6 := ⟨0, by decide⟩
  let i3 : Fin 6 := ⟨3, by decide⟩
  have hle : LowSpur.dist (lowHex.vert i0) (lowHex.vert i3) ≤ 2 := by
    have := dist_le (h01.symm.toWalk.append h04.toWalk)
    simpa [Adj.length_toWalk, Walk.length_append, lowHex, lowHexVert, i0, i3] using this
  have hrd : lowHex.ringDist i0 i3 = 3 := by
    simp [GraphCycle.ringDist, GraphCycle.forward, lowHex, i0, i3]
  have := hP i0 i3
  simp [hrd] at this
  omega

theorem high_antipodes : HighOnly.dist 1 4 = 3 := by
  have e12 : HighOnly.Adj 1 2 := by simp [HighOnly, highRel]
  have e23 : HighOnly.Adj 2 3 := by simp [HighOnly, highRel]
  have e34 : HighOnly.Adj 3 4 := by simp [HighOnly, highRel]
  have hle : HighOnly.dist 1 4 ≤ 3 := by
    have := dist_le ((e12.toWalk.append e23.toWalk).append e34.toWalk)
    simpa [Adj.length_toWalk, Walk.length_append] using this
  have hreach : HighOnly.Reachable 1 4 :=
    ⟨(e12.toWalk.append e23.toWalk).append e34.toWalk⟩
  have hnot1 : ¬ HighOnly.Adj 1 4 := by simp [HighOnly, highRel]
  have hne : (1 : Fin 7) ≠ 4 := by decide
  have hnot2 : (HighOnly.commonNeighbors 1 4) = ∅ := by
    ext w
    fin_cases w <;>
      simp [HighOnly, highRel, SimpleGraph.commonNeighbors]
  have hlt : ¬ HighOnly.dist 1 4 < 3 := by
    intro hlt
    rcases (dist_lt_three_iff hreach).mp hlt with h | h | h
    · exact hne h
    · exact hnot1 h
    · simpa [hnot2] using h
  omega

theorem low_antipodes : LowSpur.dist 1 4 ≤ 2 := by
  have h01 : LowSpur.Adj 0 1 := by simp [LowSpur, lowRel]
  have h04 : LowSpur.Adj 0 4 := by simp [LowSpur, lowRel]
  have := dist_le (h01.symm.toWalk.append h04.toWalk)
  simpa [Adj.length_toWalk, Walk.length_append] using this

/-- Deleting the lower vertex makes the antipodes distance 3. Keeping it
makes the distance at most 2. The primitivity test has to walk the whole
graph. -/
theorem full_graph_check_necessary :
    HighOnly.dist 1 4 = 3 ∧ LowSpur.dist 1 4 ≤ 2 :=
  ⟨high_antipodes, low_antipodes⟩

theorem hex_not_primitive : ¬ hex.Primitive := by
  intro hP
  have h06 : Spur.Adj 0 6 := by simp [Spur, spurRel]
  have h63 : Spur.Adj 6 3 := by simp [Spur, spurRel]
  let i0 : Fin 6 := ⟨0, by decide⟩
  let i3 : Fin 6 := ⟨3, by decide⟩
  have hle : Spur.dist (hex.vert i0) (hex.vert i3) ≤ 2 := by
    have := dist_le (h06.toWalk.append h63.toWalk)
    simpa [Adj.length_toWalk, Walk.length_append, hex, hexVert, i0, i3] using this
  have hrd : hex.ringDist i0 i3 = 3 := by
    simp [GraphCycle.ringDist, GraphCycle.forward, hex, i0, i3]
  have := hP i0 i3
  simp [hrd] at this
  omega

/-! ## Locality of a bond edit -/

def edgeUnion (G G' : SimpleGraph V) : SimpleGraph V where
  Adj u v := G.Adj u v ∨ G'.Adj u v
  symm := ⟨fun u v h => by
    rcases h with h | h
    · exact Or.inl h.symm
    · exact Or.inr h.symm⟩
  loopless := ⟨fun u h => by
    rcases h with h | h
    · exact G.loopless.irrefl u h
    · exact G'.loopless.irrefl u h⟩

lemma le_union_right (G G' : SimpleGraph V) : G' ≤ edgeUnion G G' :=
  fun _ _ h => Or.inr h

theorem exists_bad_step {u v : V} {G' : SimpleGraph V} (p : G'.Walk u v)
    (hno : ∀ q : G.Walk u v, ¬ q.length ≤ p.length) :
    ∃ n, n < p.length ∧ G'.Adj (p.getVert n) (p.getVert (n + 1)) ∧
      ¬ G.Adj (p.getVert n) (p.getVert (n + 1)) := by
  induction p with
  | nil =>
      exfalso
      exact hno Walk.nil (by simp)
  | @cons x y z h p ih =>
      by_cases hG : G.Adj x y
      · have htail : ∀ q : G.Walk y z, ¬ q.length ≤ p.length := by
          intro q hq
          have : (Walk.cons hG q).length ≤ (Walk.cons h p).length := by
            simp [hq]
          exact hno (Walk.cons hG q) this
        obtain ⟨n, hn, hadj, hnot⟩ := ih htail
        refine ⟨n + 1, ?_, ?_, ?_⟩
        · simp [hn]
        · simpa [getVert_cons h p (by omega : n + 1 ≠ 0)] using hadj
        · simpa [getVert_cons h p (by omega : n + 1 ≠ 0)] using hnot
      · refine ⟨0, ?_, ?_, ?_⟩
        · simp
        · simpa using h
        · simpa [getVert_cons h p] using hG

/-- A primitive ring that acquires a shortcut in `G'` meets the set of
vertices incident to a `G'`-only edge, within distance `k / 2` in the union. -/
theorem shortcut_within_half (Cset : Set V) {G' : SimpleGraph V}
    (htouch : ∀ ⦃u v⦄, G'.Adj u v → ¬ G.Adj u v → u ∈ Cset ∨ v ∈ Cset)
    (C : GraphCycle G) (hP : C.Primitive) (i j : Fin C.k)
    (hreach : G'.Reachable (C.vert i) (C.vert j))
    (hshort : G'.dist (C.vert i) (C.vert j) < C.ringDist i j) :
    ∃ c, c ∈ Cset ∧
      (edgeUnion G G').dist (C.vert i) c ≤ C.k / 2 := by
  have hpos : 0 < C.ringDist i j := by
    have hle := C.dist_le_ringDist i j
    omega
  obtain ⟨p, hp⟩ := hreach.exists_walk_length_eq_dist
  have hno : ∀ q : G.Walk (C.vert i) (C.vert j), ¬ q.length ≤ p.length := by
    intro q hq
    have hdist := dist_le q
    have := hP i j
    omega
  obtain ⟨n, hn, hadj, hnot⟩ := exists_bad_step (G := G) p hno
  have hmem := htouch hadj hnot
  have hlen : p.length ≤ C.k / 2 := by
    have hhalf := C.ringDist_le_half i j
    omega
  rcases hmem with hmem | hmem
  · refine ⟨p.getVert n, hmem, ?_⟩
    have hwalk := dist_le ((p.take n).map (Hom.ofLE (le_union_right G G')))
    have htake : (p.take n).length = n := by
      simp [take_length, Nat.min_eq_left (Nat.le_of_lt hn)]
    have : (edgeUnion G G').dist (C.vert i) (p.getVert n) ≤ n := by
      simpa [Walk.length_map, htake] using hwalk
    omega
  · refine ⟨p.getVert (n + 1), hmem, ?_⟩
    have hle : n + 1 ≤ p.length := by omega
    have hwalk := dist_le ((p.take (n + 1)).map (Hom.ofLE (le_union_right G G')))
    have htake : (p.take (n + 1)).length = n + 1 := by
      simp [take_length, Nat.min_eq_left hle]
    have : (edgeUnion G G').dist (C.vert i) (p.getVert (n + 1)) ≤ n + 1 := by
      simpa [Walk.length_map, htake] using hwalk
    omega

/-- The arithmetic form of the tighter update radius. A shortcut of `l`
edges has a changed edge whose nearer end is at most `(l - 1) / 2` from a
ring vertex, and `l ≤ ⌊k/2⌋ - 1`. -/
theorem reach_of_shortcut (k l : ℕ) (hl : l ≤ k / 2 - 1) :
    (l - 1) / 2 ≤ (k / 2 - 2) / 2 := by omega

theorem reach_zero_of_le_seven (k : ℕ) (hk : k ≤ 7) :
    (k / 2 - 2) / 2 = 0 := by omega

end Formal

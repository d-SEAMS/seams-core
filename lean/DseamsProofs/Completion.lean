/-
  The ring completion of the seeded cage assignment, as a closure on a
  finite vertex set.

  MODEL of `ring::ringAdjacentCompletion` (seams-core, cage_affiliation.cpp):
  vertices `V` (a `Fintype`), a finite family `R : Finset (Finset V)` of
  rings, and a labelled set `L`. One step labels every vertex `v` for which
  some ring `r ∈ R` contains `v` and has all its other vertices labelled.
  The C++ walks a frontier and revisits only the rings through a newly
  labelled vertex; the Python prototype passes over all rings once per
  round. These are the same set by the theorems here, not by test:

  * `step` is extensive and monotone;
  * iterating `step` reaches a fixed point (`closure`) within `card V`
    rounds;
  * `closure L` is the least fixed point above `L`, hence independent of
    the visiting order any implementation chooses;
  * `closure ∅ = ∅`: with nothing labelled, nothing is labelled (the
    structural zero on a liquid without a seed);
  * soundness: every vertex the closure adds beyond `L` lies on a ring
    whose other vertices are in the closure (the all-but-one rule as a
    certificate).

  The edge-sharing rule the code first used (accept a ring that shares an
  edge with a fully labelled ring) is stated too, and `edgeRule_not_sound`
  exhibits a three-vertex ring accepted from two labelled vertices, which
  is the mechanism that flooded a liquid with a nucleus in it.

  What this file does NOT establish: that the C++ frontier bookkeeping
  computes `step` (a code property covered by the Catch2 identity, refill
  and structural-zero cases), or anything about which rings a physical
  liquid contains.
-/
import Mathlib.Data.Finset.Basic
import Mathlib.Data.Finset.Card
import Mathlib.Data.Fintype.Basic
import Mathlib.Order.FixedPoints

namespace Dseams

variable {V : Type*} [DecidableEq V] [Fintype V]

/-- A ring `r` fills the vertex `v` under the labelled set `L`: `v ∈ r` and
every other vertex of `r` is labelled. -/
def fills (L : Finset V) (r : Finset V) (v : V) : Prop :=
  v ∈ r ∧ r.erase v ⊆ L

instance (L r : Finset V) (v : V) : Decidable (fills L r v) := by
  unfold fills; infer_instance

/-- One completion round over a ring family `R`. -/
def step (R : Finset (Finset V)) (L : Finset V) : Finset V :=
  L ∪ Finset.univ.filter (fun v => ∃ r ∈ R, fills L r v)

theorem subset_step (R : Finset (Finset V)) (L : Finset V) : L ⊆ step R L :=
  Finset.subset_union_left

theorem fills_mono {L L' : Finset V} (h : L ⊆ L') {r : Finset V} {v : V}
    (hf : fills L r v) : fills L' r v :=
  ⟨hf.1, hf.2.trans h⟩

theorem step_mono (R : Finset (Finset V)) {L L' : Finset V} (h : L ⊆ L') :
    step R L ⊆ step R L' := by
  intro v hv
  rcases Finset.mem_union.mp hv with hL | hF
  · exact Finset.mem_union_left _ (h hL)
  · refine Finset.mem_union_right _ ?_
    rw [Finset.mem_filter] at hF ⊢
    obtain ⟨-, r, hr, hf⟩ := hF
    exact ⟨Finset.mem_univ _, r, hr, fills_mono h hf⟩

/-- `n` rounds of completion. -/
def iter (R : Finset (Finset V)) (L : Finset V) : ℕ → Finset V
  | 0 => L
  | n + 1 => step R (iter R L n)

theorem iter_mono_succ (R : Finset (Finset V)) (L : Finset V) (n : ℕ) :
    iter R L n ⊆ iter R L (n + 1) :=
  subset_step R _

theorem iter_mono (R : Finset (Finset V)) (L : Finset V) {m n : ℕ} (h : m ≤ n) :
    iter R L m ⊆ iter R L n := by
  induction h with
  | refl => exact le_rfl
  | step _ ih => exact ih.trans (iter_mono_succ R L _)

theorem subset_iter (R : Finset (Finset V)) (L : Finset V) (n : ℕ) : L ⊆ iter R L n :=
  iter_mono R L (Nat.zero_le n)

/-- The completion: `card V` rounds, which is enough (`iter_stable`). -/
def closure (R : Finset (Finset V)) (L : Finset V) : Finset V :=
  iter R L (Fintype.card V)

/-- A round that changes nothing changes nothing afterwards. -/
theorem iter_fixed_of_eq (R : Finset (Finset V)) (L : Finset V) {n : ℕ}
    (h : iter R L (n + 1) = iter R L n) : ∀ k, iter R L (n + k) = iter R L n := by
  intro k
  induction k with
  | zero => rfl
  | succ k ih =>
    show step R (iter R L (n + k)) = iter R L n
    rw [ih]; exact h

/-- Strictly growing rounds gain at least one vertex each, so there are at
most `card V` of them. -/
theorem card_iter_ge (R : Finset (Finset V)) (L : Finset V) (n : ℕ)
    (hgrow : ∀ k < n, iter R L (k + 1) ≠ iter R L k) :
    n ≤ (iter R L n).card := by
  induction n with
  | zero => exact Nat.zero_le _
  | succ n ih =>
    have hne : iter R L (n + 1) ≠ iter R L n := hgrow n (Nat.lt_succ_self n)
    have hsub : iter R L n ⊂ iter R L (n + 1) :=
      Finset.ssubset_iff_subset_ne.mpr ⟨iter_mono_succ R L n, fun h => hne h.symm⟩
    have := Finset.card_lt_card hsub
    have ih' := ih (fun k hk => hgrow k (Nat.lt_succ_of_lt hk))
    omega

/-- Some round at or before `card V` is a fixed point. -/
theorem exists_fixed_round (R : Finset (Finset V)) (L : Finset V) :
    ∃ n ≤ Fintype.card V, iter R L (n + 1) = iter R L n := by
  by_contra hcon
  push_neg at hcon
  have h := card_iter_ge R L (Fintype.card V + 1)
    (fun k hk => hcon k (Nat.lt_succ_iff.mp hk))
  have hle : (iter R L (Fintype.card V + 1)).card ≤ Fintype.card V :=
    Finset.card_le_univ _
  omega

theorem iter_stable (R : Finset (Finset V)) (L : Finset V) (k : ℕ) :
    iter R L (Fintype.card V + k) = closure R L := by
  obtain ⟨n, hn, hfix⟩ := exists_fixed_round R L
  unfold closure
  have h1 := iter_fixed_of_eq R L hfix (Fintype.card V + k - n)
  have h2 := iter_fixed_of_eq R L hfix (Fintype.card V - n)
  rw [Nat.add_sub_cancel' (by omega)] at h1
  rw [Nat.add_sub_cancel' hn] at h2
  rw [h1, h2]

/-- The closure is a fixed point of one round. -/
theorem step_closure (R : Finset (Finset V)) (L : Finset V) :
    step R (closure R L) = closure R L := by
  have := iter_stable R L 1
  simpa [iter] using this

theorem subset_closure (R : Finset (Finset V)) (L : Finset V) : L ⊆ closure R L :=
  subset_iter R L _

/-- Every fixed point of `step` above `L` contains the closure: the closure
is the least fixed point, so no implementation order can produce a
different set. -/
theorem closure_le_of_fixed (R : Finset (Finset V)) (L M : Finset V) (hL : L ⊆ M)
    (hM : step R M = M) : closure R L ⊆ M := by
  unfold closure
  generalize Fintype.card V = n
  induction n with
  | zero => exact hL
  | succ n ih =>
    show step R (iter R L n) ⊆ M
    calc step R (iter R L n) ⊆ step R M := step_mono R ih
      _ = M := hM

theorem closure_mono (R : Finset (Finset V)) {L L' : Finset V} (h : L ⊆ L') :
    closure R L ⊆ closure R L' :=
  closure_le_of_fixed R L _ (h.trans (subset_closure R L')) (step_closure R L')

theorem closure_idem (R : Finset (Finset V)) (L : Finset V) :
    closure R (closure R L) = closure R L :=
  Finset.Subset.antisymm
    (closure_le_of_fixed R _ _ le_rfl (step_closure R L))
    (subset_closure R _)

/-- Nothing fills a vertex from the empty labelled set once rings have at
least two vertices (a ring is a cycle of length at least three, so its
other vertices form a nonempty set that cannot lie inside `∅`). -/
theorem step_empty (R : Finset (Finset V)) (hR : ∀ r ∈ R, 2 ≤ r.card) :
    step R (∅ : Finset V) = ∅ := by
  unfold step
  rw [Finset.empty_union, Finset.filter_eq_empty_iff]
  intro v _ ⟨r, hr, hv, hsub⟩
  have hc : (r.erase v).card = r.card - 1 := Finset.card_erase_of_mem hv
  have hpos : 0 < (r.erase v).card := by have := hR r hr; omega
  have : r.erase v = ∅ := Finset.subset_empty.mp hsub
  rw [this, Finset.card_empty] at hpos
  exact Nat.lt_irrefl 0 hpos

/-- Structural zero: no seed, no label. -/
theorem closure_empty (R : Finset (Finset V)) (hR : ∀ r ∈ R, 2 ≤ r.card) :
    closure R (∅ : Finset V) = ∅ := by
  unfold closure
  generalize Fintype.card V = n
  induction n with
  | zero => rfl
  | succ n ih =>
    show step R (iter R ∅ n) = ∅
    rw [show iter R (∅ : Finset V) n = ∅ from ih]
    exact step_empty R hR

/-- Soundness of one round: a vertex added by a round lies on a ring whose
other vertices were labelled before the round. -/
theorem step_sound (R : Finset (Finset V)) (L : Finset V) {v : V}
    (hv : v ∈ step R L) (hnew : v ∉ L) : ∃ r ∈ R, fills L r v := by
  rcases Finset.mem_union.mp hv with h | h
  · exact absurd h hnew
  · exact (Finset.mem_filter.mp h).2

/-- Soundness of the closure: every added vertex sits on a ring whose other
vertices are in the closure. The all-but-one rule is its own certificate. -/
theorem closure_sound (R : Finset (Finset V)) (L : Finset V) {v : V}
    (hv : v ∈ closure R L) (hnew : v ∉ L) : ∃ r ∈ R, fills (closure R L) r v := by
  -- find the first round that contains v
  unfold closure at hv ⊢
  generalize hn : Fintype.card V = n at hv
  induction n with
  | zero => exact absurd hv hnew
  | succ n ih =>
    by_cases hprev : v ∈ iter R L n
    · obtain ⟨r, hr, hf⟩ := ih (by omega) hprev
      exact ⟨r, hr, fills_mono (iter_mono_succ R L n) hf⟩
    · obtain ⟨r, hr, hf⟩ := step_sound R (iter R L n) hv hprev
      exact ⟨r, hr, fills_mono (iter_mono_succ R L n) hf⟩

/-! ### The rule that flooded a liquid

`edgeStep` accepts every vertex of a ring that shares an edge (two vertices)
with a fully labelled ring. It is monotone and extensive as well, but it is
not sound in the sense above: two labelled vertices suffice to accept a
whole ring, whatever its other vertices are. -/

/-- The edge-sharing rule: `r` shares two labelled vertices with some ring
whose vertices are all labelled. -/
def edgeAccepts (R : Finset (Finset V)) (L : Finset V) (r : Finset V) : Prop :=
  ∃ r' ∈ R, r' ⊆ L ∧ 2 ≤ (r ∩ r').card

def edgeStep (R : Finset (Finset V)) (L : Finset V) : Finset V :=
  L ∪ Finset.univ.filter (fun v => ∃ r ∈ R, v ∈ r ∧ edgeAccepts R L r)

/-- Two labelled vertices accept a ring with an unlabelled remainder:
a concrete four-vertex instance. `a b` are labelled and form a ring with
`c`; the ring `{a, b, d}` is then accepted although `d` is unlabelled and
`{a, b, d}` has two unlabelled-free vertices only. -/
theorem edgeRule_not_sound :
    let a : Fin 4 := 0; let b : Fin 4 := 1; let c : Fin 4 := 2; let d : Fin 4 := 3
    let R : Finset (Finset (Fin 4)) := {{a, b, c}, {a, b, d}}
    let L : Finset (Fin 4) := {a, b, c}
    d ∈ edgeStep R L ∧ ¬ ∃ r ∈ R, fills L r d := by
  intro a b c d R L
  refine ⟨?_, ?_⟩
  · -- d is accepted through the ring {a, b, d}, which shares a and b with {a, b, c}
    refine Finset.mem_union_right _ (Finset.mem_filter.mpr ⟨Finset.mem_univ _, {a, b, d}, ?_, ?_, ?_⟩)
    · simp [R]
    · simp [a, b, d]
    · refine ⟨{a, b, c}, by simp [R], by simp [L], ?_⟩
      decide
  · -- but no ring has all its other vertices labelled: {a, b, d} misses d only if d ∈ L, and it is not
    rintro ⟨r, hr, hv, hsub⟩
    simp only [R, Finset.mem_insert, Finset.mem_singleton] at hr
    rcases hr with rfl | rfl
    · simp [a, b, c, d] at hv
    · have : c ∈ r.erase d := by decide
      have := hsub this
      simp [L, a, b, c] at this
      revert this; decide

end Dseams

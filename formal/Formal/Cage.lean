/-
Copyright (c) d-SEAMS developers
Released under the MIT License.
-/

import Mathlib.Combinatorics.SimpleGraph.Basic
import Mathlib.Combinatorics.SimpleGraph.Connectivity.Connected
import Mathlib.Combinatorics.SimpleGraph.Paths
import Mathlib.Data.Finset.Card
import Mathlib.Data.Finset.Powerset
import Mathlib.Tactic

/-!
# Cage search along exposed edges

A cage is grown by adding one face at a time. Each addition shares an edge
that the faces already chosen use only once. The search stops at `F` faces,
so every branch is finite. A connected dual has an exit from every proper
nonempty subset, so every cage is reached from its minimum face. The parent
that deletes the greatest face is unique, so each cage is emitted once.

The tree of all subsets of an `F`-face dual has `2^F` nodes. Walking that
tree does not have delay polynomial in `F`. The canonical parent removes
duplicate discovery; it does not by itself make the delay polynomial.
-/

open Finset

namespace Formal

variable {ι : Type*} [DecidableEq ι]

/-- One growth step adds a single face. -/
def grows (prev next : Finset ι) : Prop :=
  ∃ f, f ∉ prev ∧ next = insert f prev

theorem grows_card (prev next : Finset ι) (h : grows prev next) :
    next.card = prev.card + 1 := by
  rcases h with ⟨f, hf, rfl⟩
  simpa [card_insert_of_notMem hf]

/-- After at most `F` additions the face budget is spent. -/
theorem growth_hits_budget (start steps F : ℕ) (h0 : start = 0) (hsteps : steps ≤ F) :
    start + steps ≤ F := by
  omega

/-- Deleting the greatest element is the only parent of a nonempty set. -/
theorem canonical_parent_unique [LinearOrder ι] (s : Finset ι) (hs : s.Nonempty) :
    ∃! p, ∃ m ∈ s, m = s.max' hs ∧ p = s.erase m := by
  refine ⟨s.erase (s.max' hs), ?_, ?_⟩
  · exact ⟨s.max' hs, s.max'_mem hs, rfl, rfl⟩
  · rintro p ⟨m, hm, hmax, hp⟩
    subst hp
    have : m = s.max' hs := by
      apply le_antisymm
      · exact s.le_max' m hm
      · rw [← hmax]
    simp [this]

/-- The first edge of a walk that starts in `S` and ends outside `S`. -/
theorem walk_exit {G : SimpleGraph ι} {S : Finset ι} {u v : ι}
    (p : G.Walk u v) (hu : u ∈ S) (hv : v ∉ S) :
    ∃ a ∈ S, ∃ b, b ∉ S ∧ G.Adj a b := by
  induction p with
  | nil =>
      exact absurd hu (by simpa using hv)
  | @cons a b c hadj q ih =>
      by_cases hb : b ∈ S
      · exact ih hb hv
      · exact ⟨a, hu, b, hb, hadj⟩

/-- A connected graph offers an edge out of every nonempty proper subset. -/
theorem connected_has_exit {G : SimpleGraph ι} [Fintype ι]
    (hconn : G.Connected) {S : Finset ι} (hS : S.Nonempty) (hproper : ∃ v, v ∉ S) :
    ∃ u ∈ S, ∃ v, v ∉ S ∧ G.Adj u v := by
  obtain ⟨v, hv⟩ := hproper
  obtain ⟨u, hu⟩ := hS
  obtain ⟨p⟩ := hconn u v
  exact walk_exit p hu hv

/-- Every subset of an `F`-element set, counted. This is why a search that
materialises partial face sets does not have delay polynomial in `F`. -/
theorem powerset_exponential (F : ℕ) :
    (Finset.range F).powerset.card = 2 ^ F := by
  simp

end Formal

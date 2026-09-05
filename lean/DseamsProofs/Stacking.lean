/-
  Ice-I stacking planes from TUM rings, not CHILL+ molecules.

  MODEL of `ring::stackingPlanes` and `topoparam::tumLayerStack`:

  * A ring that is the basal of a passing HC pair is a hexagonal plane.
  * A ring that passes the DDC equatorial test is a cubic plane.
  * findDDC refuses HC-affiliated rings as equatorial candidates, so the
    two plane sets are disjoint.
  * Prismatic and peripheral rings, and clathrate faces that are not
    those two predicates, are not planes.
  * A layer letter is a function of the two plane counts in that bin,
    so it does not depend on the order the rings were enumerated.

  What this file does not establish: that the C++ basalConditions test
  matches a physical ice bilayer (Catch2 on the HC fixture covers that).
-/
import Mathlib.Data.Finset.Basic
import Mathlib.Data.Nat.Basic
import Mathlib.Tactic

namespace Dseams.Stacking

variable {V : Type*} [DecidableEq V]

/-- Equatorial candidacy in findDDC: HC-affiliated rings are excluded. -/
def equatorialCandidate (hc : Finset V) (r : V) : Prop := r ∉ hc

theorem basal_equatorial_disjoint
    (basal equatorial hc : Finset V)
    (hBasal : basal ⊆ hc)
    (hEq : ∀ r ∈ equatorial, equatorialCandidate hc r) :
    basal ∩ equatorial = ∅ := by
  ext r
  simp only [Finset.mem_inter, Finset.not_mem_empty, iff_false, not_and]
  intro hb he
  exact hEq r he (hBasal hb)

inductive LayerLetter
  | H
  | C
  | M
  | empty
  deriving DecidableEq, Repr

/-- Majority letter of one layer. Empty when no plane ring votes. -/
def layerLetter (nHex nCubic : ℕ) : LayerLetter :=
  if nHex = 0 ∧ nCubic = 0 then
    LayerLetter.empty
  else if nCubic > nHex then
    LayerLetter.C
  else if nHex > nCubic then
    LayerLetter.H
  else
    LayerLetter.M

theorem layerLetter_empty : layerLetter 0 0 = LayerLetter.empty := rfl

theorem layerLetter_clathrate {nHex nCubic : ℕ} (hH : nHex = 0)
    (hC : nCubic = 0) :
    layerLetter nHex nCubic = LayerLetter.empty := by
  subst hH; subst hC; rfl

/-- A permutation of the votes is the same pair of counts, hence the
    same letter. That is the visit-order independence of tumLayerStack. -/
theorem layerLetter_perm_invariant (nHex nCubic : ℕ)
    (σH σC : ℕ) (hH : σH = nHex) (hC : σC = nCubic) :
    layerLetter σH σC = layerLetter nHex nCubic := by
  subst hH; subst hC; rfl

theorem layerLetter_hex :
    layerLetter 2 1 = LayerLetter.H := by
  native_decide

theorem layerLetter_cubic :
    layerLetter 1 3 = LayerLetter.C := by
  native_decide

end Dseams.Stacking

/-
  Dual occupancy and ion-shell rings. Model of `site::guestOccupancy`,
  `site::guestOccupancyInside`, and `site::shellRings`.
-/

namespace DseamsProofs.Occupancy

def normEdge (a b : Nat) : Nat × Nat := if a ≤ b then (a, b) else (b, a)

def isClosed (faces : List (List Nat)) : Prop :=
  faces ≠ [] ∧ ∀ f ∈ faces, f.length ≥ 3

theorem empty_not_closed : ¬ isClosed [] := by
  intro h
  exact h.1 rfl

def radiusAssigned (dist : Nat → Nat) (r c : Nat) : Prop :=
  dist c ≤ r ∧ ∀ c', dist c ≤ dist c' ∨ r < dist c'

theorem radius_within (dist : Nat → Nat) {r c : Nat}
    (h : radiusAssigned dist r c) : dist c ≤ r :=
  h.1

def insideParity (onFace : Prop) (crossings : Nat) : Prop :=
  onFace ∨ crossings % 2 = 1

theorem on_face_inside (n : Nat) : insideParity True n :=
  Or.inl trivial

theorem odd_inside (n : Nat) : insideParity False (2 * n + 1) := by
  refine Or.inr ?_
  simp

def capped (shell : Nat → Prop) (ring : List Nat) : Prop :=
  ring ≠ [] ∧ ∀ v ∈ ring, shell v

def broken (ion : Nat) (ring : List Nat) : Prop := ion ∈ ring

theorem empty_not_capped (shell : Nat → Prop) : ¬ capped shell [] := by
  intro h
  exact h.1 rfl

end DseamsProofs.Occupancy

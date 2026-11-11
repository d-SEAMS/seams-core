/-
Copyright (c) d-SEAMS developers
Released under the MIT License.
-/

import Formal.Rings
import Mathlib.Combinatorics.SimpleGraph.Metric
import Mathlib.Tactic

/-!
# A ghost skin of `⌊k/2⌋` bonds

A primitive ring of length at most `k` has diameter `⌊k/2⌋`. Every vertex of
the ring therefore lies in the closed ball of that radius about any owned
vertex of the ring. A rank that holds the induced subgraph on that ball sees
the ring, and it sees every strict shortcut as well: a vertex of a shortest
path between two ring vertices stays inside the same ball.
-/

open SimpleGraph

namespace Formal

variable {V : Type*}

theorem min_le_avg (x y : ℕ) : min x y ≤ (x + y) / 2 := by
  rw [Nat.min_def]
  split_ifs <;> omega

/-- A vertex on a shortest `u–v` path of length `l` is at most `(α+β+l)/2`
from a vertex whose distances to the ends are `α` and `β`. -/
theorem shortcut_vertex_bound (a b t l : ℕ) (ht : t ≤ l) :
    min (a + t) (b + (l - t)) ≤ (a + b + l) / 2 := by
  have hmin := min_le_avg (a + t) (b + (l - t))
  have hsum : a + t + (b + (l - t)) = a + b + l := by omega
  simpa [hsum] using hmin

/-- If the ring distance `rd` is at most both arcs out of an owned vertex,
and the cycle has length at most `2L+1`, a shortcut shorter than `rd` stays
inside the closed ball of radius `L`. -/
theorem shortcut_within_half_radius (L a b l rd n : ℕ)
    (hn : n ≤ 2 * L + 1) (hrd : rd + a + b ≤ n) (hl : l + 1 ≤ rd) :
    (a + b + l) / 2 ≤ L := by
  have hsum : a + b + l ≤ 2 * L := by omega
  omega

/-- Every vertex of a primitive ring of length at most `k` is within `⌊k/2⌋`
of every other vertex of the ring. -/
theorem primitive_within_half {G : SimpleGraph V} (C : GraphCycle G)
    (hP : C.Primitive) (i j : Fin C.k) {k : ℕ} (hk : C.k ≤ k) :
    G.Reachable (C.vert i) (C.vert j) ∧ G.dist (C.vert i) (C.vert j) ≤ k / 2 := by
  have hrd := hP i j
  have hhalf := C.ringDist_le_half i j
  have hle : C.ringDist i j ≤ k / 2 := by omega
  refine ⟨?_, ?_⟩
  · exact ⟨(C.arcWalk (C.forward i j) i).copy rfl (C.forward_steps i j)⟩
  · simpa [hrd] using hle

/-- The ghost skin. `dist` of an unreachable pair is the junk value 0, so
membership asks for a real walk. -/
def graphBall (G : SimpleGraph V) (s : V) (L : ℕ) : Set V :=
  {x | G.Reachable s x ∧ G.dist s x ≤ L}

theorem primitive_ring_subset_ball {G : SimpleGraph V} (C : GraphCycle G)
    (hP : C.Primitive) (i : Fin C.k) :
    ∀ j, C.vert j ∈ graphBall G (C.vert i) (C.k / 2) := by
  intro j
  have h := primitive_within_half C hP i j (le_refl C.k)
  simpa [graphBall] using h

/-- Triangle bound along a shortest walk: the vertex `t` steps from `u` is
at most `dist(s,u) + t` from `s`. -/
theorem dist_along_shortcut {G : SimpleGraph V} {s u v : V} (p : G.Walk u v)
    (hp : p.length = G.dist u v) (hsu : G.Reachable s u) (t : ℕ) (ht : t ≤ p.length) :
    G.dist s (p.getVert t) ≤ G.dist s u + t := by
  have hpre := prefix_of_shortest p hp t ht
  have htri := hsu.dist_triangle_left (p.getVert t)
  have hreach : G.dist u (p.getVert t) = t := hpre
  have : G.dist s (p.getVert t) ≤ G.dist s u + G.dist u (p.getVert t) := htri
  simpa [hreach] using this

end Formal

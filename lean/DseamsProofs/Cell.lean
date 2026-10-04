/-
  Restricted-triclinic edge recovery and the one-image ball.

  A LAMMPS bound span is the edge plus the tilt extent:
  `xspan = lx + xmax - xmin`, so `lx = xspan - xmax + xmin`.
  For `a = (lx, 0, 0)`, `b = (xy, ly, 0)`, `c = (xz, yz, lz)` with
  positive `lx, ly, lz`, every nonzero integer combination has
  Euclidean norm at least `min(lx, ly, lz)`. A closed ball of radius
  `c < m / 2` therefore contains at most one point of any coset of the
  lattice. Radius `m / 2` is attained by `p = -m/2` and `p + m`.
-/

import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Tactic

namespace DseamsProofs.Cell

abbrev E := EuclideanSpace ℝ (Fin 3)

/-- `!₂[x, y, z]` is `(WithLp.equiv 2 _).symm ![x, y, z]`. -/
def ofCoords (x y z : ℝ) : E := !₂[x, y, z]

@[simp]
theorem ofCoords_zero (x y z : ℝ) : ofCoords x y z 0 = x := by
  simp [ofCoords]

@[simp]
theorem ofCoords_one (x y z : ℝ) : ofCoords x y z 1 = y := by
  simp [ofCoords]

@[simp]
theorem ofCoords_two (x y z : ℝ) : ofCoords x y z 2 = z := by
  simp [ofCoords]

theorem recovered_edge (lx xmax xmin : ℝ) :
    (lx + xmax - xmin) - xmax + xmin = lx := by
  ring

def xminOf (xy xz : ℝ) : ℝ :=
  min 0 (min xy (min xz (xy + xz)))

def xmaxOf (xy xz : ℝ) : ℝ :=
  max 0 (max xy (max xz (xy + xz)))

theorem xminOf_nonpos (xy xz : ℝ) : xminOf xy xz ≤ 0 := by
  unfold xminOf
  exact min_le_left _ _

theorem xmaxOf_nonneg (xy xz : ℝ) : 0 ≤ xmaxOf xy xz := by
  unfold xmaxOf
  exact le_max_left _ _

theorem span_recovers_edge (lx xy xz : ℝ) :
    let xmin := xminOf xy xz
    let xmax := xmaxOf xy xz
    (lx + xmax - xmin) - xmax + xmin = lx := by
  simpa using recovered_edge lx (xmaxOf xy xz) (xminOf xy xz)

private lemma abs_coord_le_norm (v : E) (i : Fin 3) : |v i| ≤ ‖v‖ := by
  rw [EuclideanSpace.norm_eq]
  have hsum : 0 ≤ ∑ j : Fin 3, ‖v j‖ ^ 2 :=
    Finset.sum_nonneg fun _ _ => sq_nonneg _
  apply (Real.le_sqrt (abs_nonneg _) hsum).2
  calc
    |v i| ^ 2 = ‖v i‖ ^ 2 := by simp [Real.norm_eq_abs]
    _ ≤ ∑ j : Fin 3, ‖v j‖ ^ 2 :=
      Finset.single_le_sum (f := fun j => ‖v j‖ ^ 2)
        (fun _ _ => sq_nonneg _) (Finset.mem_univ i)

private lemma one_le_abs_intCast {n : ℤ} (h : n ≠ 0) : (1 : ℝ) ≤ |(n : ℝ)| := by
  have hnat : (1 : ℤ) ≤ |n| := Int.one_le_abs h
  exact_mod_cast hnat

structure Edges where
  lx : ℝ
  ly : ℝ
  lz : ℝ
  xy : ℝ
  xz : ℝ
  yz : ℝ
  hx : 0 < lx
  hy : 0 < ly
  hz : 0 < lz

def lattice (H : Edges) (na nb nc : ℤ) : E :=
  ofCoords (na * H.lx + nb * H.xy + nc * H.xz) (nb * H.ly + nc * H.yz) (nc * H.lz)

def minEdge (H : Edges) : ℝ :=
  min H.lx (min H.ly H.lz)

theorem minEdge_le_lx (H : Edges) : minEdge H ≤ H.lx :=
  min_le_left _ _

theorem minEdge_le_ly (H : Edges) : minEdge H ≤ H.ly :=
  le_trans (min_le_right _ _) (min_le_left _ _)

theorem minEdge_le_lz (H : Edges) : minEdge H ≤ H.lz :=
  le_trans (min_le_right _ _) (min_le_right _ _)

theorem lattice_norm_ge (H : Edges) (na nb nc : ℤ)
    (h : (na, nb, nc) ≠ (0, 0, 0)) :
    minEdge H ≤ ‖lattice H na nb nc‖ := by
  by_cases hc : nc ≠ 0
  · have habs : (1 : ℝ) ≤ |(nc : ℝ)| := one_le_abs_intCast hc
    have hmul : H.lz ≤ |(nc : ℝ)| * H.lz := by
      simpa using mul_le_mul_of_nonneg_right habs (le_of_lt H.hz)
    have hcoord : |(nc : ℝ) * H.lz| ≤ ‖lattice H na nb nc‖ := by
      simpa [lattice] using abs_coord_le_norm (lattice H na nb nc) 2
    have hz : H.lz ≤ |(nc : ℝ) * H.lz| := by
      rw [abs_mul, abs_of_pos H.hz]
      exact hmul
    exact le_trans (minEdge_le_lz H) (le_trans hz hcoord)
  · push_neg at hc
    by_cases hb : nb ≠ 0
    · have habs : (1 : ℝ) ≤ |(nb : ℝ)| := one_le_abs_intCast hb
      have hmul : H.ly ≤ |(nb : ℝ)| * H.ly := by
        simpa using mul_le_mul_of_nonneg_right habs (le_of_lt H.hy)
      have hraw : |(nb : ℝ) * H.ly + (nc : ℝ) * H.yz| ≤ ‖lattice H na nb nc‖ := by
        simpa [lattice] using abs_coord_le_norm (lattice H na nb nc) 1
      have hcoord : |(nb : ℝ) * H.ly| ≤ ‖lattice H na nb nc‖ := by
        simpa [hc, abs_mul, abs_of_pos H.hy] using hraw
      have hy : H.ly ≤ |(nb : ℝ) * H.ly| := by
        rw [abs_mul, abs_of_pos H.hy]
        exact hmul
      exact le_trans (minEdge_le_ly H) (le_trans hy hcoord)
    · push_neg at hb
      have ha : na ≠ 0 := by
        intro hna
        apply h
        simp [hna, hb, hc]
      have habs : (1 : ℝ) ≤ |(na : ℝ)| := one_le_abs_intCast ha
      have hmul : H.lx ≤ |(na : ℝ)| * H.lx := by
        simpa using mul_le_mul_of_nonneg_right habs (le_of_lt H.hx)
      have hraw :
          |(na : ℝ) * H.lx + (nb : ℝ) * H.xy + (nc : ℝ) * H.xz| ≤
            ‖lattice H na nb nc‖ := by
        simpa [lattice] using abs_coord_le_norm (lattice H na nb nc) 0
      have hcoord : |(na : ℝ) * H.lx| ≤ ‖lattice H na nb nc‖ := by
        simpa [hb, hc, abs_mul, abs_of_pos H.hx] using hraw
      have hx : H.lx ≤ |(na : ℝ) * H.lx| := by
        rw [abs_mul, abs_of_pos H.hx]
        exact hmul
      exact le_trans (minEdge_le_lx H) (le_trans hx hcoord)

theorem norm_sub_le_add (p v : E) : ‖v‖ ≤ ‖p + v‖ + ‖p‖ := by
  calc
    ‖v‖ = ‖(p + v) - p‖ := by simp
    _ ≤ ‖p + v‖ + ‖p‖ := norm_sub_le _ _

theorem ball_unique (m c : ℝ) (hc : c < m / 2) (p v : E)
    (hv : m ≤ ‖v‖) (hp : ‖p‖ ≤ c) (hpv : ‖p + v‖ ≤ c) : False := by
  have htri : ‖v‖ ≤ ‖p + v‖ + ‖p‖ := norm_sub_le_add p v
  have h2 : ‖v‖ ≤ c + c := le_trans htri (add_le_add hpv hp)
  have htwo : ‖v‖ ≤ 2 * c := by
    simpa [two_mul] using h2
  have hlt : 2 * c < m := by
    linarith
  exact lt_irrefl m (lt_of_le_of_lt (le_trans hv htwo) hlt)

private lemma norm_axis (a : ℝ) : ‖ofCoords a 0 0‖ = |a| := by
  rw [EuclideanSpace.norm_eq, Fin.sum_univ_three]
  simp [Real.norm_eq_abs, sq_abs, Real.sqrt_sq_eq_abs]

theorem half_radius_attained (m : ℝ) (hm : 0 ≤ m) :
    ‖ofCoords (-m / 2) 0 0‖ = m / 2 ∧
      ‖ofCoords (-m / 2) 0 0 + ofCoords m 0 0‖ = m / 2 ∧
      ‖ofCoords m 0 0‖ = m := by
  have hsum : ofCoords (-m / 2) 0 0 + ofCoords m 0 0 = ofCoords (m / 2) 0 0 := by
    ext i
    fin_cases i
    · simp; ring
    · simp
    · simp
  refine ⟨?_, ?_, ?_⟩
  · rw [norm_axis, abs_div, abs_neg, abs_of_pos (by norm_num : (0 : ℝ) < 2),
      abs_of_nonneg hm]
  · rw [hsum, norm_axis, abs_of_nonneg (by linarith : (0 : ℝ) ≤ m / 2)]
  · rw [norm_axis, abs_of_nonneg hm]

/-- An index row that leads with the atom contributes that entry to a
raw length and to a same-species count. Coordination drops it. -/
def degree (row : List Nat) (i : Nat) : Nat :=
  (row.filter (· ≠ i)).length

theorem degree_self_header (nbrs : List Nat) (i : Nat)
    (h : ∀ n ∈ nbrs, n ≠ i) :
    degree (i :: nbrs) i = nbrs.length := by
  unfold degree
  rw [List.filter_cons_of_neg (p := fun x => decide (x ≠ i)) (by simp)]
  exact (List.filter_length_eq_length).2 fun a ha => decide_eq_true (h a ha)

end DseamsProofs.Cell

"""Symbolic verification of Horn's absolute orientation algorithm.

Tests the mathematical correctness of the quaternion-based absolute orientation
implementation in absOrientation.cpp against the reference formulas from
Horn (1987), "Closed-form solution of absolute orientation using unit
quaternions", JOSA A 4(4):629-642.

Run with:
    pixi run -e formalanalysis pytest analysis/test_horn_symbolic.py -v
"""

import numpy as np
import pytest
from sympy import (
    Matrix,
    Rational,
    cos,
    det,
    eye,
    expand,
    pi,
    simplify,
    sqrt,
    symbols,
    trigsimp,
)


# ---------------------------------------------------------------------------
# Shared symbolic variables
# ---------------------------------------------------------------------------

Sxx, Sxy, Sxz = symbols("Sxx Sxy Sxz")
Syx, Syy, Syz = symbols("Syx Syy Syz")
Szx, Szy, Szz = symbols("Szx Szy Szz")

q0, qx, qy, qz = symbols("q0 qx qy qz")

# Unit quaternion constraint: q0^2 + qx^2 + qy^2 + qz^2 = 1
_unit_constraint = {q0**2 + qx**2 + qy**2 + qz**2: 1}


def _sub_unit(expr):
    """Substitute the unit-quaternion constraint into *expr*.

    SymPy does not handle the substitution q0^2+qx^2+qy^2+qz^2 -> 1 as a
    single rule reliably, so we eliminate q0^2 = 1 - qx^2 - qy^2 - qz^2 and
    simplify.
    """
    return simplify(expr.subs(q0**2, 1 - qx**2 - qy**2 - qz**2))


def _sub_unit_poly(expr):
    """Reduce *expr* modulo the ideal (q0^2+qx^2+qy^2+qz^2 - 1).

    Repeatedly substitutes q0^2 -> 1 - qx^2 - qy^2 - qz^2 until no
    degree > 1 terms in q0 remain, then does the same for qx, qy, qz
    with respect to the remaining variables.  This is equivalent to
    polynomial reduction modulo the constraint but avoids Poly.div issues.
    """
    e = expand(expr)
    # Iteratively reduce: replace each variable's square using the constraint.
    # q0^2 = 1 - qx^2 - qy^2 - qz^2, etc. -- but we only need one
    # elimination variable to fully reduce degree-4 expressions from R*R^T.
    for _ in range(10):
        e_new = expand(e.subs(q0**2, 1 - qx**2 - qy**2 - qz**2))
        if e_new == e:
            break
        e = e_new
    return simplify(e)


# ---------------------------------------------------------------------------
# Helpers: build the N matrix exactly as the C++ code does
# ---------------------------------------------------------------------------


def _build_N_from_code():
    """Construct the 4x4 N matrix following calcMatrixN (lines 161-207)."""
    N = Matrix.zeros(4, 4)
    # Diagonal
    N[0, 0] = Sxx + Syy + Szz
    N[1, 1] = Sxx - Syy - Szz
    N[2, 2] = -Sxx + Syy - Szz
    N[3, 3] = -Sxx - Syy + Szz
    # Row 0
    N[0, 1] = Syz - Szy
    N[0, 2] = Szx - Sxz
    N[0, 3] = Sxy - Syx
    # Row 1
    N[1, 0] = Syz - Szy
    N[1, 2] = Sxy + Syx
    N[1, 3] = Szx + Sxz
    # Row 2
    N[2, 0] = Szx - Sxz
    N[2, 1] = Sxy + Syx
    N[2, 3] = Syz + Szy
    # Row 3
    N[3, 0] = Sxy - Syx
    N[3, 1] = Szx + Sxz
    N[3, 2] = Syz + Szy
    return N


# ---------------------------------------------------------------------------
# Helpers: rotation matrices from quaternions
# ---------------------------------------------------------------------------


def _R_correct():
    """The standard Hamilton-convention quaternion-to-rotation matrix.

    For unit quaternion q = q0 + qx*i + qy*j + qz*k the rotation matrix is:
        R(0,0) = q0^2 + qx^2 - qy^2 - qz^2
        R(0,1) = 2*(qx*qy - q0*qz)
        R(0,2) = 2*(qx*qz + q0*qy)
        R(1,0) = 2*(qx*qy + q0*qz)
        R(1,1) = q0^2 - qx^2 + qy^2 - qz^2
        R(1,2) = 2*(qy*qz - q0*qx)
        R(2,0) = 2*(qx*qz - q0*qy)
        R(2,1) = 2*(qy*qz + q0*qx)
        R(2,2) = q0^2 - qx^2 - qy^2 + qz^2
    """
    return Matrix(
        [
            [
                q0**2 + qx**2 - qy**2 - qz**2,
                2 * (qx * qy - q0 * qz),
                2 * (qx * qz + q0 * qy),
            ],
            [
                2 * (qx * qy + q0 * qz),
                q0**2 - qx**2 + qy**2 - qz**2,
                2 * (qy * qz - q0 * qx),
            ],
            [
                2 * (qx * qz - q0 * qy),
                2 * (qy * qz + q0 * qx),
                q0**2 - qx**2 - qy**2 + qz**2,
            ],
        ]
    )


def _R_v1():
    """The rotation matrix as quat2RotMatrix shipped it in v1.0.

    The current C++ carries the correct per-component formulas; this
    transcription of the v1.0 release is kept as the defect record the
    verification found:
        R(0,0) = q0*q0 + qx*qx + qy*qy + qz*qz   # == ||q||^2, always 1
        R(0,1) = 2*(qx*qy - q0*qz)
        R(0,2) = 2*(qx*qz + q0*qy)
        R(1,0) = 2*(qy*qx + q0*qz)
        R(1,1) = q0*q0 + qx*qx + qy*qy + qz*qz    # == 1
        R(1,2) = 2*(qy*qz - q0*qy)                 # defect: needs q0*qx
        R(2,0) = 2*(qz*qx - q0*qy)
        R(2,1) = 2*(qz*qy + q0*qz)                 # defect: needs q0*qx
        R(2,2) = q0*q0 + qx*qx + qy*qy + qz*qz    # == 1
    """
    return Matrix(
        [
            [
                q0**2 + qx**2 + qy**2 + qz**2,
                2 * (qx * qy - q0 * qz),
                2 * (qx * qz + q0 * qy),
            ],
            [
                2 * (qy * qx + q0 * qz),
                q0**2 + qx**2 + qy**2 + qz**2,
                2 * (qy * qz - q0 * qy),
            ],
            [
                2 * (qz * qx - q0 * qy),
                2 * (qz * qy + q0 * qz),
                q0**2 + qx**2 + qy**2 + qz**2,
            ],
        ]
    )


# ===================================================================
# Tests
# ===================================================================


class TestNMatrixSymmetry:
    """Verify that the N matrix from calcMatrixN is symmetric."""

    def test_n_matrix_symmetry(self):
        N = _build_N_from_code()
        diff = N - N.T
        for i in range(4):
            for j in range(4):
                assert expand(diff[i, j]) == 0, (
                    f"N[{i},{j}] - N[{j},{i}] = {diff[i, j]} != 0"
                )


class TestCorrectQuatFormula:
    """Verify that the CORRECT quaternion-to-rotation formula is a proper rotation."""

    def test_det_is_one(self):
        """det(R_correct) = 1 for unit quaternions."""
        R = _R_correct()
        d = det(R)
        d_simplified = _sub_unit(d)
        assert simplify(d_simplified - 1) == 0, f"det(R_correct) = {d_simplified}, expected 1"

    def test_orthogonality(self):
        """R_correct * R_correct^T = I for unit quaternions."""
        R = _R_correct()
        product = R * R.T
        for i in range(3):
            for j in range(3):
                entry = _sub_unit_poly(product[i, j])
                expected = 1 if i == j else 0
                assert simplify(entry - expected) == 0, (
                    f"(R*R^T)[{i},{j}] = {entry}, expected {expected}"
                )


class TestV1DefectRecord:
    """The v1.0 formula disagrees with the correct one in exactly five elements.

    This is the defect record: the current C++ carries the correct formula
    (its regression test lives in tests/absor-test.cpp), and this class pins
    what the verification found in the release it audited.
    """

    def test_v1_defect_is_exactly_five_elements(self):
        R_c = _R_correct()
        R_v1 = _R_v1()
        discrepancies = {}
        for i in range(3):
            for j in range(3):
                diff = _sub_unit(R_v1[i, j] - R_c[i, j])
                if simplify(diff) != 0:
                    discrepancies[(i, j)] = str(simplify(diff))
        assert set(discrepancies) == {(0, 0), (1, 1), (1, 2), (2, 1), (2, 2)}, (
            "The v1.0 defect record is the three unit-norm diagonals plus the "
            f"two copy-paste off-diagonals; found {sorted(discrepancies)}"
        )


class TestNumericalBugDemo:
    """Demonstrate the bug with a concrete quaternion.

    90-degree rotation around z-axis:
        q = (cos(pi/4), 0, 0, sin(pi/4)) = (sqrt(2)/2, 0, 0, sqrt(2)/2)

    Expected rotation matrix:
        R = [[0, -1, 0],
             [1,  0, 0],
             [0,  0, 1]]
    """

    _q_vals = {q0: sqrt(2) / 2, qx: 0, qy: 0, qz: sqrt(2) / 2}
    _R_expected = Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])

    def test_correct_formula_gives_right_result(self):
        R = _R_correct().subs(self._q_vals)
        for i in range(3):
            for j in range(3):
                assert simplify(R[i, j] - self._R_expected[i, j]) == 0, (
                    f"R_correct[{i},{j}] = {R[i,j]}, expected {self._R_expected[i,j]}"
                )

    def test_v1_formula_gave_wrong_result(self):
        """The v1.0 formula produces an incorrect matrix for this rotation."""
        R = _R_v1().subs(self._q_vals)
        mismatches = []
        for i in range(3):
            for j in range(3):
                if simplify(R[i, j] - self._R_expected[i, j]) != 0:
                    mismatches.append((i, j))
        # R(0,0) and R(1,1) read 1 where 0 is expected, and R(2,1) picks up
        # 2*q0*qz = 1; the R(1,2) defect vanishes here because qy = 0
        assert sorted(mismatches) == [(0, 0), (1, 1), (2, 1)], (
            f"v1.0 quarter-turn defect signature changed: {mismatches}"
        )


class TestScaleFactorFormula:
    r"""Verify that s = sqrt(sum ||r_ref||^2 / sum ||r_target||^2) minimizes
    the scaled alignment error sum ||r_ref - s * R * r_target||^2.

    For the optimal R (which aligns directions), the objective becomes:
        E(s) = sum ||r_ref||^2 - 2*s*sum(r_ref . R*r_target) + s^2 * sum ||r_target||^2

    Taking dE/ds = 0 gives:
        s_opt = sum(r_ref . R*r_target) / sum ||r_target||^2

    Horn's paper shows that for the optimal quaternion,
        sum(r_ref . R*r_target) = sqrt(sum ||r_ref||^2 * sum ||r_target||^2)

    so that s_opt = sqrt(sum ||r_ref||^2 / sum ||r_target||^2), matching the code.
    """

    def test_scale_factor_formula(self):
        s, A, B, C = symbols("s A B C", positive=True)
        # A = sum ||r_ref||^2, B = sum ||r_target||^2,
        # C = sum(r_ref . R*r_target) at optimal rotation
        # E(s) = A - 2*s*C + s^2*B
        E = A - 2 * s * C + s**2 * B
        dE_ds = E.diff(s)
        # Solve for optimal s
        s_opt = -dE_ds.subs(s, 0)  # coefficient trick won't work; solve directly
        from sympy import solve

        s_solutions = solve(dE_ds, s)
        assert len(s_solutions) == 1
        s_opt = s_solutions[0]
        assert simplify(s_opt - C / B) == 0, f"s_opt = {s_opt}, expected C/B"

        # At optimal rotation, Horn proves C = sqrt(A*B), so s = sqrt(A/B)
        s_horn = s_opt.subs(C, sqrt(A * B))
        assert simplify(s_horn - sqrt(A / B)) == 0, (
            f"s_horn = {s_horn}, expected sqrt(A/B)"
        )


class TestRMSDFormula:
    """Verify RMSD computation with a small numerical example.

    Three points with a known 90-degree rotation around z (no scaling):
        target = [[1,0,0], [0,1,0], [0,0,1]]
        ref    = [[0,1,0], [-1,0,0], [0,0,1]]  (rotated by 90 deg around z)

    Both sets are already centered (centroid = (0, 1/3, 1/3) for target, etc.),
    so we center them first, apply the known rotation, and check RMSD = 0.
    """

    def test_rmsd_formula(self):
        target = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        # 90 deg rotation around z
        R_true = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        ref = (R_true @ target.T).T

        # Center both
        target_c = target - target.mean(axis=0)
        ref_c = ref - ref.mean(axis=0)

        scale = 1.0
        n = len(target_c)
        total_err = 0.0
        for i in range(n):
            rotated = R_true @ target_c[i]
            err_vec = ref_c[i] - scale * rotated
            total_err += np.dot(err_vec, err_vec)
        rmsd = np.sqrt(total_err / n)

        assert rmsd < 1e-14, f"RMSD = {rmsd}, expected ~0 for perfect alignment"

    def test_rmsd_nonzero_with_wrong_rotation(self):
        """RMSD should be nonzero when the wrong rotation is applied."""
        target = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        R_true = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        ref = (R_true @ target.T).T

        target_c = target - target.mean(axis=0)
        ref_c = ref - ref.mean(axis=0)

        # Use identity instead of correct rotation
        R_wrong = np.eye(3)
        scale = 1.0
        n = len(target_c)
        total_err = 0.0
        for i in range(n):
            rotated = R_wrong @ target_c[i]
            err_vec = ref_c[i] - scale * rotated
            total_err += np.dot(err_vec, err_vec)
        rmsd = np.sqrt(total_err / n)

        assert rmsd > 0.1, f"RMSD = {rmsd}, expected nonzero for wrong rotation"

    def test_rmsd_with_scale(self):
        """RMSD should be zero when scale factor is correctly applied."""
        target = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
        ref = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

        target_c = target - target.mean(axis=0)
        ref_c = ref - ref.mean(axis=0)

        # Scale factor: sqrt(sum||ref||^2 / sum||target||^2)
        sum_ref_sq = np.sum(ref_c**2)
        sum_tgt_sq = np.sum(target_c**2)
        scale = np.sqrt(sum_ref_sq / sum_tgt_sq)

        R = np.eye(3)  # No rotation needed
        n = len(target_c)
        total_err = 0.0
        for i in range(n):
            rotated = R @ target_c[i]
            err_vec = ref_c[i] - scale * rotated
            total_err += np.dot(err_vec, err_vec)
        rmsd = np.sqrt(total_err / n)

        assert rmsd < 1e-14, f"RMSD = {rmsd}, expected ~0 with correct scale"
        assert abs(scale - 0.5) < 1e-14, f"scale = {scale}, expected 0.5"

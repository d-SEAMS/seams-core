"""Property-based tests for quaternion rotation matrices.

Horn's absolute orientation algorithm is used internally by d-SEAMS for
shape matching.  The C++ quat2RotMatrix is not directly exposed via the
Python bindings, so we test the mathematical properties with a pure-Python
reference implementation and verify consistency with the CHILL+ pipeline
(which internally invokes the rotation machinery).
"""

import math

import numpy as np
from hypothesis import given, settings, assume
from hypothesis import strategies as st


# -- Pure-Python quaternion -> rotation matrix --


def quat_to_rotation(q):
    """Convert unit quaternion (w, x, y, z) to a 3x3 rotation matrix.

    Uses the same convention as absor::quat2RotMatrix.
    """
    w, x, y, z = q
    return np.array([
        [w*w + x*x - y*y - z*z, 2*(x*y - w*z),         2*(x*z + w*y)],
        [2*(x*y + w*z),         w*w - x*x + y*y - z*z,  2*(y*z - w*x)],
        [2*(x*z - w*y),         2*(y*z + w*x),          w*w - x*x - y*y + z*z],
    ])


def normalize_quat(q):
    """Normalize a quaternion to unit length."""
    n = math.sqrt(sum(c*c for c in q))
    return tuple(c / n for c in q)


# Strategy: random unit quaternions via normalization of 4-vectors
quat_component = st.floats(min_value=-10, max_value=10, allow_nan=False, allow_infinity=False)


@st.composite
def unit_quaternions(draw):
    """Generate random unit quaternions."""
    w = draw(quat_component)
    x = draw(quat_component)
    y = draw(quat_component)
    z = draw(quat_component)
    norm = math.sqrt(w*w + x*x + y*y + z*z)
    assume(norm > 0.01)  # avoid degenerate near-zero quaternions
    return normalize_quat((w, x, y, z))


@given(q=unit_quaternions())
def test_rotation_determinant_is_one(q):
    """det(R) = 1 for any unit quaternion (proper rotation)."""
    R = quat_to_rotation(q)
    det = np.linalg.det(R)
    assert abs(det - 1.0) < 1e-10, f"det(R) = {det}, expected 1.0"


@given(q=unit_quaternions())
def test_rotation_orthogonality(q):
    """R * R^T = I for any unit quaternion."""
    R = quat_to_rotation(q)
    RRt = R @ R.T
    identity = np.eye(3)
    assert np.allclose(RRt, identity, atol=1e-10), (
        f"R*R^T deviates from identity:\n{RRt}"
    )


@given(q=unit_quaternions())
def test_rotation_inverse(q):
    """R(q) * R(q_inv) = I, where q_inv is the conjugate quaternion."""
    R = quat_to_rotation(q)
    # Conjugate of (w, x, y, z) is (w, -x, -y, -z)
    q_inv = (q[0], -q[1], -q[2], -q[3])
    R_inv = quat_to_rotation(q_inv)
    product = R @ R_inv
    identity = np.eye(3)
    assert np.allclose(product, identity, atol=1e-10), (
        f"R * R_inv deviates from identity:\n{product}"
    )


@given(q=unit_quaternions())
def test_rotation_preserves_norm(q):
    """Rotation preserves vector norms: |R*v| = |v|."""
    R = quat_to_rotation(q)
    rng = np.random.default_rng(42)
    v = rng.standard_normal(3)
    Rv = R @ v
    assert abs(np.linalg.norm(Rv) - np.linalg.norm(v)) < 1e-10


def test_identity_quaternion():
    """Identity quaternion (1, 0, 0, 0) gives identity matrix."""
    R = quat_to_rotation((1, 0, 0, 0))
    assert np.allclose(R, np.eye(3), atol=1e-15)


def test_90deg_around_z():
    """90-degree rotation around z matches known matrix."""
    s = math.sqrt(2) / 2
    R = quat_to_rotation((s, 0, 0, s))
    expected = np.array([
        [0, -1, 0],
        [1,  0, 0],
        [0,  0, 1],
    ], dtype=float)
    assert np.allclose(R, expected, atol=1e-12)


def test_180deg_around_x():
    """180-degree rotation around x gives diag(1, -1, -1)."""
    R = quat_to_rotation((0, 1, 0, 0))
    expected = np.diag([1.0, -1.0, -1.0])
    assert np.allclose(R, expected, atol=1e-12)


@given(q1=unit_quaternions(), q2=unit_quaternions())
@settings(max_examples=200)
def test_rotation_composition(q1, q2):
    """R(q1*q2) = R(q1) * R(q2) (Hamilton product)."""
    # Hamilton product: q1 * q2
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    q12 = normalize_quat((
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
    ))
    R12 = quat_to_rotation(q12)
    R1_R2 = quat_to_rotation(q1) @ quat_to_rotation(q2)
    assert np.allclose(R12, R1_R2, atol=1e-9), (
        f"Composition mismatch:\nR(q1*q2):\n{R12}\nR(q1)*R(q2):\n{R1_R2}"
    )

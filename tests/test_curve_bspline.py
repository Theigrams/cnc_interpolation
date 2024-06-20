import numpy as np
import pytest

from core.curve import BSpline


@pytest.fixture
def bspline_curve():
    control_points = np.array([[5.0, 5.0], [10.0, 10.0], [20.0, 15.0], [35.0, 15.0], [45.0, 10.0], [50.0, 5.0]])
    degree = 3
    knots = [0.0, 0.0, 0.0, 0.0, 0.33, 0.66, 1.0, 1.0, 1.0, 1.0]
    curve = BSpline(control_points, degree, knots)
    return curve


@pytest.mark.parametrize(
    "param, res",
    [
        (0.0, (5.0, 5.0)),
        (0.3, (18.617, 13.377)),
        (0.5, (27.645, 14.691)),
        (0.6, (32.143, 14.328)),
        (1.0, (50.0, 5.0)),
    ],
)
def test_evaluate(bspline_curve: BSpline, param, res):
    curve = bspline_curve
    eval_pts = curve.evaluate(param)
    assert np.allclose(eval_pts, res, atol=1e-3)


@pytest.mark.parametrize(
    "param, res",
    [
        (0.0, (45.45454545, 45.45454545)),
        (0.3, (45.26671675, 13.52366642)),
        (0.5, (45.02416338, -0.25500369)),
        (0.6, (44.93369634, -7.00602307)),
        (1.0, (44.11764706, -44.11764706)),
    ],
)
def test_derivative(bspline_curve: BSpline, param, res):
    curve = bspline_curve
    eval_pts = curve.derivative(param)
    assert np.allclose(eval_pts, res, atol=1e-3)


def test_get_length(bspline_curve: BSpline):
    curve = bspline_curve
    assert np.isclose(curve.length, 50.334872752192354)


def test_arc_length_parameterize(bspline_curve: BSpline):
    curve = bspline_curve
    assert np.isclose(curve.arc_parameter(curve.length / 2), 0.5, atol=1e-1)


def test_get_points_arc_equal(bspline_curve: BSpline):
    n = 20
    curve = bspline_curve
    points = curve.get_points(n, arc_equal=True)
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    assert np.allclose(lengths, curve.length / (n - 1), atol=1e-1)

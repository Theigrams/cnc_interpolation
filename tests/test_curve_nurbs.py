import numpy as np
import pytest

from core.curve import NURBS


@pytest.fixture
def circle_curve():
    control_points = np.array(
        [[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1], [1, 0]], dtype=float
    )
    degree = 2
    weights = np.array([1, 0.707, 1, 0.707, 1, 0.707, 1, 0.707, 1])
    knots = np.array([0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1])
    curve = NURBS(control_points, degree, knots, weights)
    return curve


@pytest.mark.parametrize(
    "param, res",
    [
        (0.0, (1.0, 0.0)),
        (0.25, (0.0, 1.0)),
        (0.5, (-1.0, 0.0)),
        (0.75, (0.0, -1.0)),
        (1.0, (1.0, 0.0)),
    ],
)
def test_evaluate(circle_curve: NURBS, param, res):
    eval_pts = circle_curve.evaluate(param)
    assert np.allclose(eval_pts, res, atol=1e-3)


@pytest.mark.parametrize(
    "param, res",
    [
        (0.0, (0.0, 5.656)),
        (0.25, (-5.656, 0.0)),
        (0.5, (0.0, -5.656)),
        (0.75, (5.656, 0.0)),
        (1.0, (0.0, 5.656)),
    ],
)
def test_derivative(circle_curve: NURBS, param, res):
    eval_pts = circle_curve.derivative(param)
    assert np.allclose(eval_pts, res, atol=1e-3)


def test_get_length(circle_curve: NURBS):
    assert np.isclose(circle_curve.length, np.pi * 2, atol=1e-3)


def test_arc_length_parameterize(circle_curve: NURBS):
    assert np.isclose(circle_curve.arc_parameter(circle_curve.length / 2), 0.5, atol=1e-1)


@pytest.fixture
def nurbs_curve():
    control_points = np.array([[5.0, 5.0], [15.0, 25.0], [35.0, 25.0], [45.0, 5.0]])
    weights = np.array([1.0, 2.0, 2.0, 1.0])
    degree = 3
    knots = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0])
    curve = NURBS(control_points, degree, knots, weights)
    return curve


def test_second_derivative(nurbs_curve: NURBS):
    us = [0.0, 0.25, 0.5, 0.75, 1.0]
    eval_pts = nurbs_curve.second_derivative(us)
    expected = np.array(
        [[-240.0, -960.0], [-28.50816, -145.48992], [0.0, -78.36734694], [28.50816, -145.48992], [240.0, -960.0]]
    )
    assert np.allclose(eval_pts, expected, atol=1e-3)


def test_curvature(nurbs_curve: NURBS):
    us = [0.0, 0.25, 0.5, 0.75, 1.0]
    curvatures = nurbs_curve.curvature(us)
    expected = np.array([-0.0119257, -0.05309784, -0.06666667, -0.05309784, -0.0119257])
    assert np.allclose(curvatures, expected, atol=1e-3)

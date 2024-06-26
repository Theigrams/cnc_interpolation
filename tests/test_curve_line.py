import numpy as np
import pytest

from core.curve import Line


def test_evaluate():
    line = Line([0, 0], [1, 1])
    assert np.array_equal(line.evaluate(0.5), np.array([0.5, 0.5]))


def test_evaluate_list():
    line = Line([0, 0], [1, 1])
    assert np.array_equal(line.evaluate([0, 0.5, 1]), np.array([[0, 0], [0.5, 0.5], [1, 1]]))


def test_derivative():
    line = Line([0, 0], [1, 1])
    assert np.array_equal(line.derivative(0), np.array([1, 1]))


def test_get_length():
    line = Line([0, 0], [1, 1])
    assert np.isclose(line.get_length(), np.sqrt(2))


def test_arc_length_parameterize():
    line = Line([0, 0], [1, 1])
    assert np.isclose(line.arc_parameter(1), np.sqrt(2) / 2)

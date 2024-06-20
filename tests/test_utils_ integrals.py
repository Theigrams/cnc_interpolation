import numpy as np
import pytest

from utils.integrals import adaptive_integral


def test_adaptive_integral_simpson():
    def f(x):
        return np.exp(x)

    a, b = 0, 10
    analytic_solution = np.exp(b) - np.exp(a)
    simpson_result = adaptive_integral(f, a, b, 1e-6, method="simpson")
    assert np.isclose(simpson_result, analytic_solution, atol=1e-6)


def test_adaptive_integral_gauss():
    def f(x):
        return np.exp(x)

    a, b = 0, 10
    analytic_solution = np.exp(b) - np.exp(a)
    gauss_result = adaptive_integral(f, a, b, 1e-6, method="gauss")
    assert np.isclose(gauss_result, analytic_solution, atol=1e-6)


def test_adaptive_integral_tolerance():
    def f(x):
        return x**2

    a, b = 0, 1
    analytic_solution = 1 / 3
    result_high_tol = adaptive_integral(f, a, b, 1e-3, method="simpson")
    result_low_tol = adaptive_integral(f, a, b, 1e-6, method="simpson")
    assert result_low_tol - result_high_tol < 1e-3
    assert np.isclose(result_low_tol, analytic_solution, atol=1e-6)


def test_adaptive_integral_intervals():
    def f(x):
        return x**2

    intervals = []
    adaptive_integral(f, 0, 1, 1e-6, method="simpson", intervals=intervals)
    assert len(intervals) > 0

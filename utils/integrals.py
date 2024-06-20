import numpy as np


def simpson(f, a, b):
    c = (a + b) / 2
    return (b - a) / 6 * (f(a) + 4 * f(c) + f(b))


def gauss_3point(f, a, b):
    nodes = np.array([-np.sqrt(3 / 5), 0, np.sqrt(3 / 5)])
    weights = np.array([5 / 9, 8 / 9, 5 / 9])
    transformed_nodes = 0.5 * (b - a) * nodes + 0.5 * (b + a)
    S_ab = 0.5 * (b - a) * np.sum(weights * f(transformed_nodes))
    return S_ab


def adaptive_integral(f, a, b, tol, intervals=None, method="gauss"):
    c = (a + b) / 2
    if method == "simpson":
        integrate = simpson
    else:
        integrate = gauss_3point
    S_ab = integrate(f, a, b)
    S_ac = integrate(f, a, c)
    S_cb = integrate(f, c, b)

    if np.abs(S_ac + S_cb - S_ab) < tol:
        if intervals is not None:
            intervals.append((a, b, S_ac + S_cb))
        return S_ac + S_cb
    else:
        left = adaptive_integral(f, a, c, tol / 2, intervals, method)
        right = adaptive_integral(f, c, b, tol / 2, intervals, method)
        return left + right


if __name__ == "__main__":

    def f(x):
        return np.exp(x)

    # 解析解
    analytic_solution = np.exp(10) - 1

    # 使用辛普森法计算积分
    simpson_intervals = []
    simpson_result = adaptive_integral(f, 0, 10, 1e-7, method="simpson", intervals=simpson_intervals)
    simpson_error = np.abs(simpson_result - analytic_solution)

    # 使用高斯3点法计算积分
    gauss_intervals = []
    gauss_result = adaptive_integral(f, 0, 10, 1e-7, method="gauss", intervals=gauss_intervals)
    gauss_error = np.abs(gauss_result - analytic_solution)

    print(f"解析解: {analytic_solution}")
    print(f"辛普森法结果: {simpson_result}, 误差: {simpson_error}, 区间数: {len(simpson_intervals)}")
    print(f"高斯3点法结果: {gauss_result}, 误差: {gauss_error}, 区间数: {len(gauss_intervals)}")

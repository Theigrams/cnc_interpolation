import numpy as np
from scipy.optimize import newton

from core.toolpath import ToolPath


class BidirectionalScanner:
    def __init__(self, path: ToolPath, dt: float, v_max: float, a_max: float, j_max: float):
        """
        Perform bidirectional scanning to calculate the feedrate limitation.

        Args:
            path (ToolPath): The toolpath to be scanned.
            dt (float): The sampling interval.
            v_max (float): The maximum velocity.
            a_max (float): The maximum acceleration.
            j_max (float): The maximum jerk.
        """
        self.path = path
        self.error = path.chord_error
        self.dt = dt
        self.v_max = v_max
        self.a_max = a_max
        self.j_max = j_max
        self.v_limit = self.path.get_v_limit(dt, self.v_max, self.a_max, j_max=j_max)

    def backward_scan(self):
        self.v_limit[-1] = 0
        for i in range(self.path.N - 1, -1, -1):
            v_je = self.jerk_limit(self.v_limit[i + 1], self.path.L[i])
            v_ae = self.acceleration_limit(self.v_limit[i + 1])
            self.v_limit[i] = min(self.v_limit[i], v_je, v_ae)

    def forward_scan(self):
        self.v_limit[0] = 0
        for i in range(1, self.path.N + 1):
            v_je = self.jerk_limit(self.v_limit[i - 1], self.path.L[i - 1])
            v_ae = self.acceleration_limit(self.v_limit[i - 1])
            self.v_limit[i] = min(self.v_limit[i], v_je, v_ae)

    def acceleration_limit(self, v0):
        """在jerk拉满时, 五阶段双S型曲线能加到的最大速度"""
        v_ae = v0 + self.a_max**2 / self.j_max
        return v_ae

    def jerk_limit(self, v0, L):
        """Huan Zhao, 2013, IJAMT"""

        def f14(v, v0):
            return (v - v0) ** 3 + 4 * v0 * (v - v0) ** 2 + 4 * v0**2 * (v - v0) - L**2 * self.j_max

        def f14_prime(v, v0):
            return 3 * (v - v0) ** 2 + 8 * v0 * (v - v0) + 4 * v0**2

        f = lambda v: f14(v, v0)
        f_prime = lambda v: f14_prime(v, v0)
        v_je = newton(f, v0 + 0.1, f_prime)
        return v_je

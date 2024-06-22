from typing import List, Union

import numpy as np

import utils.geometry as geom
from core.curve import Bezier, Line
from core.toolpath import Block, ToolPath

CurveList = List[Union[Line, Bezier]]


def compute_d2(L, c1, beta, chord_error):
    """
    Compute the approximate optimal d2 values for each line segment, using the adaptive subdivision scheme, to satisfy Eq. (13).

    Args:
        L (N,): Array of shape (N,) containing the lengths of the line segments.
        c1 (float): The maximum curvature.
        beta (N-1,): The angles between two consecutive line segments.
        chord_error (float): The chord error tolerance.
    Returns:
        d2 (N-1,): Approximate optimal d2 values.
    """
    N = len(L)
    # Step 1: Initialize d2 according to Eq. (5).
    d2 = 2 * chord_error / (np.sin(beta) + 1e-6)
    # Step 2: Adjust first and last segments.
    d2[0] = min(d2[0], L[0] / (c1 + 1))
    d2[-1] = min(d2[-1], L[-1] / (c1 + 1))
    # Step 3: Filter and adjust segments.
    for i in range(N - 1):
        # meet the condition d2[i] <= L[i] and d2[i] <= L[i + 1]
        while d2[i] > L[i] or d2[i] > L[i + 1]:
            d2[i] /= 2
    for i in range(N - 3):
        # meet the conditions of Eq. (11)
        while d2[i] + d2[i + 1] > L[i + 1] / (c1 + 1):
            d2[i], d2[i + 1] = d2[i] / 2, d2[i + 1] / 2
    return d2


class SmoothedPath(ToolPath):
    def __init__(self, points, chord_error, c1):
        super().__init__(points, chord_error)
        self.c1 = c1
        self.beta = self.turning_angles / 2
        self.d2 = compute_d2(self.L, self.c1, self.beta, self.chord_error)
        self.blocks = self.generate_blocks()

    def generate_ctrlpts(self, P, d2):
        """
        Generate the control points for the Bezier curves B3 and B4, using Theorem 1.

        Args:
            P (3, dim): Points P0, P1, P2.
            d2 (float): The distance between P1 and Q1.
        Returns:
            ctrlpts_3 (4, dim): Control points for the Bezier curve B3.
            ctrlpts_4 (4, dim): Control points for the Bezier curve B4.
        """
        dir1 = geom.normalize(P[0] - P[1])
        dir2 = geom.normalize(P[2] - P[1])
        Q0 = P[1] + (1 + self.c1) * d2 * dir1
        Q1 = P[1] + d2 * dir1
        Q2 = (P[1] + Q1) / 2
        Q5 = P[1] + d2 * dir2
        Q6 = P[1] + (1 + self.c1) * d2 * dir2
        Q4 = (Q5 + P[1]) / 2
        Q3 = (Q2 + Q4) / 2
        ctrlpts_3 = np.array([Q0, Q1, Q2, Q3])
        ctrlpts_4 = np.array([Q3, Q4, Q5, Q6])
        return ctrlpts_3, ctrlpts_4

    def generate_blocks(self):
        """
        Generate the blocks including the line segments and Bezier curves.

        Returns:
            blocks (N,): Array of shape (N,) containing the blocks for each line segment.
        """
        blocks = []
        bezier_segments: CurveList = []
        for i in range(self.N - 1):
            ctrlpts_3, ctrlpts_4 = self.generate_ctrlpts(self.points[i : i + 3], self.d2[i])
            bezier_segments.append(Bezier(ctrlpts_3))
            bezier_segments.append(Bezier(ctrlpts_4))
        for i in range(self.N):
            if i == 0:
                spline2 = bezier_segments[2 * i]
                line = Line(self.points[i], spline2(0))
                blocks.append(Block([line, spline2]))
            elif i == self.N - 1:
                spline1 = bezier_segments[2 * i - 1]
                line = Line(spline1(1), self.points[i + 1])
                blocks.append(Block([spline1, line]))
            else:
                spline1, spline2 = bezier_segments[2 * i - 1], bezier_segments[2 * i]
                line = Line(spline1(1), spline2(0))
                blocks.append(Block([spline1, line, spline2]))
        return blocks

    def get_v_limit(self, Ts, v_max, a_max, j_max):
        """
        Compute the maximum velocity at each point on the path, using Eq. (16) in Zhao2013.

        Returns:
            v_limit (N+1,): Maximum velocity at each point on the path.
        """
        self.v_max = v_max
        self.a_max = a_max
        self.j_max = j_max
        curvatures = 4 * np.sin(self.beta) / (3 * self.d2 * np.cos(self.beta) ** 2) + 1e-6
        self.curvatures = curvatures
        v_chord = self.chord_error_limit(curvatures, Ts)
        v_curvature = self.curvature_limit(curvatures)
        v_limit = np.minimum(v_chord, v_curvature)
        v_limit = np.concatenate(([0], v_limit, [0]))
        v_limit = np.minimum(v_limit, v_max)
        return v_limit

    def chord_error_limit(self, curvature, Ts):
        """Yeh SS and Hsu PL, 2002, CAD"""
        a = 2 * self.chord_error / curvature - self.chord_error**2
        a = np.maximum(a, 0)
        v_chord = (2 / Ts) * np.sqrt(a)
        return v_chord

    def curvature_limit(self, curvature):
        """Lai JY et al. 2008, IJAMT"""
        v_acc = np.sqrt(self.a_max / curvature)
        v_jerk = (self.j_max / curvature**2) ** (1 / 3)
        return np.minimum(v_acc, v_jerk)


if __name__ == "__main__":
    P = np.array(
        [[200, 100, 200], [100, 200, 200], [0, 200, 100], [0, 100, 0], [100, 0, 0], [200, 0, 100], [200, 100, 200]]
    )

    chord_error = 5
    c1 = 0.5
    smoothed_path = SmoothedPath(P, chord_error, c1)
    print(smoothed_path.d2)

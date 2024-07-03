from typing import List, Union

import numpy as np

import utils.geometry as geom
from core.curve import BSpline, Line
from core.toolpath import Block, ToolPath

CurveList = List[Union[Line, BSpline]]


class CcrPath(ToolPath):
    def __init__(self, points, chord_error):
        super().__init__(points, chord_error)

        self.theta = np.pi - self.turning_angles
        self.alpha = 2 * np.arctan((1 / 3 + np.sin(self.theta / 2)) / np.cos(self.theta / 2))  # Eq. (7)
        self.N0N4_list = self.adjustment_length()
        N3N4_limit1 = 3 * chord_error / np.cos(self.alpha / 2)  # Eq. (4)
        N3N4_limit2 = self.N0N4_list / (4 * np.cos((self.alpha - self.theta) / 2) + 1)  # Eq. (11)
        self.N3N4_list = np.minimum(N3N4_limit1, N3N4_limit2)
        self.blocks = self.generate_blocks()
        self.lengths = np.array([block.length for block in self.blocks])

    def adjustment_length(self):
        L = self.L
        N = len(L)
        L_T = np.zeros(N - 1)
        for i in range(N - 1):
            if i == 0:
                L_T[i] = np.min([L[i], L[i + 1] / 2])
            elif i == N - 2:
                L_T[i] = np.min([L[i] / 2, L[i + 1]])
            else:
                L_T[i] = np.min([L[i] / 2, L[i + 1] / 2])
        return L_T * 0.9

    def generate_ctrlpts(self, P, L, i):
        """
        Generate the control points for the Bspline curve.

        Args:
            P (3, dim): Points P0, P1, P2.
            L (float): The distance between N3 and N4.
            i (int): The index of the line segment.
        Returns:
            ctrlpts (9, dim): Control points for the Bspline curve.
        """
        theta = self.theta[i]
        alpha = self.alpha[i]
        is2D = P.shape[1] == 2
        if is2D:
            P1, P2, P3 = np.hstack((P, np.zeros((3, 1))))
        else:
            P1, P2, P3 = P
        e1 = geom.normalize(P1 - P2)
        e2 = geom.normalize(P3 - P2)
        n = geom.normalize(np.cross(e1, e2))
        rot_angle = (alpha - theta) / 2
        e3 = geom.rodrigues_rotate(e1, n, -rot_angle)
        e4 = geom.rodrigues_rotate(e2, n, rot_angle)
        N4 = P2
        N3 = N4 + L * e3
        N5 = N4 + L * e4
        N1 = N2 = N4 + 4 * L * np.cos((alpha - theta) / 2) * e1
        N6 = N7 = N4 + 4 * L * np.cos((alpha - theta) / 2) * e2
        N0 = N4 + (4 * np.cos((alpha - theta) / 2) + 1) * L * e1
        N8 = N4 + (4 * np.cos((alpha - theta) / 2) + 1) * L * e2
        ctrlpts = np.array([N0, N1, N2, N3, N4, N5, N6, N7, N8])
        if is2D:
            return ctrlpts[:, :2]
        return ctrlpts

    def split_bspline(self, spline: BSpline):
        """
        Split the B-spline curve into two B-spline curves.

        Args:
            spline (BSpline): The B-spline curve.
        Returns:
            spline_1 (BSpline): The first B-spline curve.
            spline_2 (BSpline): The second B-spline curve.
        """
        N0, N1, N2, N3, N4, N5, N6, N7, N8 = spline.control_points
        mid_point = spline(0.5)
        Q1 = np.array([N0, N1, N2, N3, N4 + 1 / 3 * (N3 - N4), mid_point])
        Q2 = np.array([mid_point, N4 + 1 / 3 * (N5 - N4), N5, N6, N7, N8])
        U1 = np.array([0, 0, 0, 0, 1 / 3, 2 / 3, 1, 1, 1, 1])
        U2 = np.array([0, 0, 0, 0, 1 / 3, 2 / 3, 1, 1, 1, 1])
        spline_1 = BSpline(Q1, 3, U1)
        spline_2 = BSpline(Q2, 3, U2)
        return spline_1, spline_2

    def generate_blocks(self):
        """
        Generate the blocks including the line segments and Bezier curves.

        Returns:
            blocks (N,): Array of shape (N,) containing the blocks for each line segment.
        """
        blocks = []
        # Compute transition B-spline
        bsplines: CurveList = []
        U0 = [0, 0, 0, 0, 1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6, 1, 1, 1, 1]
        for i in range(self.N - 1):
            ctrlpts = self.generate_ctrlpts(self.points[i : i + 3], self.N3N4_list[i], i)
            bsplines.append(BSpline(ctrlpts, 3, U0))
        self.bsplines = bsplines

        # Compute curvature peaks
        self.curvature_peaks = np.array([np.abs(spline.curvature(0.5)) for spline in bsplines])

        # Compute split B-splines
        bspline_segments: CurveList = []
        for i in range(self.N - 1):
            # spline_1, spline_2 = self.split_bspline(bsplines[i])
            spline_1, spline_2 = bsplines[i].split(0.5)
            bspline_segments.append(spline_1)
            bspline_segments.append(spline_2)

        # Generate blocks
        for i in range(self.N):
            if i == 0:
                spline_2 = bspline_segments[2 * i]
                line = Line(self.points[i], spline_2(0))
                blocks.append(Block([line, spline_2]))
            elif i == self.N - 1:
                spline_1 = bspline_segments[2 * i - 1]
                line = Line(spline_1(1), self.points[i + 1])
                blocks.append(Block([spline_1, line]))
            else:
                spline_1, spline_2 = bspline_segments[2 * i - 1], bspline_segments[2 * i]
                line = Line(spline_1(1), spline_2(0))
                blocks.append(Block([spline_1, line, spline_2]))
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
        v_chord = self.chord_error_limit(self.curvature_peaks, Ts)
        v_curvature = self.curvature_limit(self.curvature_peaks)
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
    ccr_path = CcrPath(P, chord_error)
    print(ccr_path.lengths)

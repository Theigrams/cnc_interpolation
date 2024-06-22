from typing import List, Union

import numpy as np

import utils.geometry as geom
from core.curve import Bezier, BSpline, Line

CurveList = List[Union[Line, Bezier, BSpline]]


class Block:
    """进给率规划单元"""

    def __init__(self, curves: CurveList):
        self.curves = curves
        self.lengths = np.array([curve.length for curve in curves])
        self.cum_lengths = np.cumsum(self.lengths)
        self.length = self.cum_lengths[-1]


class ToolPath:
    def __init__(self, points, chord_error):
        self.points = points
        self.chord_error = chord_error
        self.N = len(points) - 1
        self.tangents, self.L = geom.compute_tangents(points)
        self.turning_angles = geom.compute_turning_angles(self.tangents)
        self.blocks: List[Block] = []

    def generate_blocks(self):
        raise NotImplementedError

    def get_v_limit(self, Ts, v_max, a_max, **kwargs):
        """
        Compute the maximum velocity at each point on the path.

        Returns:
            v_limit (N+1,): Maximum velocity at each point on the path.
        """
        raise NotImplementedError


class LinearPath(ToolPath):
    def __init__(self, points, chord_error):
        super().__init__(points, chord_error)
        self.beta = self.turning_angles / 2
        self.blocks = self.generate_blocks()

    def generate_blocks(self):
        blocks = []
        for i in range(self.N):
            line = Line(self.points[i], self.points[i + 1])
            blocks.append(Block([line]))
        return blocks

    def get_v_limit(self, Ts, v_max, a_max, **kwargs):
        """
        Compute the maximum velocity at each point on the path, using Eq. (25) in Zhao2013.

        Returns:
            v_limit (N+1,): Maximum velocity at each point on the path.
        """
        v_limit = a_max * Ts / (2 * np.sin(self.beta) + 1e-6)
        v_limit = np.concatenate(([0], v_limit, [0]))
        v_limit = np.minimum(v_limit, v_max)
        return v_limit

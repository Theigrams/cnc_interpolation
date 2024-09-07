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

    def get_u_from_length(self, s: float):
        """根据曲线段长度, 返回参数 u 和对应的曲线段索引"""
        idx = np.searchsorted(self.cum_lengths, s, side="left")
        if idx >= len(self.cum_lengths):
            raise ValueError(f"Invalid length: {s}, exceeds total length {self.length}")
        arc_length = s - self.cum_lengths[idx - 1] if idx > 0 else s
        u = self.curves[idx].get_u_from_length(arc_length)
        return u, idx

    def __str__(self):
        return f"Block(curves={[str(curve) for curve in self.curves]})"

    def __repr__(self):
        return f"Block(curves={self.curves})"


class ToolPath:
    def __init__(self, points, chord_error):
        self.points = points
        self.chord_error = chord_error
        self.N = len(points) - 1
        self.tangents, self.L = geom.compute_tangents(points)
        self.turning_angles = geom.compute_turning_angles(self.tangents)
        self.blocks: List[Block] = []
        self.lengths = self.L

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
        self.lengths = np.array([block.length for block in self.blocks])

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


if __name__ == "__main__":
    from papers.Zhao2013.algorithm import SmoothedPath

    points = np.array([[2, 0], [4, 2], [2, 4], [0, 2], [2, 0]])
    chord_error = 0.2
    smoothed_path = SmoothedPath(points, chord_error, 0.5)
    for block in smoothed_path.blocks:
        print(block)

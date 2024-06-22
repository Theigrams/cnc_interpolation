from typing import List

import matplotlib.pyplot as plt
import numpy as np

from core.curve import Line
from core.toolpath import Block, LinearPath, ToolPath


def plot_path(path: LinearPath, save=False):
    plt.figure()
    plt.plot(*path.points.T, "k-")
    plt.plot(*path.points.T, "ro")
    plt.axis("equal")
    if save:
        plt.savefig("path.pdf", bbox_inches="tight")
    plt.show()


def plot_toolpath(path: ToolPath, save=False):
    plt.figure()
    blocks: List[Block] = path.blocks
    for block in blocks:
        for curve in block.curves:
            pts = curve.get_points(30)
            style = "b-" if isinstance(curve, Line) else "r-"
            plt.plot(pts[:, 0], pts[:, 1], style)
    plt.axis("equal")
    if save:
        plt.savefig("toolpath.pdf", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    from papers.Zhao2013.algorithm import SmoothedPath

    points = np.array([[2, 0], [4, 2], [2, 4], [0, 2], [2, 0]])
    chord_error = 0.2
    linear_path = LinearPath(points, chord_error)
    plot_path(linear_path)
    smoothed_path = SmoothedPath(points, chord_error, 0.5)
    plot_toolpath(smoothed_path)

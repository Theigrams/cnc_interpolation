import numpy as np

from core.feedrate_profile import FivePhaseProfile as FeedrateProfile
from core.feedrate_scheduler import FeedrateScheduler
from core.toolpath import Block, ToolPath


class Interpolator:
    def __init__(self, path: ToolPath, scheduler: FeedrateScheduler, Ts: float):
        self.path = path
        self.scheduler = scheduler
        self.Ts = Ts
        self.N = len(path.blocks)
        self.blocks = path.blocks
        self.profiles = self.scheduler.profiles

    def interpolate(self):
        interpolated_points = []
        remaining_time = 0
        for i in range(self.N):
            block = self.blocks[i]
            profile = self.profiles[i]
            block_points, remaining_time = self.interpolate_block(block, profile, remaining_time)
            interpolated_points += block_points
        return np.array(interpolated_points)

    def interpolate_block(self, block: Block, profile: FeedrateProfile, remaining_time: float):
        total_time = profile.total_time
        interpolated_points = []
        t = remaining_time
        while t < total_time:
            s, v, a, j = profile.get_motion_state(t)
            u, idx = block.get_u_from_length(s)
            point = block.curves[idx](u)
            interpolated_points.append(point)
            t += self.Ts
        remaining_time = t - total_time

        return interpolated_points, remaining_time


if __name__ == "__main__":
    import utils.visualization as vis
    from core.look_ahead import BidirectionalScanner
    from papers.Zhao2013.algorithm import SmoothedPath

    Ts = 0.0005
    V_MAX = 100
    A_MAX = 3000
    J_MAX = 60000
    chord_error = 0.2
    points = np.array([[2, 0], [4, 2], [2, 4], [0, 2], [2, 0]])

    smooth_path = SmoothedPath(points, chord_error, 0.5)
    vis.plot_toolpath(smooth_path)
    scanner_smooth = BidirectionalScanner(smooth_path, Ts, V_MAX, A_MAX, J_MAX)
    smooth_path_lengths = [block.length for block in smooth_path.blocks]
    v_lim_smooth = scanner_smooth.v_limit
    smooth_scheduler = FeedrateScheduler(smooth_path_lengths, v_lim_smooth, V_MAX, A_MAX, J_MAX, Ts)
    smooth_profiles = smooth_scheduler.profiles
    interpolator = Interpolator(smooth_path, smooth_scheduler, Ts)
    interpolated_points = interpolator.interpolate()
    print(len(interpolated_points))

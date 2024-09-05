from typing import List

import numpy as np

from core.feedrate_profiles import FivePhaseProfile as FeedrateProfile


class FeedrateScheduler:
    def __init__(self, seg_lengths, v_lim, v_max, a_max, j_max, Ts):
        self.seg_lengths = seg_lengths
        self.v_lim = v_lim
        self.v_max = v_max
        self.a_max = a_max
        self.j_max = j_max
        self.Ts = Ts
        self.profiles: List[FeedrateProfile] = []
        self.schedule()
        self.total_time = np.sum([profile.total_time for profile in self.profiles])

    def schedule(self):
        N = len(self.seg_lengths)
        for i in range(N):
            v_str = self.v_lim[i]
            v_end = self.v_lim[i + 1]
            seg_length = self.seg_lengths[i]
            profile = FeedrateProfile(v_str, v_end, seg_length, self.v_max, self.a_max, self.j_max, self.Ts)
            profile.generate_profile()
            self.profiles.append(profile)

    def get_profile_data(self, dt=0.01):
        dtype = [("T", "f8"), ("S", "f8"), ("V", "f8"), ("A", "f8"), ("J", "f8")]
        profile_data = np.empty((0,), dtype=dtype)
        t_curr = 0
        for profile in self.profiles:
            ts = np.arange(0, profile.total_time, dt)
            for t in ts:
                s, v, a, j = profile.get_motion_state(t)
                motion_data = np.array([(t + t_curr, s, v, a, j)], dtype=dtype)
                profile_data = np.append(profile_data, motion_data)
            t_curr += profile.total_time
        return profile_data

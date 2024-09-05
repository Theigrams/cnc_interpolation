import numpy as np


def double_s_trajectory(h, v0, v1, v_max, a_max, j_max):
    # 初始化输出参数
    v_max_reached, a_max_reached, a_min_reached = False, False, False

    # 判断轨迹可行性
    T_j_star = min(np.sqrt(np.abs(v1 - v0) / j_max), a_max / j_max)
    if T_j_star < a_max / j_max:
        if h < T_j_star * (v0 + v1):
            print(f"h: {h}, T_j_star: {T_j_star}, v0: {v0}, v1: {v1}")
            raise ValueError("Trajectory is not feasible: h > T_j_star * (v0 + v1)")
    else:
        if h < (v0 + v1) * (T_j_star + np.abs(v1 - v0) / a_max) / 2:
            raise ValueError("Trajectory is not feasible: h > (v0 + v1) * (T_j_star + np.abs(v1 - v0) / a_max) / 2")

    # 假设v_max可达
    a_max_reached = (v_max - v0) * j_max >= a_max**2
    if a_max_reached:
        # a_max可达, 使用式(3.22)
        T_j1 = a_max / j_max
        T_a = T_j1 + (v_max - v0) / a_max
    else:
        # a_max不可达，使用式(3.21)
        T_j1 = np.sqrt(np.abs(v_max - v0) / j_max)
        T_a = 2 * T_j1

    a_min_reached = (v_max - v1) * j_max >= a_max**2
    if a_min_reached:
        # a_max可达, 使用式(3.24)
        T_j2 = a_max / j_max
        T_d = T_j2 + (v_max - v1) / a_max
    else:
        # a_max不可达，使用式(3.23)
        T_j2 = np.sqrt(np.abs(v_max - v1) / j_max)
        T_d = 2 * T_j2
    # 计算恒速段时间, 使用式(3.25)
    h_acc = h - (v0 + v_max) / 2 * T_a - (v1 + v_max) / 2 * T_d
    T_v = h_acc / v_max

    if T_v > 0:
        # v_max可以达到
        v_max_reached = True
    else:
        # v_max不可达，使用公式(3.26a),(3.26b),(3.26c)和(3.27)计算T_a和T_d
        v_max_reached = False
        T_v = 0
        for _ in range(200):
            delta_v = (
                (a_max**4 / j_max**2) + 2 * (v0**2 + v1**2) + a_max * (4 * h - 2 * (a_max / j_max) * (v0 + v1))
            )  # (3.27)
            T_j1 = T_j2 = a_max / j_max  # (3.26a)
            T_a = (a_max**2 / j_max - 2 * v0 + np.sqrt(delta_v)) / (2 * a_max)  # (3.26b)
            T_d = (a_max**2 / j_max - 2 * v1 + np.sqrt(delta_v)) / (2 * a_max)  # (3.26c)

            # 检查a_max是否可达
            if T_a > 2 * T_j1 and T_d > 2 * T_j2:
                break

            a_max *= 0.99
            if a_max <= 1e-6:
                break
        if T_a < 0 or T_d < 0:
            if v0 > v1:
                T_d = 2 * h / (v1 + v0)  # (3.28a)
                T_j2 = (j_max * h - np.sqrt(j_max * (j_max * h**2 + (v1 + v0) ** 2 * (v1 - v0)))) / (
                    j_max * (v1 + v0)
                )  # (3.28b)
                T_a, T_j1 = 0, 0
            else:
                T_a = 2 * h / (v1 + v0)  # (3.29a)
                T_j1 = (j_max * h - np.sqrt(j_max * (j_max * h**2 - (v1 + v0) ** 2 * (v1 - v0)))) / (
                    j_max * (v1 + v0)
                )  # (3.29b)
                T_d, T_j2 = 0, 0

    a_lim_a = j_max * T_j1
    a_lim_d = -j_max * T_j2
    v_lim = v0 + a_lim_a * (T_a - T_j1)
    a_max_reached = a_max - a_lim_a < 1e-6
    a_min_reached = a_max + a_lim_d < 1e-6

    return {
        "T_j1": T_j1,
        "T_j2": T_j2,
        "T_a1": T_a - 2 * T_j1,
        "T_a2": T_d - 2 * T_j2,
        "T_a": T_a,
        "T_v": T_v,
        "T_d": T_d,
        "T": T_a + T_v + T_d,
        "a_max_reached": a_max_reached,
        "a_min_reached": a_min_reached,
        "v_max_reached": v_max_reached,
        "a_lim_a": a_lim_a,
        "a_lim_d": a_lim_d,
        "v_lim": v_lim,
    }


class SevenPhaseProfile:
    def __init__(self, v_str, v_end, arc_length, v_max, a_max, j_max, Ts):
        self.v_str = v_str
        self.v_end = v_end
        self.arc_length = arc_length
        self.v_max = v_max
        self.a_max = a_max
        self.j_max = j_max
        self.Ts = Ts

    def generate_profile(self):
        result = double_s_trajectory(self.arc_length, self.v_str, self.v_end, self.v_max, self.a_max, self.j_max)
        self.jerk_profile = np.array([self.j_max, 0, -self.j_max, 0, -self.j_max, 0, self.j_max])
        self.time_profile = np.array(
            [
                result["T_j1"],
                result["T_a1"],
                result["T_j1"],
                result["T_v"],
                result["T_j2"],
                result["T_a2"],
                result["T_j2"],
            ]
        )
        self.total_time = np.sum(self.time_profile)
        self.result = result

    def get_motion_state(self, t):
        T_j1, T_a1, _, T_v, T_j2, T_a2, _ = self.time_profile
        T_a = 2 * T_j1 + T_a1
        T = self.total_time

        self.a_lim_a, self.a_lim_d = self.result["a_lim_a"], self.result["a_lim_d"]
        self.v_lim = self.result["v_lim"]

        if 0 <= t < T_j1:
            return self.run_phase_1(t)
        elif T_j1 <= t < T_j1 + T_a1:
            return self.run_phase_2(t)
        elif T_j1 + T_a1 <= t < T_a:
            return self.run_phase_3(t)
        elif T_a <= t < T_a + T_v:
            return self.run_phase_4(t)
        elif T_a + T_v <= t < T_a + T_v + T_j2:
            return self.run_phase_5(t)
        elif T_a + T_v + T_j2 <= t < T_a + T_v + T_j2 + T_a2:
            return self.run_phase_6(t)
        elif T - T_j2 <= t <= T:
            return self.run_phase_7(t)
        else:
            raise ValueError("Invalid time value")

    def run_phase_1(self, t):
        J, v_str = self.j_max, self.v_str
        a = J * t
        v = v_str + J * t**2 / 2
        s = self.v_str * t + J * t**3 / 6
        return s, v, a, J

    def run_phase_2(self, t):
        J, v_str, a_lim_a = self.j_max, self.v_str, self.a_lim_a
        T_j1 = self.time_profile[0]
        a = a_lim_a
        v = v_str + a_lim_a * (t - T_j1 / 2)
        s = self.v_str * t + a_lim_a / 6 * (3 * t**2 - 3 * T_j1 * t + T_j1**2)
        return s, v, a, 0

    def run_phase_3(self, t):
        J, v_lim = self.j_max, self.v_lim
        T_a = self.time_profile[0] + self.time_profile[1] + self.time_profile[2]
        a = J * (T_a - t)
        v = v_lim - J * (T_a - t) ** 2 / 2
        s = (self.v_str + v_lim) * T_a / 2 - v_lim * (T_a - t) + J * (T_a - t) ** 3 / 6
        return s, v, a, -J

    def run_phase_4(self, t):
        v_lim = self.v_lim
        T_a = sum(self.time_profile[:3])
        s = (self.v_str + v_lim) * T_a / 2 + v_lim * (t - T_a)
        return s, v_lim, 0, 0

    def run_phase_5(self, t):
        J, v_lim = self.j_max, self.v_lim
        T = self.total_time
        T_d = self.time_profile[4] + self.time_profile[5] + self.time_profile[6]
        a = -J * (t - T + T_d)
        v = v_lim - J * (t - T + T_d) ** 2 / 2
        s = self.arc_length - (v_lim + self.v_end) * T_d / 2 + v_lim * (t - T + T_d) - J * (t - T + T_d) ** 3 / 6
        return s, v, a, -J

    def run_phase_6(self, t):
        J, v_lim, a_lim_d = self.j_max, self.v_lim, self.a_lim_d
        T = self.total_time
        T_d = self.time_profile[4] + self.time_profile[5] + self.time_profile[6]
        T_j2 = self.time_profile[4]
        a = a_lim_d
        v = v_lim + a_lim_d * (t - T + T_d - T_j2 / 2)
        s = (
            self.arc_length
            - (v_lim + self.v_end) * T_d / 2
            + v_lim * (t - T + T_d)
            + a_lim_d / 6 * (3 * (t - T + T_d) ** 2 - 3 * T_j2 * (t - T + T_d) + T_j2**2)
        )
        return s, v, a, 0

    def run_phase_7(self, t):
        J, v_end = self.j_max, self.v_end
        T = self.total_time
        a = -J * (T - t)
        v = v_end + J * (T - t) ** 2 / 2
        s = self.arc_length - v_end * (T - t) - J * (T - t) ** 3 / 6
        return s, v, a, J


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    p = SevenPhaseProfile(1, 0, 10, 5, 10, 30, 0.002)
    p.generate_profile()
    total_t = p.total_time
    S, V, A, J = [], [], [], []
    T = np.linspace(0, total_t, 1000)
    for t in T:
        s, v, a, j = p.get_motion_state(t)
        S.append(s)
        V.append(v)
        A.append(a)
        J.append(j)

    plt.figure(figsize=(10, 8))
    plt.subplot(2, 2, 1)
    plt.plot(T, S)
    plt.xlabel("Time (s)")
    plt.ylabel("Position (mm)")
    plt.title("Position vs. Time")
    plt.grid(True)

    plt.subplot(2, 2, 2)
    plt.plot(T, V)
    plt.xlabel("Time (s)")
    plt.ylabel("Velocity (mm/s)")
    plt.title("Velocity vs. Time")
    plt.grid(True)

    plt.subplot(2, 2, 3)
    plt.plot(T, A)
    plt.xlabel("Time (s)")
    plt.ylabel("Acceleration (mm/s^2)")
    plt.title("Acceleration vs. Time")
    plt.grid(True)

    plt.subplot(2, 2, 4)
    plt.plot(T, J)
    plt.xlabel("Time (s)")
    plt.ylabel("Jerk (mm/s^3)")
    plt.title("Jerk vs. Time")
    plt.grid(True)

    plt.tight_layout()
    plt.show()

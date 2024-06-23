import numpy as np
from scipy.optimize import newton


class FivePhaseProfile:
    ACC_TYPES = {
        "CV": 1,
        "ACC": 2,
        "DEC": 3,
        "ACC_CV": 4,
        "CV_DEC": 5,
        "ACC_DEC": 6,
        "ACC_CV_DEC": 7,
    }

    def __init__(self, v_str, v_end, arc_length, v_max, a_max, j_max, Ts):
        self.v_str = v_str
        self.v_end = v_end
        self.arc_length = arc_length
        self.v_max = v_max
        self.a_max = a_max
        self.j_max = j_max
        self.Ts = Ts  # 插补周期

    def generate_profile(self):
        self.acc_type = self.determine_type(self.v_str, self.v_end, self.arc_length)
        self.jerk_profile = np.array([self.j_max, -self.j_max, 0, -self.j_max, self.j_max])
        self.time_profile = self.get_profile(self.acc_type, self.v_str, self.v_end, self.arc_length)
        self.total_time = np.sum(self.time_profile)

    def determine_type(self, v_str, v_end, arc_length):
        if v_str > self.v_max or v_end > self.v_max:
            raise Exception("Invalid velocity values")

        J = self.j_max
        L_r1 = 2 * v_str * self.Ts
        L_r2 = (v_str + v_end) * np.sqrt(np.abs(v_end - v_str) / J)
        L_r3 = (v_str + self.v_max) * np.sqrt((self.v_max - v_str) / J) + (self.v_max + v_end) * np.sqrt(
            (self.v_max - v_end) / J
        )
        if v_str < v_end:
            L_r4 = (v_str + v_end) * np.sqrt((v_end - v_str) / J)
        else:
            L_r4 = (v_end + self.v_max) * np.sqrt((self.v_max - v_end) / J)

        if v_str == self.v_max and v_end == self.v_max:
            return self.ACC_TYPES["CV"]
        elif v_str < self.v_max and v_end < self.v_max:
            if arc_length <= L_r1:
                return self.ACC_TYPES["CV"]
            elif arc_length <= L_r2:
                if v_str < v_end:
                    return self.ACC_TYPES["ACC"]
                else:
                    return self.ACC_TYPES["DEC"]
            elif arc_length <= L_r3:
                return self.ACC_TYPES["ACC_DEC"]
            else:
                return self.ACC_TYPES["ACC_CV_DEC"]
        elif v_str == self.v_max or v_end == self.v_max:
            if arc_length <= L_r1:
                return self.ACC_TYPES["CV"]
            elif arc_length <= L_r4:
                if v_str <= v_end:
                    return self.ACC_TYPES["ACC"]
                else:
                    return self.ACC_TYPES["DEC"]
            else:
                if v_str <= v_end:
                    return self.ACC_TYPES["ACC_CV"]
                else:
                    return self.ACC_TYPES["CV_DEC"]
        else:
            raise Exception("Invalid velocity values")

    def get_profile(self, acc_type, v_str, v_end, arc_length):
        if acc_type == self.ACC_TYPES["CV"]:
            return self._cv_profile(v_str, v_end, arc_length)
        elif acc_type == self.ACC_TYPES["ACC"]:
            return self._acc_profile(v_str, v_end, arc_length)
        elif acc_type == self.ACC_TYPES["DEC"]:
            return self._dec_profile(v_str, v_end, arc_length)
        elif acc_type == self.ACC_TYPES["ACC_CV"]:
            return self._acc_cv_profile(v_str, v_end, arc_length)
        elif acc_type == self.ACC_TYPES["CV_DEC"]:
            return self._cv_dec_profile(v_str, v_end, arc_length)
        elif acc_type == self.ACC_TYPES["ACC_DEC"]:
            return self._acc_dec_profile(v_str, v_end, arc_length)
        elif acc_type == self.ACC_TYPES["ACC_CV_DEC"]:
            return self._acc_cv_dec_profile(v_str, v_end, arc_length)
        else:
            raise Exception("Invalid acceleration type")

    def get_motion_state(self, t):
        t_str, _, t_c, t_end, _ = self.time_profile
        if t <= t_str:
            return self.run_phase_1(t)
        elif t <= 2 * t_str:
            return self.run_phase_2(t - t_str, t_str)
        elif t <= 2 * t_str + t_c:
            return self.run_phase_3(t - 2 * t_str, t_str)
        elif t <= 2 * t_str + t_c + t_end:
            return self.run_phase_4(t - 2 * t_str - t_c, t_str, t_c)
        elif t <= 2 * t_str + t_c + 2 * t_end:
            return self.run_phase_5(t - 2 * t_str - t_c - t_end, t_str, t_c, t_end)
        else:
            raise Exception("Invalid time value")

    def _cv_profile(self, v_str, v_end, arc_length):
        t_c = 2 * arc_length / (v_str + v_end)
        time_profile = np.array([0, 0, t_c, 0, 0])
        return time_profile

    def _acc_profile(self, v_str, v_end, arc_length):

        def f_type2(t):
            return self.j_max * t**3 + 2 * v_str * t - arc_length

        def f_prime_type2(t):
            return 3 * self.j_max * t**2 + 2 * v_str

        t_guess = arc_length / (v_str + v_end)
        t_str = newton(f_type2, t_guess, f_prime_type2)
        time_profile = np.array([t_str, t_str, 0, 0, 0])
        return time_profile

    def _dec_profile(self, v_str, v_end, arc_length):

        def f_type3(t):
            return self.j_max * t**3 + 2 * v_end * t - arc_length

        def f_prime_type3(t):
            return 3 * self.j_max * t**2 + 2 * v_end

        t_guess = arc_length / (v_str + v_end)
        t_end = newton(f_type3, t_guess, f_prime_type3)
        time_profile = np.array([0, 0, 0, t_end, t_end])
        return time_profile

    def _acc_cv_profile(self, v_str, v_end, arc_length):
        t_str = np.sqrt((v_end - v_str) / self.j_max)
        L_r4 = (v_str + v_end) * t_str
        t_c = (arc_length - L_r4) / v_end
        time_profile = np.array([t_str, t_str, t_c, 0, 0])
        return time_profile

    def _cv_dec_profile(self, v_str, v_end, arc_length):
        t_end = np.sqrt((v_str - v_end) / self.j_max)
        L_r4 = (v_str + v_end) * t_end
        t_c = (arc_length - L_r4) / v_str
        time_profile = np.array([0, 0, t_c, t_end, t_end])
        return time_profile

    def _acc_dec_profile(self, v_str, v_end, arc_length):

        def f_type6(t):
            return (
                self.j_max * (v_end - v_str) * t**4
                + 2 * self.j_max * arc_length * t**3
                - (v_str - v_end) ** 2 * t**2
                + 4 * arc_length * v_str * t
                + (v_str + v_end) ** 2 * (v_str - v_end) / self.j_max
                - arc_length**2
            )

        def f_prime_type6(t):
            return (
                4 * self.j_max * (v_end - v_str) * t**3
                + 6 * self.j_max * arc_length * t**2
                - 2 * (v_str - v_end) ** 2 * t
                + 4 * arc_length * v_str
            )

        t_guess = arc_length / (v_str + v_end)
        t_str = newton(f_type6, t_guess, f_prime_type6)
        delta = t_str**2 + (v_str - v_end) / self.j_max
        delta = np.maximum(delta, 0)
        t_end = np.sqrt(delta)
        time_profile = np.array([t_str, t_str, 0, t_end, t_end])
        return time_profile

    def _acc_cv_dec_profile(self, v_str, v_end, arc_length):

        t_str = np.sqrt((self.v_max - v_str) / self.j_max)
        t_end = np.sqrt((self.v_max - v_end) / self.j_max)
        L_acc = (v_str + self.v_max) * t_str
        L_dec = (self.v_max + v_end) * t_end
        t_c = (arc_length - L_acc - L_dec) / self.v_max
        time_profile = np.array([t_str, t_str, t_c, t_end, t_end])
        return time_profile

    def run_phase_1(self, t):
        J, v_str = self.j_max, self.v_str
        a = J * t
        v = v_str + J * t**2 / 2
        s = t * (J * t**2 / 6 + v_str)
        return s, v, a, J

    def run_phase_2(self, t, t_str):
        J, v_str = self.j_max, self.v_str
        a = J * (t_str - t)
        v = v_str + J * t_str**2 / 2 + J * t_str * t - J * t**2 / 2
        s = J * t_str**3 / 6 + J * t_str * t**2 / 2 - J * t**3 / 6 + t_str * v_str + t * (J * t_str**2 + 2 * v_str) / 2
        return s, v, a, -J

    def run_phase_3(self, t, t_str):
        J, v_str = self.j_max, self.v_str
        a = 0
        v = J * t_str**2 + v_str
        s = J * t_str**3 + J * t_str**2 * t + 2 * t_str * v_str + v_str * t
        return s, v, a, 0

    def run_phase_4(self, t, t_str, t_c):
        J, v_str = self.j_max, self.v_str
        a = -J * t
        v = v_str + J * t_str**2 - J * t**2 / 2
        s = (
            J * t_c * t_str**2
            + J * t_str**3
            + J * t_str**2 * t
            - J * t**3 / 6
            + t_c * v_str
            + 2 * t_str * v_str
            + v_str * t
        )
        return s, v, a, -J

    def run_phase_5(self, t, t_str, t_c, t_end):
        J, v_str = self.j_max, self.v_str
        a = J * (-t_end + t)
        v = v_str - J * t_end**2 / 2 - J * t_end * t + J * t_str**2 + J * t**2 / 2
        s = (
            J * t_c * t_str**2
            - J * t_end**3 / 6
            - J * t_end**2 * t / 2
            + J * t_end * t_str**2
            - J * t_end * t**2 / 2
            + J * t_str**3
            + J * t_str**2 * t
            + J * t**3 / 6
            + t_c * v_str
            + t_end * v_str
            + 2 * t_str * v_str
            + v_str * t
        )
        return s, v, a, J


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    t_str = 1
    t_c = 2
    t_end = 1
    p1 = FivePhaseProfile(0, 0, 8, 10, 50, 2, 0.002)
    p1.time_profile = np.array([t_str, t_str, t_c, t_end, t_end])
    total_t = 2 * t_str + t_c + 2 * t_end
    S, V, A, J = [], [], [], []
    T = np.linspace(0, total_t, 100)
    for t in T:
        s, v, a, j = p1.get_motion_state(t)
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

import numpy as np
import pytest

from core.feedrate_profiles.seven_phase_profile import SevenPhaseProfile


def test_example_3_9():
    profile = SevenPhaseProfile(v_str=1, v_end=0, arc_length=10, v_max=5, a_max=10, j_max=30, Ts=0.001)
    profile.generate_profile()
    result = profile.result

    assert np.isclose(result["T_j1"], 0.3333, atol=1e-4)
    assert np.isclose(result["T_j2"], 0.3333, atol=1e-4)
    assert np.isclose(result["T_a"], 0.7333, atol=1e-4)
    assert np.isclose(result["T_v"], 1.1433, atol=1e-4)
    assert np.isclose(result["T_d"], 0.8333, atol=1e-4)
    assert result["v_max_reached"]
    assert result["a_max_reached"]
    assert result["a_min_reached"]


def test_example_3_10():
    profile = SevenPhaseProfile(v_str=1, v_end=0, arc_length=10, v_max=10, a_max=10, j_max=30, Ts=0.001)
    profile.generate_profile()
    result = profile.result

    assert np.isclose(result["T_j1"], 0.3333, atol=1e-4)
    assert np.isclose(result["T_j2"], 0.3333, atol=1e-4)
    assert np.isclose(result["T_a"], 1.0747, atol=1e-4)
    assert np.isclose(result["T_v"], 0, atol=1e-4)
    assert np.isclose(result["T_d"], 1.1747, atol=1e-4)
    assert not result["v_max_reached"]
    assert result["a_max_reached"]
    assert result["a_min_reached"]
    assert np.isclose(result["v_lim"], 8.4136, atol=1e-4)


def test_example_3_11():
    profile = SevenPhaseProfile(v_str=7, v_end=0, arc_length=10, v_max=10, a_max=10, j_max=30, Ts=0.001)
    profile.generate_profile()
    result = profile.result

    assert np.isclose(result["T_j1"], 0.2321, atol=1e-4)
    assert np.isclose(result["T_j2"], 0.2321, atol=1e-4)
    assert np.isclose(result["T_a"], 0.4666, atol=1e-4)
    assert np.isclose(result["T_v"], 0, atol=1e-4)
    assert np.isclose(result["T_d"], 1.4718, atol=1e-4)
    assert not result["v_max_reached"]
    assert np.isclose(result["v_lim"], 8.6329, atol=1e-4)
    assert np.isclose(result["a_lim_a"], 6.9641, atol=1e-4)
    assert np.isclose(result["a_lim_d"], -6.9641, atol=1e-4)


def test_example_3_12():
    profile = SevenPhaseProfile(v_str=7.5, v_end=0, arc_length=10, v_max=10, a_max=10, j_max=30, Ts=0.001)
    profile.generate_profile()
    result = profile.result

    assert np.isclose(result["T_j1"], 0, atol=1e-4)
    assert np.isclose(result["T_j2"], 0.0973, atol=1e-4)
    assert np.isclose(result["T_a"], 0, atol=1e-4)
    assert np.isclose(result["T_v"], 0, atol=1e-4)
    assert np.isclose(result["T_d"], 2.6667, atol=1e-4)
    assert not result["v_max_reached"]
    assert np.isclose(result["v_lim"], 7.5, atol=1e-4)
    assert np.isclose(result["a_lim_a"], 0, atol=1e-4)
    assert np.isclose(result["a_lim_d"], -2.9190, atol=1e-4)


# 添加一个新的测试用例来检查运动状态
def test_motion_state():
    profile = SevenPhaseProfile(v_str=1, v_end=0, arc_length=10, v_max=5, a_max=10, j_max=30, Ts=0.001)
    profile.generate_profile()

    # 检查初始状态
    s, v, a, j = profile.get_motion_state(0)
    assert np.isclose(s, 0)
    assert np.isclose(v, 1)
    assert np.isclose(a, 0)
    assert np.isclose(j, 30)

    # 检查最终状态
    s, v, a, j = profile.get_motion_state(profile.total_time)
    assert np.isclose(s, 10)
    assert np.isclose(v, 0)
    assert np.isclose(a, 0)
    assert np.isclose(j, 30)

    # 检查中间状态（例如，在恒速阶段）
    mid_time = profile.time_profile[0] + profile.time_profile[1] + profile.time_profile[2] + profile.time_profile[3] / 2
    s, v, a, j = profile.get_motion_state(mid_time)
    assert 0 < s < 10
    assert np.isclose(v, 5)
    assert np.isclose(a, 0)
    assert np.isclose(j, 0)

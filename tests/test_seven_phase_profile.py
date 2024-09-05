import numpy as np
import pytest
from seven_phase_profile import double_s_trajectory, SevenPhaseProfile


def test_example_3_9():
    result = double_s_trajectory(h=10, v0=1, v1=0, v_max=5, a_max=10, j_max=30)

    assert np.isclose(result["T_j1"], 0.3333, atol=1e-4)
    assert np.isclose(result["T_j2"], 0.3333, atol=1e-4)
    assert np.isclose(result["T_a"], 0.7333, atol=1e-4)
    assert np.isclose(result["T_v"], 1.1433, atol=1e-4)
    assert np.isclose(result["T_d"], 0.8333, atol=1e-4)
    assert result["v_max_reached"]
    assert result["a_max_reached"]
    assert result["a_min_reached"]


def test_example_3_10():
    result = double_s_trajectory(h=10, v0=1, v1=0, v_max=10, a_max=10, j_max=30)

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
    result = double_s_trajectory(h=10, v0=7, v1=0, v_max=10, a_max=10, j_max=30)

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
    result = double_s_trajectory(h=10, v0=7.5, v1=0, v_max=10, a_max=10, j_max=30)

    assert np.isclose(result["T_j1"], 0, atol=1e-4)
    assert np.isclose(result["T_j2"], 0.0973, atol=1e-4)
    assert np.isclose(result["T_a"], 0, atol=1e-4)
    assert np.isclose(result["T_v"], 0, atol=1e-4)
    assert np.isclose(result["T_d"], 2.6667, atol=1e-4)
    assert not result["v_max_reached"]
    assert np.isclose(result["v_lim"], 7.5, atol=1e-4)
    assert np.isclose(result["a_lim_a"], 0, atol=1e-4)
    assert np.isclose(result["a_lim_d"], -2.9190, atol=1e-4)


def test_seven_phase_profile():
    profile = SevenPhaseProfile(v_str=1, v_end=0, arc_length=10, v_max=5, a_max=10, j_max=30, Ts=0.002)
    profile.generate_profile()

    assert np.isclose(profile.total_time, 2.7099, atol=1e-4)
    
    # Test motion states at different time points
    t_points = [0, 0.3, 0.7, 1.2, 1.8, 2.2, 2.5, profile.total_time]
    expected_states = [
        (0, 1, 0, 30),
        (0.1350, 2.7000, 9.0000, 0),
        (0.4900, 4.2333, 3.0000, -30),
        (2.7500, 5.0000, 0, 0),
        (6.7500, 5.0000, 0, 0),
        (8.7667, 3.7667, -7.0000, 0),
        (9.6650, 1.3000, -9.0000, 30),
        (10, 0, 0, 30)
    ]

    for t, expected in zip(t_points, expected_states):
        s, v, a, j = profile.get_motion_state(t)
        assert np.isclose(s, expected[0], atol=1e-4)
        assert np.isclose(v, expected[1], atol=1e-4)
        assert np.isclose(a, expected[2], atol=1e-4)
        assert np.isclose(j, expected[3], atol=1e-4)

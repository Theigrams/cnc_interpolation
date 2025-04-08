import numpy as np

EPSILON = 1e-16  # Avoid division by zero


def normalize(vectors):
    vectors = np.array(vectors)
    norm = np.linalg.norm(vectors, axis=-1, keepdims=True)
    unit_vectors = vectors / (norm + EPSILON)
    return unit_vectors


def compute_tangents(points):
    path_tangents = np.diff(points, axis=0)
    lengths = np.linalg.norm(path_tangents, axis=1)
    tangents = path_tangents / (lengths.reshape(-1, 1) + EPSILON)
    return tangents, lengths


def compute_normals(vectors):
    normals = np.empty_like(vectors)
    normals[:, 0] = -vectors[:, 1]
    normals[:, 1] = vectors[:, 0]
    normals = normalize(normals)
    return normals


def compute_angles_with_x_axis(vectors):
    x, y = vectors[:, 0], vectors[:, 1]
    angles_with_x_axis = np.arctan2(y, x)
    return angles_with_x_axis


def compute_turning_angles(tangents):
    """
    Calculates the turning angles between consecutive tangent vectors.
    """
    dot_product = np.sum(tangents[:-1] * tangents[1:], axis=1)
    turning_angles = np.arccos(np.clip(dot_product, -1, 1))
    return turning_angles


def compute_bisector(normals):
    bisectors = normals[1:] + normals[:-1]
    bisectors = normalize(bisectors)
    return bisectors


def compute_rotation_matrix(angle):
    return np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])


import numpy as np


def rodrigues_rotate(v, k, theta):
    """
    Rotate a vector `v` around an axis `k` by an angle `theta` using Rodrigues' rotation formula.

    Parameters:
    v (numpy.ndarray): The vector to be rotated.
    k (numpy.ndarray): The axis of rotation.
    theta (float): The angle of rotation in radians.

    Returns:
    numpy.ndarray: The rotated vector.
    """
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    k_dot_v = np.dot(k, v)
    k_cross_v = np.cross(k, v)

    return v * cos_theta + k_cross_v * sin_theta + k * k_dot_v * (1 - cos_theta)


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    t = np.linspace(0, 6 * np.pi, 20)
    x = t * np.cos(t)
    y = t * np.sin(t)
    pts = np.vstack((x, y)).T
    tangents, lengths = compute_tangents(pts)
    normals = compute_normals(tangents)
    bisector = compute_bisector(tangents)
    width = 2
    left_edge = pts[1:-1] + bisector * (width / 2)
    right_edge = pts[1:-1] - bisector * (width / 2)
    plt.plot(pts[:, 0], pts[:, 1], "-r")
    plt.plot(left_edge[:, 0], left_edge[:, 1], "-b")
    plt.plot(right_edge[:, 0], right_edge[:, 1], "-g")
    plt.savefig("gemo.pdf", bbox_inches="tight")
    plt.show()

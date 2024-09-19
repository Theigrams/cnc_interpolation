import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def visualize(g01_points, interpolated_points, deviations=None):
    """
    Visualize the G01 path and the interpolated path.

    Parameters:
    g01_points (np.ndarray): The G01 points.
    interpolated_points (np.ndarray): The interpolated points.
    deviations (np.ndarray, optional): The path deviations. Defaults to None.
    """
    plt.figure(figsize=(12, 8))

    # Plot G01 path
    g01_x, g01_y = g01_points[:, 0], g01_points[:, 1]
    plt.plot(g01_x, g01_y, "b-", label="G01 Path")
    plt.plot(g01_x, g01_y, "bo", markersize=5)

    # Plot interpolated path
    interp_x, interp_y = interpolated_points[:, 0], interpolated_points[:, 1]
    plt.plot(interp_x, interp_y, "r-", label="Interpolated Path", linewidth=1)

    # Plot deviation vectors
    if deviations is not None:
        for i, dev in enumerate(deviations):
            plt.plot([interp_x[i], interp_x[i] + dev[0]], [interp_y[i], interp_y[i] + dev[1]], "g-", alpha=0.5)

    plt.legend()
    plt.title("Comparison of CNC Interpolation and G01 Path")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid(True)
    plt.axis("equal")
    plt.show()


def compute_path_deviation(G01_points: np.ndarray, interpolated_points: np.ndarray) -> np.ndarray:
    """
    Compute the path deviation from each interpolated point to the nearest G01 line segment.

    Parameters:
    G01_points: numpy array of shape (N, D)
        The G01 points defining the line segments.
    interpolated_points: numpy array of shape (M, D)
        The interpolated points.

    Returns:
    deviations: numpy array of shape (M, D)
        The minimal path deviation vectors from each interpolated point to the G01 line segments.
    min_deviations: numpy array of shape (M,)
        The minimal path deviation distance from each interpolated point to the G01 line segments.
    """
    # Get the line segments between consecutive G01 points
    segments_start = G01_points[:-1]  # Starting points of segments (N-1, D)
    segments_end = G01_points[1:]  # Ending points of segments (N-1, D)
    segments = segments_end - segments_start  # Direction vectors of segments (N-1, D)

    # Vectorization setup
    p = interpolated_points[:, np.newaxis, :]  # Shape (M, 1, D)
    a = segments_start[np.newaxis, :, :]  # Shape (1, N-1, D)
    ab = segments[np.newaxis, :, :]  # Shape (1, N-1, D)

    # Compute vectors from segment starts to points
    ap = p - a  # Shape (M, N-1, D)

    # Compute the projection scalar 't' for each point onto each segment
    ab_length_squared = np.einsum("ijk,ijk->ij", ab, ab)  # Squared lengths of segments (1, N-1)
    dot_product = np.einsum("ijk,ijk->ij", ap, ab)  # Dot product of ap and ab (M, N-1)
    t = dot_product / ab_length_squared  # Projection scalar (M, N-1)
    t = np.clip(t, 0, 1)  # Clamp t to the [0,1] interval

    # Compute the closest points on the segments to the interpolated points
    closest_points = a + ab * t[..., np.newaxis]  # Shape (M, N-1, D)

    # Compute the deviations from the interpolated points to the closest points
    deviations = p - closest_points  # Shape (M, N-1, D)
    min_deviations_idx = np.argmin(np.linalg.norm(deviations, axis=2), axis=1)  # Shape (M,)

    # Find the minimum deviation for each interpolated point
    min_deviations = deviations[np.arange(len(interpolated_points)), min_deviations_idx]  # Shape (M, D)

    return min_deviations, np.linalg.norm(min_deviations, axis=1)


def compute_chordal_height_error(G01_points: np.ndarray, interpolated_points: np.ndarray) -> np.ndarray:
    """
    Compute the chordal height error from each G01 point to the nearest interpolated point.

    Parameters:
    G01_points: numpy array of shape (N, D)
        The G01 points defining the original path.
    interpolated_points: numpy array of shape (M, D)
        The interpolated points.

    Returns:
    errors: numpy array of shape (N,)
        The chordal height error from each G01 point to the nearest interpolated point.
    """
    # Compute the distance from each G01 point to all interpolated points
    distances = np.linalg.norm(G01_points[:, np.newaxis, :] - interpolated_points[np.newaxis, :, :], axis=2)

    # Find the minimum distance for each G01 point
    errors = np.min(distances, axis=1)

    return errors


def compute_curvature(points: np.ndarray) -> np.ndarray:
    """
    Calculate the curvature for a given sequence of points.

    Parameters:
    points: numpy array of shape (N, D)
        Sequence of points representing the trajectory.

    Returns:
    curvatures: numpy array of shape (N,)
        Curvature values for each internal point.
    """
    points = np.array(points)
    v1 = points[1:-1] - points[:-2]  # (N-2, D)
    v2 = points[2:] - points[1:-1]  # (N-2, D)

    # Calculate "pseudo" cross product for 2D case
    if points.shape[1] == 2:
        cross = v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]
        cross = np.abs(cross)
    else:
        cross = np.linalg.norm(np.cross(v1, v2), axis=1)

    numerators = 2 * cross
    denominators = (
        np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1) * np.linalg.norm(points[2:] - points[:-2], axis=1)
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        curvatures = np.where(denominators != 0, numerators / denominators, 0)

    # Interpolate curvature for start and end points
    full_curvature = np.zeros(len(points))
    full_curvature[1:-1] = curvatures
    full_curvature[0] = full_curvature[1]
    full_curvature[-1] = full_curvature[-2]

    return full_curvature


def compute_direction_changes(points: np.ndarray) -> np.ndarray:
    """
    Calculate the direction change rate between adjacent segments.

    Parameters:
    points: numpy array, shape (N, D)
        Sequence of points representing the trajectory.

    Returns:
    direction_changes: numpy array, shape (N-2,)
        Direction changes between adjacent segments (in radians).
    """
    # Calculate vectors between adjacent points
    vectors = points[1:] - points[:-1]  # Shape: (N-1, D)

    # Calculate unit vectors
    unit_vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]

    # Calculate dot products of adjacent unit vectors
    dot_products = np.sum(unit_vectors[:-1] * unit_vectors[1:], axis=1)

    # Calculate angle changes using arccos function
    direction_changes = np.arccos(np.clip(dot_products, -1.0, 1.0))

    return direction_changes


def compute_rms(values: np.ndarray) -> float:
    """
    Compute the root mean square (RMS) of given values.

    Parameters:
    values: numpy array
        The numerical sequence to compute the RMS for.

    Returns:
    rms: float
        The computed RMS value.
    """
    return np.sqrt(np.mean(values**2))


if __name__ == "__main__":
    # Load data
    data_name = "rhombic"
    g01_points = np.array([[2, 0], [4, 2], [2, 4], [0, 2], [2, 0]])
    interp_data = np.load(f"experiments/output/{data_name}_interp_data.npy", allow_pickle=True)

    interpolated_points = interp_data["P"]
    T, V, A, J, D = interp_data["T"], interp_data["V"], interp_data["A"], interp_data["J"], interp_data["D"]

    # Compute metrics
    deviation_vectors, distances = compute_path_deviation(g01_points, interpolated_points)
    errors = compute_chordal_height_error(g01_points, interpolated_points)
    curvatures = compute_curvature(interpolated_points)
    delta_curvatures = np.diff(curvatures)
    direction_changes = compute_direction_changes(interpolated_points)
    a_rms = compute_rms(A)
    j_rms = compute_rms(J)

    # Generate report
    metrics = {
        "平均路径偏差": np.mean(distances),
        "最大路径偏差": np.max(distances),
        "平均弦高误差": np.mean(errors),
        "最大弦高误差": np.max(errors),
        "平均曲率": np.mean(curvatures),
        "最大曲率": np.max(curvatures),
        "平均曲率变化率": np.mean(np.abs(delta_curvatures)),
        "最大曲率变化率": np.max(np.abs(delta_curvatures)),
        "平均方向变化率": np.mean(direction_changes),
        "最大方向变化率": np.max(direction_changes),
        "RMS加速度": a_rms,
        "RMS Jerk": j_rms,
        "插补时间": T[-1],
    }

    metrics_df = pd.DataFrame(list(metrics.items()), columns=["指标", "数值"])
    print(metrics_df)

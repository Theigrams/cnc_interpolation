import numpy as np
from scipy.interpolate import BSpline as scipy_BSpline
from scipy.interpolate import interp1d

from utils import adaptive_integral, handle_dimension


class CurveSegment:
    """曲线段基类, 定义了所有曲线段的通用接口"""

    def __init__(self):
        """初始化曲线段, 计算曲线段的长度并进行弧长参数化"""
        self.length = self.get_length(u=1)
        self.arc_length_parameterize()

    def evaluate(self, u) -> np.ndarray:
        """计算曲线段上参数 u 对应的点坐标"""
        raise NotImplementedError

    def derivative(self, u) -> np.ndarray:
        """计算曲线对参数 u 的导数"""
        raise NotImplementedError

    def arc_length_derivative(self, u) -> float:
        """计算弧长导数"""
        derivative = self.derivative(u)
        return np.linalg.norm(derivative, axis=-1)

    def get_length(self, u=1) -> float:
        """计算曲线段从0到u的长度"""
        self._intervals = []
        length = adaptive_integral(self.arc_length_derivative, 0, u, 1e-7, self._intervals)
        return length

    def arc_length_parameterize(self):
        """对曲线进行弧长参数化"""
        # 提取计算弧长时保存的区间和对应弧长值
        parameter_values = [start_u for start_u, end_u, _ in self._intervals]
        arc_lengths = np.cumsum([segment_length for _, _, segment_length in self._intervals])
        parameter_values.append(self._intervals[-1][1])  # 添加最后一个参数值
        arc_lengths[-1] = self.length  # 修正最后一个弧长值为总长度
        arc_lengths = np.insert(arc_lengths, 0, 0)  # 在开始处插入弧长 0
        # 创建弧长到参数的插值函数
        self.arc_parameter = interp1d(arc_lengths, parameter_values, kind="linear")

    def get_u_from_length(self, arc_length):
        """根据弧长返回对应的参数u"""
        target_u = self.arc_parameter(arc_length)
        return target_u

    def get_points(self, n_points, arc_equal=False):
        """返回曲线上的n_points个点坐标, arc_equal为True时, 返回弧长等分的点"""
        if arc_equal:
            arc_lengths = np.linspace(0, self.length, n_points)
            us = self.get_u_from_length(arc_lengths)
        else:
            us = np.linspace(0, 1, n_points)
        return self.evaluate(us)

    def __call__(self, u):
        """调用曲线段时, 返回参数 u 对应的点坐标"""
        return self.evaluate(u)

    def __str__(self):
        return "CurveSegment"

    def __repr__(self):
        return self.__str__()


class Line(CurveSegment):
    def __init__(self, start, end):
        self.start = np.array(start)
        self.end = np.array(end)
        super().__init__()

    @handle_dimension
    def evaluate(self, u):
        if not isinstance(u, (int, float)):
            u = np.array(u).reshape(-1, 1)
        return self.start + (u * (self.end - self.start))

    def derivative(self, u):
        if not isinstance(u, (int, float)):
            u = np.array(u).reshape(-1, 1)
            return np.tile(self.end - self.start, (len(u), 1))
        return self.end - self.start

    def get_length(self, u=1):
        return np.linalg.norm(self.end - self.start) * u

    def arc_length_parameterize(self):
        self.arc_parameter = lambda s: s / self.length

    @handle_dimension
    def curvature(self, u):
        return u * 0

    def __str__(self):
        return f"Line segment: length {self.length}"

    def __repr__(self):
        return f"Line segment: {self.start} -> {self.end}"


class BSpline(CurveSegment):
    def __init__(self, control_points, degree, knots=None):
        self.control_points = np.array(control_points)
        self.degree = degree
        self.build_curve(control_points, degree, knots)
        super().__init__()

    def build_curve(self, control_points, degree, knots=None):
        if knots is not None:
            self.knots = np.array(knots)
            if np.max(self.knots) > 1:
                print("\033[93mWarning: Normalizing knots to [0, 1]\033[0m")
                self.knots /= np.max(self.knots)
        else:
            self.knots = self.generate_knots()
        if len(self.knots) != len(control_points) + degree + 1:
            raise ValueError(f"Knots length ({len(knots)}) must be ({len(control_points) + degree + 1})")
        self.curve = scipy_BSpline(self.knots, control_points, degree)

    def evaluate(self, u):
        return self.curve(u)

    def derivative(self, u):
        return self.curve.derivative(nu=1)(u)

    def curvature(self, u):
        """计算B样条曲线在参数u处的曲率"""
        du = self.curve.derivative(nu=1)(u)
        ddu = self.curve.derivative(nu=2)(u)
        cross_product = np.cross(du, ddu)
        if du.shape[-1] >= 3:
            cross_product = np.linalg.norm(cross_product, axis=-1)
        kappa = cross_product / np.linalg.norm(du, axis=-1) ** 3
        return kappa

    def generate_knots(self):
        """From https://github.com/orbingol/NURBS-Python/blob/5.x/geomdl/knotvector.py"""
        num_repeat, num_ctrlpts = self.degree, len(self.control_points)
        num_segments = num_ctrlpts - (self.degree + 1)
        knot_vector = [0.0 for _ in range(0, num_repeat)]
        knot_vector += np.linspace(0.0, 1.0, num_segments + 2).tolist()
        knot_vector += [1.0 for _ in range(0, num_repeat)]
        return np.array(knot_vector)

    def __str__(self):
        return f"BSpline curve: length {self.length}"

    def __repr__(self):
        str1 = f"BSpline curve: {self.control_points.shape[0]} control points, degree {self.degree}"
        str2 = f"  control points: {self.control_points.tolist()}"
        str3 = f"  knots: {self.knots}"
        str4 = f"  length: {self.length}"
        return f"{str1}\n{str2}\n{str3}\n{str4}"


class Bezier(BSpline):
    def __init__(self, control_points):
        # Bezier曲线的度数 n 等于控制点数减一
        degree = len(control_points) - 1
        # Bezier曲线的节点向量为[0,...,0,1,...,1]，其中0和1各重复n+1次
        knots = np.concatenate(([0] * (degree + 1), [1] * (degree + 1)))
        super().__init__(control_points, degree, knots)

    def __str__(self):
        return f"Bezier curve: length {self.length}"

    def __repr__(self):
        str1 = f"Bezier curve: {len(self.control_points)} control points, degree {self.degree}"
        str2 = f"  control points: {self.control_points.tolist()}"
        str3 = f"  length: {self.length}"
        return f"{str1}\n{str2}\n{str3}"


class NURBS(CurveSegment):
    def __init__(self, control_points, degree, knots, weights=None):
        self.control_points = np.array(control_points)
        self.degree = degree
        self.knots = np.array(knots)
        if weights is None:
            self.weights = np.ones(len(control_points))
        else:
            self.weights = np.array(weights).reshape(-1, 1)
        self.build_curve()
        super().__init__()

    def build_curve(self):
        if np.max(self.knots) > 1:
            print("\033[93mWarning: Normalizing knots to [0, 1]\033[0m")
            self.knots /= np.max(self.knots)
        if len(self.knots) != len(self.control_points) + self.degree + 1:
            raise ValueError(f"Knots length ({len(self.knots)}) must be ({len(self.control_points) + self.degree + 1})")
        # 将控制点和权重合并为齐次坐标
        self.cptsw = self.control_points * self.weights
        control_points_homo = np.concatenate((self.cptsw, self.weights), axis=1)
        self.curve = scipy_BSpline(self.knots, control_points_homo, self.degree)

    @handle_dimension
    def evaluate(self, u):
        point_homo = self.curve(u)
        return point_homo[:, :-1] / point_homo[:, -1][:, None]

    @handle_dimension
    def derivative(self, u):
        C_home, dC_home = self.curve(u), self.curve.derivative(1)(u)
        N, W = C_home[:, :-1], C_home[:, -1][:, None]
        dN, dW = dC_home[:, :-1], dC_home[:, -1][:, None]
        derivative = dN / W - N * dW / W**2
        return derivative

    @handle_dimension
    def second_derivative(self, u):
        C_home, dC_home, ddC_home = self.curve(u), self.curve.derivative(1)(u), self.curve.derivative(2)(u)
        N, W = C_home[:, :-1], C_home[:, -1][:, None]
        dN, dW = dC_home[:, :-1], dC_home[:, -1][:, None]
        ddN, ddW = ddC_home[:, :-1], ddC_home[:, -1][:, None]
        second_derivative = -N * ddW / W**2 + 2 * N * dW**2 / W**3 + ddN / W - 2 * dN * dW / W**2
        return second_derivative

    @handle_dimension
    def curvature(self, u):
        du = self.derivative(u)
        ddu = self.second_derivative(u)
        cross_product = np.cross(du, ddu)
        if du.shape[-1] >= 3:
            cross_product = np.linalg.norm(cross_product, axis=-1)
        kappa = cross_product / np.linalg.norm(du, axis=-1) ** 3
        return kappa

    def __str__(self):
        return f"NURBS curve: length {self.length}"

    def __repr__(self):
        str1 = f"NURBS curve: {self.control_points.shape[0]} control points, degree {self.degree}"
        str2 = f"  control points: {self.control_points.tolist()}"
        str3 = f"  weights: {self.weights[:,0].tolist()}"
        str4 = f"  knots: {self.knots}"
        str5 = f"  length: {self.length}"
        return f"{str1}\n{str2}\n{str3}\n{str4}\n{str5}"


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # B样条曲线
    control_points = np.array([[5.0, 5.0], [10.0, 10.0], [20.0, 15.0], [35.0, 15.0], [45.0, 10.0], [50.0, 5.0]])
    degree = 3
    knotvector = [0.0, 0.0, 0.0, 0.0, 0.33, 0.66, 1.0, 1.0, 1.0, 1.0]
    curve = BSpline(control_points, degree, knotvector)
    print(curve)
    points = curve.get_points(20)
    plt.plot(control_points[:, 0], control_points[:, 1], "ro-", label="Control points")
    plt.plot(points[:, 0], points[:, 1], "b.-", label="B-spline curve")
    plt.legend()
    plt.show()

    # Bezier曲线
    # control_points = np.array([[5.0, 5.0], [20.0, 15.0], [35.0, 15.0], [50.0, 5.0]])
    # curve = Bezier(control_points)
    # print(curve)
    # points = curve.get_points(20)
    # plt.plot(control_points[:, 0], control_points[:, 1], "ro-", label="Control points")
    # plt.plot(points[:, 0], points[:, 1], "b.-", label="Bezier curve")
    # plt.legend()
    # plt.show()

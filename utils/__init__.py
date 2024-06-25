import numpy as np

from utils.geometry import normalize
from utils.integrals import adaptive_integral


def handle_dimension(func):
    """如果u是否为标量,如则将其升维,计算后再将结果降维"""

    def wrapper(*args, **kwargs):
        # 假设u是第一个位置参数
        u = args[1]
        # 检查u是否为标量
        if isinstance(u, (int, float)):
            # 升维
            u = [u]
            # # 更新args中的u
            args = args[:1] + (u,) + args[2:]
            # 调用原始函数
            result = func(*args, **kwargs)
            # 对结果降维
            return result.squeeze()
        else:
            return func(*args, **kwargs)

    return wrapper

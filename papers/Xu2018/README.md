# A circumscribed corner rounding method based on double cubic B-splines for a five-axis linear tool path

## 摘要

这篇文章提出了一种基于双三次B样条的外切圆角方法，用于五轴线性工具路径的圆角处理。传统的内切圆角方法由于几何不连续性，导致进给率波动和过度加速，从而影响加工精度和质量。本文提出的方法通过配置控制点和使用双三次B样条平滑工具轨迹，能够有效减少曲率，提高角落区域的进给率，保证工具方向的连续变化，而无需迭代计算。仿真结果表明，该方法在提高加工速度和效率方面具有显著优势。

## 算法流程

1. **初始配置与输入**
   - 收集五轴NC加工中的线性段数据，包括工具路径的底部轨迹和顶部轨迹。
   - 确定角点和线性段的控制点。

2. **双三次B样条的构建**
   - 定义双三次B样条的基本函数和节点向量。
   - 配置控制点以生成双三次B样条，保证过渡曲线具有所需的对称性和连续性。

3. **G2连续性的实现**
   - 确保过渡曲线与线性段在接合处具有相同的切线单位和曲率向量，以实现G2连续性。
   - 计算最大逼近误差，并根据逼近误差调整控制点的位置。

4. **比例调整策略**
   - 如果过渡长度无法满足长度约束，则采用比例调整策略调整过渡长度。
   - 使用调整后的过渡长度重新配置控制点。

5. **工具方向平滑**
   - 通过双三次B样条分别平滑底部和顶部轨迹。
   - 确保工具方向在过渡曲线的入口和出口处连续变化。
   - 计算工具方向的几何导数，调整控制点以保证第一和第二几何导数的连续性。

6. **仿真与验证**
   - 进行仿真测试，验证所提方法在提高加工速度和效率方面的有效性。
   - 比较外切圆角方法与内切圆角方法在进给率和加速度方面的性能。

## Reference

```bibtex
@article{xu2018CircumscribedCornerRounding,
  title = {A Circumscribed Corner Rounding Method Based on Double Cubic {{B-splines}} for a Five-Axis Linear Tool Path},
  author = {Xu, Fuyang and Sun, Yuwen},
  year = {2018},
  month = jan,
  journal = {The International Journal of Advanced Manufacturing Technology},
  volume = {94},
  number = {1-4},
  pages = {451--462},
  issn = {0268-3768, 1433-3015},
  doi = {10.1007/s00170-017-0869-x},
  url = {http://link.springer.com/10.1007/s00170-017-0869-x}
}
```

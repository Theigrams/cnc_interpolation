# A real-time look-ahead interpolation methodology with curvature-continuous B-spline transition scheme for CNC machining of short line segments

## 算法流程

### 1. 路径平滑

- 采用 B 样条曲线对拐角处的路径进行平滑过渡，确保路径的曲率连续性
- B 样条和直线段前后连接，保证路径的光滑性
- 使用adaptive subdivision scheme来保证 B 样条不会过长

### 2. 进给率规划

- 采用五阶段双S曲线来规划进给率，确保进给率的平滑性

### 3. 前瞻规划

- 根据 chord error, curvature, normal acceleration and jerk 来计算最大进给率
- 使用双向扫描算法来计算连接点处的进给率限制

### 4. 路径插补

- 根据进给率调度结果，计算每个采样时刻的刀具位置，从而实现连续的加工路径。

## Note

- Lin2007 中 Table 3 中 Type(III) 的公式有误，应该为

  $$
  f\left(T_{\text{end}}\right)=-J_{\text{max}} T_{\text{end}}^3+2 V_{\text{str}} T_{\text{end}}-L_{\text {seg }}^m=0
  $$

## Reference

```bibtex
@article{zhao2013RealtimeLookaheadInterpolation,
  title = {A Real-Time Look-Ahead Interpolation Methodology with Curvature-Continuous {{B-spline}} Transition Scheme for {{CNC}} Machining of Short Line Segments},
  author = {Zhao, Huan and Zhu, LiMin and Ding, Han},
  year = {2013},
  month = feb,
  journal = {International Journal of Machine Tools and Manufacture},
  volume = {65},
  pages = {88--98},
  issn = {08906955},
  doi = {10.1016/j.ijmachtools.2012.10.005},
  url = {https://linkinghub.elsevier.com/retrieve/pii/S0890695512001885}
}

@article{lin2007DevelopmentDynamicsbasedNURBS,
  title = {Development of a Dynamics-Based {{NURBS}} Interpolator with Real-Time Look-Ahead Algorithm},
  author = {Lin, Ming-Tzong and Tsai, Meng-Shiun and Yau, Hong-Tzong},
  year = {2007},
  month = dec,
  journal = {International Journal of Machine Tools and Manufacture},
  volume = {47},
  number = {15},
  pages = {2246--2262},
  issn = {08906955},
  doi = {10.1016/j.ijmachtools.2007.06.005},
  url = {https://linkinghub.elsevier.com/retrieve/pii/S0890695507001174},
  copyright = {https://www.elsevier.com/tdm/userlicense/1.0/}
}
```

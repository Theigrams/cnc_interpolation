# A real-time look-ahead interpolation methodology with curvature-continuous B-spline transition scheme for CNC machining of short line segments

## 算法流程

### 1. 路径平滑

- 采用 B 样条曲线对拐角处的路径进行平滑过渡，确保路径的曲率连续性
- B 样条和直线段前后连接，保证路径的光滑性
- 使用adaptive subdivision scheme来保证 B 样条不会过长

### 2. 进给率规划

- 采用五阶段双S曲线来规划进给率，确保进给率的平滑性
- 但是可能出现加速度超出限制的情况

### 3. 前瞻规划

- 根据 chord error, curvature, normal acceleration and jerk 来计算最大进给率
- 使用双向扫描算法来计算连接点处的进给率限制

### 4. 路径插补

- 根据进给率调度结果，计算每个采样时刻的刀具位置，从而实现连续的加工路径。

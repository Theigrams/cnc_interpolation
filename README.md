# CNC 插补算法项目

这是一个使用 Python 实现的 CNC (计算机数控) 插补算法项目。项目的主要功能包括曲线插补、进给率规划、前瞻算法等，旨在提升数控机床的插补精度与平滑性。


## 项目介绍

本项目实现了CNC插补相关的各种算法，适用于高精度数控系统。主要的应用包括直线与曲线插补、基于运动学的进给率规划、前瞻算法以及刀具路径的平滑处理。

## 主要功能

- **曲线表示与插补**: 支持插补多种参数曲线，如直线、B样条曲线、Bezier曲线、NURBS曲线等。
- **进给率规划**: 实现了基于五阶段和七阶段模型的运动学规划。
- **前瞻算法**: 使用双向扫描实现速度限制计算，优化插补过程。
- **刀具路径平滑**: 实现路径转折处的平滑过渡，减少机械抖动。
- **插补器**: 基于运动学规划生成离散的插补点序列。

## 项目结构

```
cnc_interpolation/
├── core/                       # 核心算法模块
│   ├── curve.py                # 各种参数曲线的定义与计算
│   ├── feedrate_profiles       # 进给率曲线的实现
│   ├── look_ahead.py           # 前瞻算法实现
│   ├── toolpath.py             # 刀具路径的定义与处理
│   └── interpolator.py         # 插补器实现
├── docs/                       # 原理文档与算法说明
├── examples/                   # 示例代码和Jupyter笔记本
├── experiments/                # 实验与性能测试代码
├── papers/                     # 第三方论文的复现实现
├── tests/                      # 项目测试模块
└── utils/                      # 工具函数模块
```

## 安装说明

1. 克隆项目到本地：

   ```bash
   git clone https://github.com/Theigrams/cnc_interpolation.git
   cd cnc_interpolation
   ```

2. 安装依赖：

   ```bash
   pip install -r requirements.txt
   ```

3. 添加项目路径到 Python 环境变量：

   ```bash
   export PYTHONPATH=$PYTHONPATH:$(pwd)
   ```

## 使用说明

1. **示例程序**: 查看 `examples/rhombic_interpolation.ipynb` 了解基础的插补算法。

2. **运行测试**: 使用 `pytest` 运行单元测试：

   ```bash
   pytest
   ```

3. **修改进给率规划算法**: 可以通过修改 `core/feedrate_scheduler.py` 来切换不同的进给率规划模型：
   
   - 使用七阶段模型：

     ```python
     from core.feedrate_profiles import SevenPhaseProfile as FeedrateProfile
     ```

   - 使用五阶段模型：

     ```python
     from core.feedrate_profiles import FivePhaseProfile as FeedrateProfile
     ```

## 主要模块

- **curve.py**: 定义了支持的各种参数曲线，包括 B样条、Bezier 曲线和 NURBS。
- **feedrate_profiles**: 实现了不同阶段的进给率曲线模型，包括五阶段和七阶段模型。
- **look_ahead.py**: 前瞻算法，优化路径规划中的速度限制。
- **toolpath.py**: 定义和处理刀具路径。
- **interpolator.py**: 负责生成离散的插补点序列。

## 参考论文

- **Xu2018**: 基于双三次 B 样条的外切圆角方法。
- **Zhao2013**: 基于 B 样条过渡的实时前瞻插补方法。
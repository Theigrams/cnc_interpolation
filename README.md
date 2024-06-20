# CNC 插补算法项目

## 项目结构

```bash
cnc_interpolation/
├── core/                      # 核心算法模块
│   ├── __init__.py
│   ├── curve.py              # 曲线类定义 (直线、B样条等)
│   ├── feedrate_planner.py   # 进给率规划类
│   └── interpolator.py       # 插补器类
├── docs/                      # 原理文档
├── examples/                  # 示例代码
├── utils/                     # 工具函数模块
├── tests/                     # 测试模块
├── papers/                    # 第三方论文复现
└── README.md                  # 项目说明
```

## 使用说明

## 添加路径

```bash
export export PYTHONPATH=$PYTHONPATH:$(pwd)
```

### 测试

```bash
pytest
```

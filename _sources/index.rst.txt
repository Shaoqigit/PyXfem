.. PyAcoustiX documentation master file

=============================================
PyAcoustiX / SAcouS 文档
=============================================

**PyAcoustiX** (Simple Acoustic Simulator - SAcouS) 是一个基于纯Python实现的有限元方法(FEM)声学仿真库。

该库主要用于声学/振动声学仿真，但也可扩展用于其他物理场问题。它旨在为研究、教育或工业原型开发提供最简单、清晰且易于理解的FEM工具。

主要特性
========

- **高阶形函数**: 除经典的线性/二次Lagrange多项式外，还支持高阶Lobatto形函数 (最高p=4)
- **变阶单元**: 每个单元可以独立设置插值阶数
- **扩展单元技术**: 使用双节点技术处理不连续性问题
- **丰富的材料模型**: 支持多种声学材料，包括刚性多孔材料、Limp模型和Biot-UP模型
- **FEM-TMM耦合**: 支持传递矩阵法(TMM)与FEM耦合，用于复杂几何中的多层结构建模
- **模态降阶**: 支持模态域降阶方法(MoR)加速计算

快速开始
========

安装
----

.. code-block:: bash

   pip install pyacoustix

基本使用
--------

.. code-block:: python

   from SAcouS import *
   
   # 定义材料
   air = Air('air')
   
   # 创建网格
   mesh = Mesh1D(nodes, connectivity)
   
   # 求解...

命令行工具
----------

.. code-block:: bash

   # 使用输入文件运行仿真
   sacous -i input_file.axi
   
   # 进入交互式控制台
   sacous

文档目录
========

.. toctree::
   :maxdepth: 2
   :caption: 用户指南

   user_guide/installation
   user_guide/quickstart
   user_guide/tutorial_1d
   user_guide/tutorial_2d
   user_guide/tutorial_3d
   user_guide/cli_usage

.. toctree::
   :maxdepth: 2
   :caption: API参考

   api_reference/index

.. toctree::
   :maxdepth: 2
   :caption: 理论基础

   theory/fem_basics
   theory/helmholtz
   theory/jca_model
   theory/biot_theory
   theory/tmm_theory

.. toctree::
   :maxdepth: 2
   :caption: 示例

   examples/index

.. toctree::
   :maxdepth: 2
   :caption: 开发者指南

   developer_guide/architecture
   developer_guide/contributing
   developer_guide/testing

索引
====

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

项目链接
========

- `GitHub仓库 <https://github.com/Shaoqigit/PyXfem>`_
- `PyPI页面 <https://pypi.org/project/PyAcoustiX/>`_
- `Docker镜像 <https://hub.docker.com/r/shaoqiwu/pyacoustix>`_

许可证
======

本项目采用 MIT 许可证。详见 `LICENSE <https://github.com/Shaoqigit/PyXfem/blob/main/LICENSE>`_ 文件。

作者
====

- **Shaoqi WU** - shaoqiwu@outlook.com

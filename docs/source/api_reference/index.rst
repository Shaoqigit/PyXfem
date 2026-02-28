API参考
=======

本部分提供PyAcoustiX所有模块的完整API参考。

.. toctree::
   :maxdepth: 2

   mesh
   materials
   postprocess
   acxfem
   acxtmm
   acxmor
   interface

快速索引
--------

核心模块
~~~~~~~~

- :mod:`SAcouS.Mesh`: 网格处理 (1D/2D/3D)
- :mod:`SAcouS.Materials`: 材料模型定义
- :mod:`SAcouS.PostProcess`: 后处理和可视化

FEM模块
~~~~~~~

- :mod:`SAcouS.acxfem.Basis`: 基函数实现
- :mod:`SAcouS.acxfem.Assembly`: 矩阵组装
- :mod:`SAcouS.acxfem.Solver`: 线性求解器
- :mod:`SAcouS.acxfem.DofHandler`: 自由度管理
- :mod:`SAcouS.acxfem.BCsImpose`: 边界条件
- :mod:`SAcouS.acxfem.Quadratures`: 数值积分

高级功能
~~~~~~~~

- :mod:`SAcouS.acxtmm`: 传递矩阵法 (TMM)
- :mod:`SAcouS.acxmor`: 模态降阶 (MoR)
- :mod:`SAcouS.interface`: 输入/输出接口

使用示例
--------

导入模块
~~~~~~~~

.. code-block:: python

   # 导入整个包
   from SAcouS import *
   
   # 或导入特定模块
   from SAcouS.Mesh import Mesh1D, Mesh2D
   from SAcouS.Materials import Air, EquivalentFluid
   from SAcouS.acxfem import HelmholtzAssembler

类型注解
~~~~~~~~

PyAcoustiX使用类型注解提高代码可读性:

.. code-block:: python

   def solve(mesh: Mesh2D, material: BaseMaterial, freq: float) -> np.ndarray:
       ...

常见类型:

- ``np.ndarray``: NumPy数组
- ``complex128``: 复数类型
- ``Dict[str, Any]``: 字典
- ``Optional[T]``: 可选类型

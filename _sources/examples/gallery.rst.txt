示例画廊
========

以下是PyAcoustiX的一些典型应用示例。

1D Helmholtz问题
----------------

**问题描述**: 声波在Kundt管中的传播

**材料**: 空气 + JCA多孔材料

**结果**: 压力分布

.. image:: ../_static/example_1d.png
   :width: 600px
   :align: center

.. code-block:: python

   # 关键代码
   from SAcouS.Materials import Air, EquivalentFluid
   from SAcouS.Mesh import Mesh1D
   
   air = Air('air')
   porous = EquivalentFluid('foam', 0.98, 3750, 1.17, 742e-6, 110e-6)
   
   # 创建网格和求解...

2D腔体问题
----------

**问题描述**: 2D矩形腔体的声学响应

**边界条件**: 左侧速度激励，其他刚性壁

**结果**: 声压级分布

.. image:: ../_static/example_2d.png
   :width: 600px
   :align: center

多层吸声结构
------------

**问题描述**: 多层材料的吸声系数

**方法**: TMM + FEM耦合

**结果**: 频率相关的吸声系数

.. image:: ../_static/example_absorption.png
   :width: 600px
   :align: center

Biot多孔弹性材料
----------------

**问题描述**: 饱和土壤中的波传播

**材料**: Biot多孔弹性材料

**结果**: 快波和慢波传播

频率响应分析
------------

**问题描述**: 扬声器频率响应

**方法**: 频率扫描 + MoR加速

**结果**: SPL vs 频率

.. image:: ../_static/example_frf.png
   :width: 600px
   :align: center

3D复杂几何
----------

**问题描述**: 汽车内部声场

**网格**: Gmsh生成的四面体网格

**结果**: 3D声压分布

更多示例
--------

查看 ``tests/`` 目录获取更多示例代码。

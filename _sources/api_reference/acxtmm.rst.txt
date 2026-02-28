acxtmm模块
==========

acxtmm模块实现了传递矩阵法(Transfer Matrix Method, TMM)，用于多层声学结构的快速计算。

.. module:: SAcouS.acxtmm

TMM核心
-------

.. autoclass:: SAcouS.acxtmm.Tmm.TMMFluid
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: SAcouS.acxtmm.Tmm.TMMElastic
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: SAcouS.acxtmm.Tmm.TMMPoroElastic1
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: SAcouS.acxtmm.Tmm.TMMPoroElastic2
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: SAcouS.acxtmm.Tmm.TMMPoroElastic3
   :members:
   :undoc-members:
   :show-inheritance:

导纳组装器
----------

.. autoclass:: SAcouS.acxtmm.AdmAssembler.AdmAssembler
   :members:
   :undoc-members:
   :show-inheritance:

导纳基函数
----------

.. autoclass:: SAcouS.acxtmm.AdmBasis.AdmFluid
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: SAcouS.acxtmm.AdmBasis.AdmElastic
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: SAcouS.acxtmm.AdmBasis.AdmPoroElastic
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: SAcouS.acxtmm.AdmBasis.AdmPoroElastic2
   :members:
   :undoc-members:
   :show-inheritance:

使用示例
--------

多层结构分析
~~~~~~~~~~~~

.. code-block:: python

   from SAcouS.acxtmm import TMMFluid
   from SAcouS.Materials import Air, EquivalentFluid
   
   # 定义材料层
   air = Air('air')
   foam = EquivalentFluid('foam', 0.98, 3750, 1.17, 742e-6, 110e-6)
   
   # 层厚度 [m]
   thicknesses = [0.01, 0.05, 0.01]  # 空气-泡沫-空气
   materials = [air, foam, air]
   
   # 创建TMM模型
   freq = 1000
   tmm = TMMFluid(materials, thicknesses, freq)
   
   # 计算吸声系数
   absorption = tmm.absorption_coefficient()
   print(f"吸声系数: {absorption:.3f}")
   
   # 计算表面阻抗
   Zs = tmm.surface_impedance()
   print(f"表面阻抗: {Zs:.2f}")

FEM-TMM耦合
~~~~~~~~~~~

.. code-block:: python

   from SAcouS.acxtmm import AdmAssembler
   
   # 在FEM模型中使用TMM作为边界条件
   assembler = AdmAssembler(
       dof_handler, mesh.subdomains, tmm_model
   )

频率扫描
~~~~~~~~

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   
   freqs = np.linspace(100, 5000, 100)
   absorption = []
   
   for freq in freqs:
       tmm = TMMFluid(materials, thicknesses, freq)
       absorption.append(tmm.absorption_coefficient())
   
   plt.plot(freqs, absorption)
   plt.xlabel('Frequency (Hz)')
   plt.ylabel('Absorption Coefficient')
   plt.show()

理论背景
--------

传递矩阵法基于平面波假设，每层材料用2×2传递矩阵表示:

.. math::

   \begin{bmatrix} p \\ v \end{bmatrix}_{out} = 
   \begin{bmatrix} T_{11} & T_{12} \\ T_{21} & T_{22} \end{bmatrix}
   \begin{bmatrix} p \\ v \end{bmatrix}_{in}

对于多层结构，总传递矩阵是各层矩阵的乘积:

.. math::

   T_{total} = T_n \cdot T_{n-1} \cdots T_2 \cdot T_1

从总传递矩阵可以计算:

- 表面阻抗: :math:`Z_s = \frac{T_{11}}{T_{21}}`
- 吸声系数: :math:`\alpha = 1 - |R|^2`，其中 :math:`R = \frac{Z_s - Z_0}{Z_s + Z_0}`

注意事项
--------

- TMM假设平面波传播，适用于垂直入射
- 对于斜入射，需要修正传递矩阵
- 对于复杂几何，建议使用FEM-TMM耦合

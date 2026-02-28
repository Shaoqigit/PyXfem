acxmor模块
==========

acxmor模块实现了模态降阶(Modal Order Reduction, MoR)方法，用于加速大规模声学问题的求解。

.. module:: SAcouS.acxmor

模态降阶
--------

.. autoclass:: SAcouS.acxmor.ModalReduction.EigenSolver
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: SAcouS.acxmor.ModalReduction.ModalReduction
   :members:
   :undoc-members:
   :show-inheritance:

降阶基求解器
------------

.. autoclass:: SAcouS.acxmor.myRBSolver.RBSolver_fromResidual
   :members:
   :undoc-members:
   :show-inheritance:

经验插值方法
------------

.. autoclass:: SAcouS.acxmor.myEIM.nonIntrusiveEIMV2
   :members:
   :undoc-members:
   :show-inheritance:

使用示例
--------

基本模态降阶
~~~~~~~~~~~~

.. code-block:: python

   from SAcouS.acxmor import ModalReduction
   from SAcouS.acxfem import HelmholtzAssembler
   
   # 创建原始FEM模型
   assembler = HelmholtzAssembler(dof_handler, mesh.subdomains)
   
   # 创建模态降阶模型
   mor = ModalReduction(assembler)
   
   # 计算模态基 (前20阶)
   mor.compute_modes(n_modes=20)
   
   # 在特定频率求解
   omega = 2 * np.pi * 1000
   sol_reduced = mor.solve(omega)
   
   # 重构完整解
   sol_full = mor.reconstruct(sol_reduced)

频率扫描加速
~~~~~~~~~~~~

.. code-block:: python

   import numpy as np
   
   # 频率范围
   freqs = np.linspace(100, 5000, 500)
   
   # 使用MoR加速
   results = []
   for freq in freqs:
       omega = 2 * np.pi * freq
       sol_r = mor.solve(omega)
       sol_f = mor.reconstruct(sol_r)
       results.append(sol_f)

与直接求解对比
~~~~~~~~~~~~~~

.. code-block:: python

   import time
   
   # 直接求解
   t0 = time.time()
   for freq in freqs:
       omega = 2 * np.pi * freq
       solver = LinearSolver(dof_handler=dof_handler)
       solver.solve(K, F)
   t_direct = time.time() - t0
   
   # MoR求解
   t0 = time.time()
   mor.compute_modes(n_modes=20)
   for freq in freqs:
       omega = 2 * np.pi * freq
       mor.solve(omega)
   t_mor = time.time() - t0
   
   print(f"直接求解: {t_direct:.2f}s")
   print(f"MoR求解: {t_mor:.2f}s")
   print(f"加速比: {t_direct/t_mor:.1f}x")

理论背景
--------

模态降阶基于以下思想:

1. **模态展开**: 解可以表示为模态基的线性组合

   .. math::

      u(\omega) \approx \sum_{i=1}^{N} \alpha_i(\omega) \phi_i

2. **降阶系统**: 将大规模系统投影到低维模态空间

   .. math::

      \hat{K}(\omega) \alpha = \hat{F}

   其中 :math:`\hat{K} = \Phi^T K \Phi`，:math:`\Phi = [\phi_1, \phi_2, ..., \phi_N]`

3. **快速求解**: 降阶系统规模小，可快速求解

模态基的选择
~~~~~~~~~~~~

- **特征模态**: 系统的固有振型
- **Krylov子空间**: 适用于宽频带问题
- **POD模态**: 基于快照的最优基

适用场景
--------

模态降阶特别适用于:

- 大规模频率扫描问题
- 实时仿真需求
- 参数化研究
- 优化设计

限制
----

- 需要预计算模态基
- 对于强非线性问题效果有限
- 模态数量需要足够以捕捉物理现象

教程1: 1D Helmholtz问题
=======================

本教程将详细介绍如何使用PyAcoustiX求解一维Helmholtz方程。我们将模拟声波在Kundt管中的传播，这是声学中的经典问题。

问题描述
--------

考虑一个长度为 :math:`L` 的管子，其中包含两种材料:

- 前半段 (:math:`0 \leq x \u003c L/2`): 空气
- 后半段 (:math:`L/2 \leq x \leq L`): 多孔吸声材料

在左端 (:math:`x=0`) 施加速度边界条件，右端 (:math:`x=L`) 为自由出射边界。

控制方程
~~~~~~~~

一维Helmholtz方程:

.. math::

   \frac{d^2 p}{dx^2} + k^2 p = 0

其中 :math:`k = \omega/c` 是波数，:math:`\omega = 2\pi f` 是角频率，:math:`c` 是声速。

完整代码
--------

.. code-block:: python

   import sys
   import numpy as np
   import matplotlib.pyplot as plt
   
   from SAcouS.acxfem import Lobbato1DElement
   from SAcouS.Materials import Air, EquivalentFluid
   from SAcouS.PostProcess import PostProcessField
   from SAcouS.Mesh import Mesh1D
   from SAcouS.acxfem import GeneralDofHandler1D, FESpace
   from SAcouS.acxfem import HelmholtzAssembler
   from SAcouS.acxfem import ApplyBoundaryConditions
   from SAcouS.acxfem import check_material_compability, plot_matrix_partten
   from SAcouS.acxfem import LinearSolver
   from analytical.fluid_sol import DoubleleLayerKundltTube

   def tutor_case_1():
       # ====================== 物理问题设置 ======================
       print("设置物理问题...")
       
       # 定义材料
       air = Air('classical air')
       
       # JCA多孔材料参数
       phi = 0.98           # 孔隙率
       sigma = 3.75e3       # 流阻 [Pa·s/m²]
       alpha = 1.17         # 曲折度
       Lambda_prime = 742e-6  # 热特征长度 [m]
       Lambda = 110e-6      # 粘滞特征长度 [m]
       
       xfm = EquivalentFluid('xfm', phi, sigma, alpha, 
                             Lambda_prime, Lambda)
       
       # 频率设置
       freq = 2000          # 频率 [Hz]
       omega = 2 * np.pi * freq
       
       # ====================== 网格定义 ======================
       print("创建网格...")
       num_elem = 100       # 单元数
       num_nodes = num_elem + 1
       nodes = np.linspace(-1, 1, num_nodes)
       
       elem_connec1 = np.arange(0, num_elem)
       elem_connec2 = np.arange(1, num_nodes)
       connectivity = np.vstack((elem_connec1, elem_connec2)).T
       
       mesh = Mesh1D(nodes, connectivity)
       elements2node = mesh.get_mesh_coordinates()
       
       # 设置子域
       air_elements = np.arange(0, int(num_elem / 2))
       xfm_elements = np.arange(int(num_elem / 2), num_elem)
       mesh.set_subdomains({air: air_elements, xfm: xfm_elements})
       
       # 检查材料兼容性
       check_material_compability(mesh.subdomains)
       
       print(f"创建了 {num_elem} 个单元的网格")
       
       # ====================== FEM基函数 ======================
       print("设置FEM基函数...")
       order = 2            # 全局阶数
       
       Pf_bases = []
       for mat, elems in mesh.subdomains.items():
           if mat.TYPE == 'Fluid':
               Pf_bases += [
                   Lobbato1DElement('Pf', order, elements2node[elem]) 
                   for elem in elems
               ]
       
       # ====================== 自由度处理 ======================
       Helmholtz_dof_handler = GeneralDofHandler1D(['Pf'], Pf_bases)
       print(f"全局自由度数量: {Helmholtz_dof_handler.nb_dofs}")
       
       # ====================== 矩阵组装 ======================
       Helmholtz_assember = HelmholtzAssembler(
           Helmholtz_dof_handler, mesh.subdomains, dtype=np.complex128
       )
       Helmholtz_assember.assembly_global_matrix(Pf_bases, 'Pf', omega)
       left_hand_matrix = Helmholtz_assember.get_global_matrix()
       
       # 显示矩阵稀疏模式
       plot_matrix_partten(left_hand_matrix)
       
       # ====================== 边界条件 ======================
       print("施加边界条件...")
       fe_space = FESpace(mesh, Pf_bases)
       right_hand_vec = np.zeros(Helmholtz_assember.nb_global_dofs, 
                                 dtype=np.complex128)
       
       # 自然边界条件 (速度边界)
       nature_bcs = {
           'type': 'fluid_velocity',
           'value': 1 * np.exp(-1j * omega),
           'position': -1
       }
       
       BCs_applier = ApplyBoundaryConditions(mesh, fe_space, 
                                             left_hand_matrix, 
                                             right_hand_vec, omega)
       BCs_applier.apply_nature_bc(nature_bcs, var='Pf')
       
       # ====================== 求解 ======================
       print("求解线性系统...")
       linear_solver = LinearSolver(dof_handler=Helmholtz_dof_handler)
       linear_solver.solve(left_hand_matrix, right_hand_vec)
       sol = linear_solver.u
       
       # ====================== 解析解对比 ======================
       print("计算解析解...")
       kundlt_tube = DoubleleLayerKundltTube(1, 1, air, xfm, omega, nature_bcs)
       ana_sol = kundlt_tube.sol_on_mesh(mesh, sol_type='pressure')
       
       # ====================== 后处理 ======================
       post_processer = PostProcessField(mesh.nodes, r'1D Helmholtz (2000Hz)')
       post_processer.plot_sol(
           (np.real(sol), 'FEM (p=2)', 'solid'),
           (np.real(ana_sol), 'Analytical', 'dashed')
       )
       
       # 计算误差
       error = post_processer.compute_error(sol, ana_sol)
       print(f"相对L2误差: {error:.2%}")
       
       plt.show()
       
       return sol

   if __name__ == "__main__":
       result = tutor_case_1()

详细解释
--------

材料定义
~~~~~~~~

**空气材料**

.. code-block:: python

   air = Air('classical air')

``Air`` 类预定义了标准大气条件下的空气属性:

- 密度 :math:`\rho = 1.213` kg/m³
- 声速 :math:`c = 343` m/s
- 温度 :math:`T = 293.15` K
- 压力 :math:`P = 101325` Pa

**JCA等效流体**

.. code-block:: python

   xfm = EquivalentFluid('xfm', phi, sigma, alpha, 
                         Lambda_prime, Lambda)

Johnson-Champoux-Allard (JCA) 模型是描述多孔吸声材料的标准模型。参数含义:

- **phi** (:math:`\phi`): 孔隙率，材料中孔隙体积与总体积之比
- **sigma** (:math:`\sigma`): 流阻，描述流体流经材料时的阻力
- **alpha** (:math:`\alpha`): 曲折度，描述孔隙路径的复杂程度
- **Lambda_prime** (:math:`\Lambda'`): 热特征长度，影响热损耗
- **Lambda** (:math:`\Lambda`): 粘滞特征长度，影响粘滞损耗

网格创建
~~~~~~~~

.. code-block:: python

   nodes = np.linspace(-1, 1, num_nodes)
   connectivity = np.vstack((elem_connec1, elem_connec2)).T
   mesh = Mesh1D(nodes, connectivity)

对于1D问题:

- ``nodes``: 节点坐标数组，形状为 ``(n_nodes,)``
- ``connectivity``: 单元连接关系，形状为 ``(n_elements, 2)``

每个单元由两个节点定义，节点编号从0开始。

子域设置
~~~~~~~~

.. code-block:: python

   mesh.set_subdomains({air: air_elements, xfm: xfm_elements})

子域用于定义不同材料占据的区域。``air_elements`` 和 ``xfm_elements`` 是NumPy数组，包含单元编号。

基函数
~~~~~~

.. code-block:: python

   Pf_bases = [Lobbato1DElement('Pf', order, elements2node[elem]) 
               for elem in elems]

``Lobbato1DElement`` 是高阶谱元法使用的基函数:

- **name**: 变量名 ('Pf' 表示声压)
- **order**: 多项式阶数 (1=线性, 2=二次, 3=三次, 4=四次)
- **nodes**: 单元节点坐标

Lobatto多项式的特点是节点包含端点，有利于施加边界条件。

边界条件
~~~~~~~~

**速度边界条件**

.. code-block:: python

   nature_bcs = {
       'type': 'fluid_velocity',
       'value': 1 * np.exp(-1j * omega),
       'position': -1
   }

- **type**: 边界条件类型
- **value**: 边界值 (复数)
- **position**: 边界位置 (-1表示右端)

速度边界条件对应Neumann边界条件:

.. math::

   \frac{\partial p}{\partial n} = -j\omega\rho_0 v_n

结果分析
~~~~~~~~

代码计算了FEM解与解析解的相对L2误差:

.. code-block:: python

   error = post_processer.compute_error(sol, ana_sol)

通常，对于足够细的网格和适当的阶数，误差应小于1%。

练习
----

1. **改变频率**: 修改 ``freq`` 变量，观察解的变化
2. **改变网格密度**: 修改 ``num_elem``，观察收敛性
3. **改变多项式阶数**: 修改 ``order``，比较高阶和低阶的精度
4. **改变材料参数**: 调整JCA参数，观察对吸声的影响

下一步
------

- :doc:`tutorial_2d`: 学习2D声学问题
- :doc:`../theory/helmholtz`: 深入了解Helmholtz方程理论
- :doc:`../theory/jca_model`: 学习JCA模型理论

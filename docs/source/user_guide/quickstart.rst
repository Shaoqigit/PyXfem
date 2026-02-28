快速开始
========

本指南将帮助你快速上手PyAcoustiX，通过一个简单的1D声学问题示例。

示例：1D Helmholtz问题
----------------------

我们将求解一个一维Helmholtz方程，模拟声波在两种材料(空气和多孔材料)中的传播。

完整代码
~~~~~~~~

.. code-block:: python

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

   # 1. 定义物理问题
   freq = 2000  # 频率 (Hz)
   omega = 2 * np.pi * freq  # 角频率
   
   # 定义材料
   air = Air('air')
   
   # JCA多孔材料参数
   phi = 0.98           # 孔隙率
   sigma = 3.75e3       # 流阻
   alpha = 1.17         # 曲折度
   Lambda_prime = 742e-6  # 热特征长度
   Lambda = 110e-6      # 粘滞特征长度
   
   porous = EquivalentFluid('porous', phi, sigma, alpha, 
                            Lambda_prime, Lambda)

   # 2. 创建网格
   num_elem = 100
   num_nodes = num_elem + 1
   nodes = np.linspace(-1, 1, num_nodes)
   
   elem_connec1 = np.arange(0, num_elem)
   elem_connec2 = np.arange(1, num_nodes)
   connectivity = np.vstack((elem_connec1, elem_connec2)).T
   
   mesh = Mesh1D(nodes, connectivity)
   
   # 设置子域 (前一半是空气，后一半是多孔材料)
   air_elements = np.arange(0, int(num_elem / 2))
   porous_elements = np.arange(int(num_elem / 2), num_elem)
   mesh.set_subdomains({air: air_elements, porous: porous_elements})
   
   # 检查材料兼容性
   check_material_compability(mesh.subdomains)

   # 3. 设置FEM基函数
   order = 2  # 多项式阶数
   elements2node = mesh.get_mesh_coordinates()
   
   Pf_bases = []
   for mat, elems in mesh.subdomains.items():
       if mat.TYPE == 'Fluid':
           Pf_bases += [
               Lobbato1DElement('Pf', order, elements2node[elem]) 
               for elem in elems
           ]

   # 4. 自由度处理
   Helmholtz_dof_handler = GeneralDofHandler1D(['Pf'], Pf_bases)
   print(f"全局自由度数量: {Helmholtz_dof_handler.nb_dofs}")

   # 5. 组装全局矩阵
   Helmholtz_assember = HelmholtzAssembler(
       Helmholtz_dof_handler, mesh.subdomains, dtype=np.complex128
   )
   Helmholtz_assember.assembly_global_matrix(Pf_bases, 'Pf', omega)
   left_hand_matrix = Helmholtz_assember.get_global_matrix()

   # 6. 施加边界条件
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

   # 7. 求解
   linear_solver = LinearSolver(dof_handler=Helmholtz_dof_handler)
   linear_solver.solve(left_hand_matrix, right_hand_vec)
   sol = linear_solver.u

   # 8. 后处理
   post_processer = PostProcessField(mesh.nodes, r'1D Helmholtz (2000Hz)')
   post_processer.plot_sol((np.real(sol), 'FEM (p=2)', 'solid'))
   plt.show()

分步解释
--------

步骤1: 定义物理问题
~~~~~~~~~~~~~~~~~~~

首先定义仿真参数和材料属性。PyAcoustiX提供了多种预定义材料:

- ``Air``: 标准空气
- ``Fluid``: 自定义流体
- ``EquivalentFluid``: JCA等效流体模型
- ``LimpPorousMaterial``: Limp多孔材料
- ``PoroElasticMaterial``: Biot多孔弹性材料

步骤2: 创建网格
~~~~~~~~~~~~~~

网格是FEM仿真的基础。对于1D问题，我们使用 ``Mesh1D`` 类:

.. code-block:: python

   mesh = Mesh1D(nodes, connectivity)

``nodes`` 是节点坐标数组，``connectivity`` 定义了单元与节点的连接关系。

步骤3: 设置基函数
~~~~~~~~~~~~~~~~~

PyAcoustiX支持多种基函数:

- ``Lobbato1DElement``: 高阶Lobatto多项式
- ``Lagrange1DElement``: Lagrange多项式

步骤4-5: 自由度处理和矩阵组装
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

自由度处理器 ``GeneralDofHandler1D`` 负责管理全局自由度编号。

组装器 ``HelmholtzAssembler`` 根据物理问题类型组装全局刚度矩阵。

步骤6: 边界条件
~~~~~~~~~~~~~~~

支持多种边界条件类型:

- 本质边界条件 (Dirichlet)
- 自然边界条件 (Neumann)
- 阻抗边界条件

步骤7: 求解
~~~~~~~~~~~

``LinearSolver`` 封装了多种线性求解器，自动选择最适合的方法。

步骤8: 后处理
~~~~~~~~~~~~~

``PostProcessField`` 提供了绘图和结果分析功能。

下一步
------

- 查看 :doc:`tutorial_1d` 了解更详细的1D问题教程
- 查看 :doc:`tutorial_2d` 学习2D声学问题
- 查看 :doc:`../api_reference/index` 了解完整的API参考

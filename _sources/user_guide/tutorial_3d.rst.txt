教程3: 3D声学问题
=================

本教程介绍如何使用PyAcoustiX求解三维声学问题。3D仿真适用于真实世界的复杂声学场景。

问题描述
--------

考虑一个三维立方体腔体，尺寸为 :math:`L_x \times L_y \times L_z`。我们将求解其中的声场分布。

控制方程
~~~~~~~~

三维Helmholtz方程:

.. math::

   \nabla^2 p + k^2 p = 0

其中 :math:`\nabla^2 = \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} + \frac{\partial^2}{\partial z^2}`

完整代码
--------

.. code-block:: python

   import numpy as np
   import meshio
   
   from SAcouS.Mesh import Mesh3D, MeshReader
   from SAcouS.Materials import Air
   from SAcouS.acxfem import (
       Lobbato3DElement, GeneralDofHandler3D,
       HelmholtzAssembler3D, ApplyBoundaryConditions3D,
       LinearSolver, FESpace3D
   )

   def solve_3d_helmholtz():
       # ====================== 物理参数 ======================
       freq = 500  # Hz (低频，适合粗网格)
       omega = 2 * np.pi * freq
       
       air = Air('air')
       
       # ====================== 网格 (使用Gmsh) ======================
       # 建议：对于3D问题，使用Gmsh生成网格
       mesh_reader = MeshReader('cube.msh', dim=3, order=1)
       mesh = mesh_reader.get_mesh()
       mesh.set_subdomains({air: np.arange(mesh.nb_elems)})
       
       print(f"3D网格: {mesh.nb_nodes} 节点, {mesh.nb_elems} 单元")
       
       # ====================== FEM设置 ======================
       order = 1  # 3D问题通常使用线性单元以控制计算量
       elements2node = mesh.get_mesh_coordinates()
       
       Pf_bases = [
           Lobbato3DElement('Pf', order, elements2node[elem])
           for elem in range(mesh.nb_elems)
       ]
       
       dof_handler = GeneralDofHandler3D(['Pf'], Pf_bases)
       print(f"自由度数量: {dof_handler.nb_dofs}")
       
       # ====================== 组装 ======================
       assembler = HelmholtzAssembler3D(
           dof_handler, mesh.subdomains, dtype=np.complex128
       )
       assembler.assembly_global_matrix(Pf_bases, 'Pf', omega)
       K = assembler.get_global_matrix()
       
       print(f"矩阵大小: {K.shape}")
       print(f"非零元素: {K.nnz}")
       
       # ====================== 边界条件 ======================
       fe_space = FESpace3D(mesh, Pf_bases)
       F = np.zeros(assembler.nb_global_dofs, dtype=np.complex128)
       
       # 获取边界表面
       inlet_surface = mesh_reader.get_facet_by_physical('inlet')
       
       # 速度边界条件
       velocity_bc = {
           'type': 'fluid_velocity',
           'value': 1.0,
           'facets': inlet_surface
       }
       
       bc_applier = ApplyBoundaryConditions3D(mesh, fe_space, K, F, omega)
       bc_applier.apply_nature_bc(velocity_bc, var='Pf')
       
       # ====================== 求解 ======================
       print("求解中...")
       solver = LinearSolver(dof_handler=dof_handler)
       solver.solve(K, F)
       sol = solver.u
       
       # ====================== 后处理 ======================
       # 导出到VTK格式 (可用ParaView打开)
       meshio.write('results_3d.vtu', meshio.Mesh(
           points=mesh.nodes,
           cells=[('tetra', mesh.connectivity)],
           point_data={'Pressure': np.abs(sol), 
                      'Phase': np.angle(sol)}
       ))
       
       print("结果已保存到 results_3d.vtu")
       print(f"压力幅值范围: [{np.abs(sol).min():.4f}, {np.abs(sol).max():.4f}]")
       
       return sol

   if __name__ == "__main__":
       solution = solve_3d_helmholtz()

Gmsh网格生成
------------

对于3D问题，建议使用Gmsh生成四面体网格:

**cube.geo**

.. code-block:: cpp

   // 3D立方体腔体
   L = 1.0;  // 边长
   
   // 定义点
   Point(1) = {0, 0, 0};
   Point(2) = {L, 0, 0};
   Point(3) = {L, L, 0};
   Point(4) = {0, L, 0};
   Point(5) = {0, 0, L};
   Point(6) = {L, 0, L};
   Point(7) = {L, L, L};
   Point(8) = {0, L, L};
   
   // 定义线
   Line(1) = {1, 2};
   Line(2) = {2, 3};
   Line(3) = {3, 4};
   Line(4) = {4, 1};
   Line(5) = {5, 6};
   Line(6) = {6, 7};
   Line(7) = {7, 8};
   Line(8) = {8, 5};
   Line(9) = {1, 5};
   Line(10) = {2, 6};
   Line(11) = {3, 7};
   Line(12) = {4, 8};
   
   // 定义面
   Line Loop(1) = {1, 2, 3, 4};
   Plane Surface(1) = {1};
   Line Loop(2) = {5, 6, 7, 8};
   Plane Surface(2) = {2};
   Line Loop(3) = {1, 10, -5, -9};
   Plane Surface(3) = {3};
   Line Loop(4) = {2, 11, -6, -10};
   Plane Surface(4) = {4};
   Line Loop(5) = {3, 12, -7, -11};
   Plane Surface(5) = {5};
   Line Loop(6) = {4, 9, -8, -12};
   Plane Surface(6) = {6};
   
   // 定义体
   Surface Loop(1) = {1, 2, 3, 4, 5, 6};
   Volume(1) = {1};
   
   // 物理组
   Physical Volume("domain") = {1};
   Physical Surface("inlet") = {6};   // x=0面
   Physical Surface("outlet") = {4};  // x=L面
   Physical Surface("walls") = {1, 2, 3, 5};
   
   // 网格尺寸
   Mesh.CharacteristicLengthMax = 0.1;
   Mesh 3;
   Save "cube.msh";

计算资源考虑
------------

3D问题的计算成本显著高于2D:

**网格密度**

- 每波长至少6-10个线性单元
- 或使用高阶单元减少网格数

**内存需求**

自由度数量估计:

.. math::

   N_{dof} \approx N_{nodes} \times N_{variables}

对于3D四面体网格:

.. math::

   N_{nodes} \approx \frac{N_{elements}}{5}

**求解时间**

- 直接求解器: 适用于 < 100,000 DOF
- 迭代求解器: 适用于更大规模问题

并行计算
--------

对于大规模3D问题，可以使用并行计算:

.. code-block:: python

   from SAcouS.acxfem import ParallelSolver
   
   solver = ParallelSolver(dof_handler=dof_handler, n_processes=4)
   solver.solve(K, F)

后处理与可视化
--------------

**ParaView**

推荐使用 `ParaView <https://www.paraview.org/>`_ 可视化3D结果:

1. 打开 ``results_3d.vtu`` 文件
2. 选择变量 (Pressure, Phase)
3. 使用切片、等值面等工具分析

**Gmsh**

也可以直接在Gmsh中查看:

.. code-block:: python

   from SAcouS.PostProcess import save_plot
   save_plot(mesh, sol, 'Pressure', 'results.msh', engine='gmsh')

实际应用
--------

3D声学仿真的典型应用:

1. **汽车内部噪声**: 发动机噪声在车厢内的传播
2. **建筑声学**: 音乐厅、会议室的声学设计
3. **消声器设计**: 排气系统的降噪优化
4. **扬声器设计**: 音箱腔体的声学特性

练习
----

1. **网格收敛性**: 使用不同密度的网格，观察结果变化
2. **频率扫描**: 在多个频率下求解，分析共振频率
3. **添加吸声材料**: 在腔体壁面添加多孔材料
4. **复杂几何**: 尝试非规则几何形状

下一步
------

- :doc:`../examples/index`: 查看更多示例
- :doc:`../theory/helmholtz`: 深入学习声学理论
- :doc:`../api_reference/index`: 查看完整API文档

教程2: 2D声学问题
=================

本教程介绍如何使用PyAcoustiX求解二维声学问题。我们将模拟声波在一个2D腔体中的传播。

问题描述
--------

考虑一个二维矩形腔体，尺寸为 :math:`L_x \times L_y`。在左边界施加速度激励，其他边界为刚性壁(速度为零)。

几何与边界条件
~~~~~~~~~~~~~~

.. code-block:: text

   +------------------+
   |                  |  ^ y
   |  空气腔体         |  |
   |                  |  |
   +------------------+  +----> x
   
   左边界: 速度激励 v = v₀
   其他边界: 刚性壁 (∂p/∂n = 0)

控制方程
~~~~~~~~

二维Helmholtz方程:

.. math::

   \nabla^2 p + k^2 p = 0

其中 :math:`\nabla^2 = \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2}` 是拉普拉斯算子。

完整代码
--------

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   
   from SAcouS.Mesh import Mesh2D, MeshReader
   from SAcouS.Materials import Air
   from SAcouS.acxfem import (
       Lobbato2DElement, GeneralDofHandler2D, 
       HelmholtzAssembler2D, ApplyBoundaryConditions2D,
       LinearSolver, FESpace2D
   )
   from SAcouS.PostProcess import plot_field

   def solve_2d_helmholtz():
       # ====================== 物理参数 ======================
       freq = 1000  # Hz
       omega = 2 * np.pi * freq
       
       # 材料
       air = Air('air')
       
       # ====================== 网格 ======================
       # 方法1: 使用Gmsh文件
       # mesh_reader = MeshReader('cavity.msh', dim=2, order=1)
       # mesh = mesh_reader.get_mesh()
       
       # 方法2: 手动创建简单网格
       Lx, Ly = 1.0, 0.5  # 腔体尺寸
       nx, ny = 20, 10    # 网格划分
       
       # 创建节点
       x = np.linspace(0, Lx, nx + 1)
       y = np.linspace(0, Ly, ny + 1)
       xv, yv = np.meshgrid(x, y)
       nodes = np.vstack([xv.ravel(), yv.ravel()]).T
       
       # 创建三角形单元连接关系
       elem_connect = []
       for j in range(ny):
           for i in range(nx):
               n1 = j * (nx + 1) + i
               n2 = n1 + 1
               n3 = n1 + nx + 1
               n4 = n3 + 1
               # 两个三角形组成一个矩形
               elem_connect.append([n1, n2, n3])
               elem_connect.append([n2, n4, n3])
       elem_connect = np.array(elem_connect)
       
       # 定义边界边
       edge_connect = []
       # 底边
       for i in range(nx):
           edge_connect.append([i, i + 1])
       # 右边
       for j in range(ny):
           edge_connect.append([(j + 1) * (nx + 1) - 1, (j + 2) * (nx + 1) - 1])
       # 顶边
       for i in range(nx):
           edge_connect.append([(ny + 1) * (nx + 1) - 1 - i, (ny + 1) * (nx + 1) - 2 - i])
       # 左边
       for j in range(ny):
           edge_connect.append([(ny - j) * (nx + 1), (ny - j - 1) * (nx + 1)])
       edge_connect = np.array(edge_connect)
       
       mesh = Mesh2D(nodes, elem_connect, edge_connect)
       mesh.set_subdomains({air: np.arange(len(elem_connect))})
       
       print(f"网格: {mesh.nb_nodes} 节点, {mesh.nb_elems} 单元")
       
       # ====================== FEM设置 ======================
       order = 2
       elements2node = mesh.get_mesh_coordinates()
       
       Pf_bases = [
           Lobbato2DElement('Pf', order, elements2node[elem])
           for elem in range(mesh.nb_elems)
       ]
       
       dof_handler = GeneralDofHandler2D(['Pf'], Pf_bases)
       print(f"自由度数量: {dof_handler.nb_dofs}")
       
       # ====================== 组装 ======================
       assembler = HelmholtzAssembler2D(
           dof_handler, mesh.subdomains, dtype=np.complex128
       )
       assembler.assembly_global_matrix(Pf_bases, 'Pf', omega)
       K = assembler.get_global_matrix()
       
       # ====================== 边界条件 ======================
       fe_space = FESpace2D(mesh, Pf_bases)
       F = np.zeros(assembler.nb_global_dofs, dtype=np.complex128)
       
       # 找到左边界节点
       left_boundary_nodes = np.where(np.abs(nodes[:, 0]) < 1e-10)[0]
       
       # 速度边界条件
       velocity_bc = {
           'type': 'fluid_velocity',
           'value': 1.0,  # 单位速度
           'nodes': left_boundary_nodes
       }
       
       bc_applier = ApplyBoundaryConditions2D(mesh, fe_space, K, F, omega)
       bc_applier.apply_nature_bc(velocity_bc, var='Pf')
       
       # ====================== 求解 ======================
       solver = LinearSolver(dof_handler=dof_handler)
       solver.solve(K, F)
       sol = solver.u
       
       # ====================== 后处理 ======================
       # 绘制压力场
       plot_field(mesh, np.real(sol), 'Pressure Field (Real Part)', 
                  quantity='Pressure', unit='Pa')
       
       # 计算声压级
       p_ref = 20e-6  # 参考声压
       SPL = 20 * np.log10(np.abs(sol) / p_ref)
       plot_field(mesh, SPL, 'Sound Pressure Level', 
                  quantity='SPL', unit='dB')
       
       plt.show()
       
       return sol

   if __name__ == "__main__":
       solution = solve_2d_helmholtz()

使用Gmsh网格
------------

对于复杂几何，建议使用 `Gmsh <http://gmsh.info/>`_ 生成网格:

1. 创建 ``.geo`` 文件:

.. code-block:: cpp

   // cavity.geo
   Lx = 1.0;
   Ly = 0.5;
   
   Point(1) = {0, 0, 0};
   Point(2) = {Lx, 0, 0};
   Point(3) = {Lx, Ly, 0};
   Point(4) = {0, Ly, 0};
   
   Line(1) = {1, 2};
   Line(2) = {2, 3};
   Line(3) = {3, 4};
   Line(4) = {4, 1};
   
   Line Loop(1) = {1, 2, 3, 4};
   Plane Surface(1) = {1};
   
   Physical Surface("domain") = {1};
   Physical Line("left") = {4};
   Physical Line("other") = {1, 2, 3};
   
   Mesh 2;
   Save "cavity.msh";

2. 在PyAcoustiX中读取:

.. code-block:: python

   mesh_reader = MeshReader('cavity.msh', dim=2, order=1)
   mesh = mesh_reader.get_mesh()
   
   # 获取物理组
   left_edges = mesh_reader.get_facet_by_physical('left')
   other_edges = mesh_reader.get_facet_by_physical('other')

高阶单元
--------

PyAcoustiX支持高阶单元以获得更高精度:

.. code-block:: python

   # 使用3阶单元
   order = 3
   Pf_bases = [
       Lobbato2DElement('Pf', order, elements2node[elem])
       for elem in range(mesh.nb_elems)
   ]

高阶单元的优势:

- 相同网格数下精度更高
- 更适合处理高频问题
- 更好的波数分辨率

但计算成本也更高，需要权衡。

边界条件类型
------------

PyAcoustiX支持多种2D边界条件:

**速度边界**

.. code-block:: python

   bc = {'type': 'fluid_velocity', 'value': 1.0, 'facets': edge_indices}

**压力边界 (Dirichlet)**

.. code-block:: python

   bc = {'type': 'pressure', 'value': 1.0, 'nodes': node_indices}

**阻抗边界**

.. code-block:: python

   bc = {'type': 'impedance', 'value': 416.5, 'facets': edge_indices}

可视化
------

PyAcoustiX提供了多种可视化选项:

**使用Matplotlib**

.. code-block:: python

   from SAcouS.PostProcess import plot_field
   plot_field(mesh, sol, 'Title')

**导出到Gmsh**

.. code-block:: python

   from SAcouS.PostProcess import save_plot
   save_plot(mesh, sol, 'Pressure', 'results.msh', engine='gmsh')

**导出到VTK (ParaView)**

.. code-block:: python

   import meshio
   meshio.write('results.vtk', meshio.Mesh(
       points=mesh.nodes,
       cells=[('triangle', mesh.connectivity)],
       point_data={'Pressure': sol}
   ))

练习
----

1. **改变频率**: 观察不同频率下的压力分布模式
2. **添加多孔材料**: 在腔体的一部分添加吸声材料
3. **使用不同网格**: 比较粗网格和细网格的结果
4. **改变单元阶数**: 比较线性单元和高阶单元的精度

下一步
------

- :doc:`tutorial_3d`: 学习3D声学问题
- :doc:`../theory/helmholtz`: 深入了解Helmholtz方程
- :doc:`../api_reference/index`: 查看完整的API文档

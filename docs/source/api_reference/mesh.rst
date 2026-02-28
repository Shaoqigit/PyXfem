Mesh模块
========

Mesh模块提供了1D、2D和3D网格的创建、读取和操作功能。

.. module:: SAcouS.Mesh

基类
----

.. autoclass:: SAcouS.Mesh.BaseMesh
   :members:
   :undoc-members:
   :show-inheritance:

1D网格
------

.. autoclass:: SAcouS.Mesh.Mesh1D
   :members:
   :undoc-members:
   :show-inheritance:

2D网格
------

.. autoclass:: SAcouS.Mesh.Mesh2D
   :members:
   :undoc-members:
   :show-inheritance:

3D网格
------

.. autoclass:: SAcouS.Mesh.Mesh3D
   :members:
   :undoc-members:
   :show-inheritance:

网格读取器
----------

.. autoclass:: SAcouS.Mesh.MeshReader
   :members:
   :undoc-members:
   :show-inheritance:

工具函数
--------

.. autofunction:: SAcouS.Mesh.mesh_constructor

使用示例
--------

创建1D网格
~~~~~~~~~~

.. code-block:: python

   import numpy as np
   from SAcouS.Mesh import Mesh1D
   
   # 定义节点
   nodes = np.linspace(0, 1, 101)  # 100个单元
   
   # 定义连接关系
   elem_connect = np.array([[i, i+1] for i in range(100)])
   
   # 创建网格
   mesh = Mesh1D(nodes, elem_connect)
   
   # 设置子域
   air_elements = np.arange(0, 50)
   porous_elements = np.arange(50, 100)
   mesh.set_subdomains({
       air: air_elements,
       porous: porous_elements
   })

读取Gmsh网格
~~~~~~~~~~~~

.. code-block:: python

   from SAcouS.Mesh import MeshReader
   
   # 读取2D网格
   reader = MeshReader('mesh.msh', dim=2, order=1)
   mesh = reader.get_mesh()
   
   # 获取物理组
   inlet_edges = reader.get_facet_by_physical('inlet')
   wall_edges = reader.get_facet_by_physical('wall')

网格操作
~~~~~~~~

.. code-block:: python

   # 获取单元坐标
   elem_coords = mesh.get_mesh_coordinates()
   
   # 网格细化
   mesh.refine_mesh(times=1)
   
   # 绘制网格
   mesh.plotmesh(withnode=True, withnodeid=True)
   
   # 计算边界法向量
   normal = mesh.compute_normal(edge)

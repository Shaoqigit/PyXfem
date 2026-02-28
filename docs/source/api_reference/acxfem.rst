acxfem模块
==========

acxfem模块是PyAcoustiX的核心FEM实现，提供了基函数、组装、求解等完整功能。

.. module:: SAcouS.acxfem

基函数 (Basis)
--------------

.. autoclass:: SAcouS.acxfem.Basis.Base1DElement
   :members:
   :undoc-members:

1D基函数
~~~~~~~~

.. autoclass:: SAcouS.acxfem.Basis.Lobbato1DElement
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.Basis.Helmholtz1DElement
   :members:
   :undoc-members:

2D/3D基函数
~~~~~~~~~~~

.. autoclass:: SAcouS.acxfem.Basis.BaseNDElement
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.Basis.Lagrange2DTriElement
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.Basis.Helmholtz2DElement
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.Basis.Lagrange2DQuadElement
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.Basis.Lagrange3DTetraElement
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.Basis.Helmholtz3DElement
   :members:
   :undoc-members:

自由度处理 (DofHandler)
-----------------------

.. autoclass:: SAcouS.acxfem.DofHandler.DofHandler1D
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.DofHandler.DofHandler1DMutipleVariable
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.DofHandler.GeneralDofHandler1D
   :members:
   :undoc-members:

有限元空间
~~~~~~~~~~

.. autoclass:: SAcouS.acxfem.DofHandler.FESpace
   :members:
   :undoc-members:

组装 (Assembly)
---------------

.. autoclass:: SAcouS.acxfem.Assembly.Assembler
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.Assembly.Assembler4Biot
   :members:
   :undoc-members:

物理组装器 (PhysicAssembler)
----------------------------

.. autoclass:: SAcouS.acxfem.PhysicAssembler.BaseAssembler
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.PhysicAssembler.HelmholtzAssembler
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.PhysicAssembler.BiotAssembler
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.PhysicAssembler.CouplingAssember
   :members:
   :undoc-members:

求解器 (Solver)
---------------

.. autoclass:: SAcouS.acxfem.Solver.BaseSolver
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.Solver.LinearSolver
   :members:
   :undoc-members:

.. autoclass:: SAcouS.acxfem.Solver.AdmittanceSolver
   :members:
   :undoc-members:

边界条件 (BCsImpose)
--------------------

.. autoclass:: SAcouS.acxfem.BCsImpose.ApplyBoundaryConditions
   :members:
   :undoc-members:

数值积分 (Quadratures)
----------------------

.. autofunction:: SAcouS.acxfem.Quadratures.gauss_legendre_1d_points

.. autofunction:: SAcouS.acxfem.Quadratures.gauss_legendre_1d_weights

.. autofunction:: SAcouS.acxfem.Quadratures.gauss_legendre_2d_tri_points

.. autofunction:: SAcouS.acxfem.Quadratures.gauss_legendre_2d_tri_weights

.. autofunction:: SAcouS.acxfem.Quadratures.gauss_legendre_3d_tetra_points

.. autofunction:: SAcouS.acxfem.Quadratures.gauss_legendre_3d_tetra_weights

.. autofunction:: SAcouS.acxfem.Quadratures.get_quadrature_points_weights

工具函数
--------

.. autofunction:: SAcouS.acxfem.Utilities.check_material_compability

.. autofunction:: SAcouS.acxfem.Utilities.display_matrix_in_array

.. autofunction:: SAcouS.acxfem.Utilities.plot_matrix_partten

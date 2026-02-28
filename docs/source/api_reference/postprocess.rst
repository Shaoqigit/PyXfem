PostProcess模块
===============

PostProcess模块提供了结果可视化和分析工具。

.. module:: SAcouS.PostProcess

基类
----

.. autoclass:: SAcouS.PostProcess.BasePostProcess
   :members:
   :undoc-members:
   :show-inheritance:

场后处理
--------

.. autoclass:: SAcouS.PostProcess.PostProcessField
   :members:
   :undoc-members:
   :show-inheritance:

频率响应后处理
--------------

.. autoclass:: SAcouS.PostProcess.PostProcessFRF
   :members:
   :undoc-members:
   :show-inheritance:

工具函数
--------

.. autofunction:: SAcouS.PostProcess.plot_field

.. autofunction:: SAcouS.PostProcess.save_gmsh

.. autofunction:: SAcouS.PostProcess.save_plot

.. autofunction:: SAcouS.PostProcess.read_solution

使用示例
--------

1D场可视化
~~~~~~~~~~

.. code-block:: python

   from SAcouS.PostProcess import PostProcessField
   import matplotlib.pyplot as plt
   
   # 创建后处理器
   post = PostProcessField(mesh.nodes, 'Pressure Distribution')
   
   # 绘制实部
   post.plot_sol(
       (np.real(sol), 'Real part', 'solid'),
       (np.imag(sol), 'Imaginary part', 'dashed')
   )
   plt.show()

频率响应
~~~~~~~~

.. code-block:: python

   from SAcouS.PostProcess import PostProcessFRF
   
   # 频率扫描结果
   freqs = np.linspace(100, 2000, 100)
   results = [...]  # 每个频率的结果
   
   # 创建FRF后处理器
   post_frf = PostProcessFRF(freqs, 'Frequency Response', 
                             acoustic_indicator='SPL')
   
   # 绘制声压级
   post_frf.plot_sol(
       (results, 'FEM', 'solid'),
       (reference, 'Reference', 'dashed')
   )
   
   # 保存结果
   post_frf.save_sol((results, 'FEM'), file_name='frf_results.txt')

2D/3D可视化
~~~~~~~~~~~

.. code-block:: python

   from SAcouS.PostProcess import plot_field
   
   # 2D场可视化
   plot_field(mesh, np.real(sol), 'Pressure Field',
              quantity='Pressure', unit='Pa')
   
   # 计算并绘制SPL
   p_ref = 20e-6
   SPL = 20 * np.log10(np.abs(sol) / p_ref)
   plot_field(mesh, SPL, 'Sound Pressure Level',
              quantity='SPL', unit='dB')

导出结果
~~~~~~~~

.. code-block:: python

   import meshio
   
   # 导出到VTK (ParaView)
   meshio.write('results.vtu', meshio.Mesh(
       points=mesh.nodes,
       cells=[('triangle', mesh.connectivity)],
       point_data={
           'Pressure': sol,
           'SPL': 20 * np.log10(np.abs(sol) / 20e-6)
       }
   ))
   
   # 导出到Gmsh
   from SAcouS.PostProcess import save_plot
   save_plot(mesh, sol, 'Pressure', 'results.msh', engine='gmsh')

误差计算
~~~~~~~~

.. code-block:: python

   # 计算与解析解的误差
   error = post.compute_error(fem_sol, analytical_sol)
   print(f"相对L2误差: {error:.2%}")

interface模块
=============

interface模块提供了输入文件解析和求解设置功能。

.. module:: SAcouS.interface

输入文件解析
------------

.. autoclass:: SAcouS.interface.parser.BaseParser
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: SAcouS.interface.parser.AcoustiXPaser
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: SAcouS.interface.parser.ParserFactory
   :members:
   :undoc-members:
   :show-inheritance:

求解设置
--------

.. autoclass:: SAcouS.interface.sol_setup.PyAcoustiXSetuper
   :members:
   :undoc-members:
   :show-inheritance:

使用示例
--------

解析输入文件
~~~~~~~~~~~~

.. code-block:: python

   from SAcouS.interface.sol_setup import PyAcoustiXSetuper
   
   # 创建设置器
   setuper = PyAcoustiXSetuper()
   
   # 解析输入文件
   setuper.parse_input('simulation.axi')
   
   # 获取求解信息
   info = setuper.sol_info
   print(f"物理类型: {info['physics']}")
   print(f"频率: {info['frequency']} Hz")
   
   # 显示欢迎信息
   setuper.welcome()

输入文件格式
------------

``.axi`` 文件采用INI格式:

.. code-block:: ini

   [Mesh]
   file = mesh.msh
   dimension = 2
   order = 1
   
   [Materials]
   air = Air, air
   foam = EquivalentFluid, 0.98, 3750, 1.17, 742e-6, 110e-6
   
   [Physics]
   type = Helmholtz
   frequency = 1000
   
   [BoundaryConditions]
   inlet = velocity, 1.0
   outlet = impedance, 416.5
   wall = rigid
   
   [Solver]
   type = direct
   
   [Output]
   file = results.vtk
   format = vtk

节说明
~~~~~~

**[Mesh]**

- ``file``: 网格文件路径
- ``dimension``: 空间维度 (1, 2, 3)
- ``order``: 网格阶数 (1或2)

**[Materials]**

格式: ``名称 = 类型, 参数...``

支持的类型:

- ``Air``: 标准空气
- ``Fluid``: 自定义流体 (密度, 声速)
- ``EquivalentFluid``: JCA模型 (孔隙率, 流阻, 曲折度, ...)
- ``LimpPorousMaterial``: Limp模型
- ``PoroElasticMaterial``: Biot模型

**[Physics]**

- ``type``: 物理问题类型 (Helmholtz, Biot, ...)
- ``frequency``: 频率 (Hz)

**[BoundaryConditions]**

格式: ``边界名称 = 类型, 值``

支持的类型:

- ``velocity``: 速度边界
- ``pressure``: 压力边界
- ``impedance``: 阻抗边界
- ``rigid``: 刚性壁

**[Solver]**

- ``type``: 求解器类型 (direct, iterative)
- ``tolerance``: 迭代容差 (迭代求解器)
- ``max_iter``: 最大迭代次数

**[Output]**

- ``file``: 输出文件名
- ``format``: 输出格式 (vtk, msh, txt)

命令行接口
----------

.. code-block:: bash

   # 运行输入文件
   sacous -i input.axi
   
   # 或
   sacous --input_file=input.axi

编程接口
--------

自定义解析器
~~~~~~~~~~~~

.. code-block:: python

   from SAcouS.interface.parser import BaseParser
   
   class MyParser(BaseParser):
       def parse_materials(self, section):
           # 自定义材料解析
           pass
       
       def parse_boundary_conditions(self, section):
           # 自定义边界条件解析
           pass
   
   parser = MyParser()
   config = parser.parse('input.axi')

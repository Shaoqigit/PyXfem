命令行使用
==========

PyAcoustiX提供了命令行工具 ``sacous``，可以方便地运行仿真和进行交互式计算。

基本用法
--------

安装完成后，``sacous`` 命令会自动添加到系统路径。

交互式控制台
~~~~~~~~~~~~

不带参数运行 ``sacous`` 会启动交互式Python控制台:

.. code-block:: bash

   $ sacous
   Welcome to PyAcoustiX console!
   PyAcoustiX 0.9.13
   base on python 3.10.13 (main, Sep 11 2023, 13:44:35) [GCC 11.2.0] on linux
   Type "help", "copyright" or "license" for more information.
   Author: Shaoqiwu@outlook.com
   (PyAcoustiXInteractiveConsole)
   >>> import numpy as np
   >>> from SAcouS import *

在交互式控制台中，你可以:

- 直接导入 ``SAcouS`` 模块
- 使用NumPy、SciPy等科学计算库
- 使用上下箭头浏览历史命令
- 使用Tab键自动补全

运行输入文件
~~~~~~~~~~~~

使用 ``-i`` 或 ``--input_file`` 参数指定输入文件:

.. code-block:: bash

   sacous -i simulation.axi

输入文件格式 (.axi)
--------------------

``.axi`` 文件是PyAcoustiX的标准输入文件格式，采用类似INI的语法:

示例输入文件
~~~~~~~~~~~~

.. code-block:: ini

   [Mesh]
   file = mesh.msh
   dimension = 2
   
   [Materials]
   air = Air, air_properties
   porous = EquivalentFluid, 0.98, 3750, 1.17, 742e-6, 110e-6
   
   [Physics]
   type = Helmholtz
   frequency = 1000
   
   [BoundaryConditions]
   inlet = velocity, 1.0
   outlet = impedance, 416.5
   
   [Solver]
   type = direct
   
   [Output]
   file = results.vtk

节说明
~~~~~~

**[Mesh]**
定义网格文件和维度。

**[Materials]**
定义材料属性。格式为: ``名称 = 类型, 参数...``

**[Physics]**
定义物理问题类型和参数。

**[BoundaryConditions]**
定义边界条件。

**[Solver]**
定义求解器设置。

**[Output]**
定义输出文件。

高级用法
--------

批处理模式
~~~~~~~~~~

可以编写shell脚本批量运行多个频率:

.. code-block:: bash

   #!/bin/bash
   for freq in 100 200 500 1000 2000; do
       sed "s/frequency = .*/frequency = $freq/" template.axi > temp.axi
       sacous -i temp.axi
       mv results.vtk results_${freq}Hz.vtk
   done

与Python脚本集成
~~~~~~~~~~~~~~~~

也可以在Python脚本中调用 ``sacous`` 的功能:

.. code-block:: python

   from SAcouS.interface.sol_setup import PyAcoustiXSetuper
   
   sol_setuper = PyAcoustiXSetuper()
   sol_setuper.parse_input('input.axi')
   # ... 进一步处理

故障排除
--------

命令未找到
~~~~~~~~~~

如果提示 ``sacous: command not found``:

1. 确保PyAcoustiX已正确安装
2. 检查Python脚本目录是否在PATH中:
   
   .. code-block:: bash
   
      pip show pyacoustix | grep Location
      
3. 尝试使用Python模块方式运行:
   
   .. code-block:: bash
   
      python -m SAcouS

输入文件错误
~~~~~~~~~~~~

如果输入文件解析失败:

1. 检查文件语法是否正确
2. 确保所有必需的节都存在
3. 查看错误信息中的行号提示

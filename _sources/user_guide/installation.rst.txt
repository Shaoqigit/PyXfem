安装指南
========

系统要求
--------

- Python 3.8 或更高版本
- NumPy, SciPy, Matplotlib
- meshio (用于网格I/O)

通过pip安装
------------

最简单的安装方式是通过PyPI:

.. code-block:: bash

   pip install pyacoustix

这将自动安装所有必需的依赖项。

从源码安装
----------

如果你想获取最新开发版本或贡献代码，可以从GitHub克隆并安装:

.. code-block:: bash

   git clone https://github.com/Shaoqigit/PyXfem.git
   cd PyXfem
   pip install -e .

使用 ``-e`` 标志以"可编辑模式"安装，这样你可以修改代码而无需重新安装。

Docker安装
----------

我们也提供了Docker镜像，方便在不同环境中使用:

.. code-block:: bash

   docker pull shaoqiwu/pyacoustix
   docker run -it shaoqiwu/pyacoustix

依赖项
------

必需依赖
~~~~~~~~

- **numpy** (>=1.23.4): 数值计算基础
- **scipy** (>=1.10.0): 科学计算和稀疏矩阵
- **matplotlib** (>=3.7.0): 绘图和可视化
- **meshio** (>=5.0.0): 网格文件I/O

可选依赖
~~~~~~~~

- **pymls**: 用于特定材料模型
- **mpi4py**: 用于并行计算
- **petsc4py**: 用于大规模问题的求解器

验证安装
--------

安装完成后，可以通过以下方式验证:

.. code-block:: python

   >>> from SAcouS import *
   >>> print("PyAcoustiX 安装成功!")

或者运行测试套件:

.. code-block:: bash

   python run_test.py

常见问题
--------

ImportError: No module named 'SAcouS'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

确保你已经正确安装了包。如果使用的是虚拟环境，请激活它:

.. code-block:: bash

   source venv/bin/activate  # Linux/Mac
   venv\Scripts\activate      # Windows

版本冲突
~~~~~~~~

如果遇到依赖版本冲突，建议使用虚拟环境:

.. code-block:: bash

   python -m venv venv
   source venv/bin/activate
   pip install pyacoustix

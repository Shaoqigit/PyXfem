测试指南
========

本章节介绍PyAcoustiX的测试策略和实践。

测试结构
--------

测试文件位于 ``tests/`` 目录:

.. code-block:: text

   tests/
   ├── test_2D_air.py          # 2D空气问题
   ├── test_biot_equation.py   # Biot方程
   ├── test_tmm_biot.py        # TMM测试
   ├── test_mor_porous.py      # MoR测试
   └── ...

运行测试
--------

运行所有测试
~~~~~~~~~~~~

.. code-block:: bash

   python run_test.py

运行特定测试
~~~~~~~~~~~~

.. code-block:: bash

   python tests/test_2D_air.py

使用pytest
~~~~~~~~~~

.. code-block:: bash

   pip install pytest
   pytest tests/

测试类型
--------

单元测试
~~~~~~~~

测试单个函数或类:

.. code-block:: python

   def test_material_creation():
       air = Air('air')
       assert air.rho == 1.213
       assert air.c == 343.0

集成测试
~~~~~~~~

测试多个组件的协同工作:

.. code-block:: python

   def test_helmholtz_solver():
       # 创建网格
       mesh = Mesh1D(nodes, connectivity)
       
       # 设置材料
       air = Air('air')
       mesh.set_subdomains({air: elements})
       
       # 组装和求解
       assembler = HelmholtzAssembler(...)
       assembler.assembly_global_matrix(...)
       
       # 验证结果
       assert np.allclose(sol, expected, rtol=1e-3)

验证测试
~~~~~~~~

与解析解对比:

.. code-block:: python

   def test_kundt_tube():
       # FEM解
       fem_sol = solve_fem(...)
       
       # 解析解
       ana_sol = analytical_solution(...)
       
       # 计算误差
       error = compute_error(fem_sol, ana_sol)
       assert error < 0.01  # 1%误差容限

编写测试
--------

测试模板
~~~~~~~~

.. code-block:: python

   import numpy as np
   from SAcouS import *
   
   def test_feature_name():
       """测试功能描述"""
       # 准备
       input_data = ...
       expected_output = ...
       
       # 执行
       actual_output = function_to_test(input_data)
       
       # 验证
       assert np.allclose(actual_output, expected_output)
       
       print("Test passed!")
   
   if __name__ == "__main__":
       test_feature_name()

测试命名
~~~~~~~~

- 测试文件: ``test_*.py``
- 测试函数: ``test_*``
- 描述性名称，说明测试内容

测试数据
~~~~~~~~

- 使用简单、可重现的测试数据
- 避免硬编码大数组
- 使用 ``np.testing`` 进行数组比较

断言方法
~~~~~~~~

.. code-block:: python

   # 数值比较
   assert np.allclose(a, b, rtol=1e-5, atol=1e-8)
   
   # 数组相等
   np.testing.assert_array_equal(a, b)
   
   # 近似相等
   np.testing.assert_almost_equal(a, b, decimal=5)
   
   # 引发异常
   with pytest.raises(ValueError):
       invalid_operation()

持续集成
--------

GitHub Actions配置:

.. code-block:: yaml

   name: Tests
   
   on: [push, pull_request]
   
   jobs:
     test:
       runs-on: ubuntu-latest
       steps:
         - uses: actions/checkout@v2
         - name: Set up Python
           uses: actions/setup-python@v2
           with:
             python-version: 3.9
         - name: Install dependencies
           run: |
             pip install -e .
             pip install pytest
         - name: Run tests
           run: pytest tests/

测试覆盖率
----------

生成覆盖率报告:

.. code-block:: bash

   pip install pytest-cov
   pytest --cov=SAcouS tests/

目标覆盖率: > 80%

调试测试
--------

使用pdb调试:

.. code-block:: python

   import pdb; pdb.set_trace()

使用pytest的调试模式:

.. code-block:: bash

   pytest --pdb tests/test_example.py

最佳实践
--------

1. **独立性**: 每个测试应该独立运行
2. **可重复性**: 使用固定随机种子
3. **快速**: 单个测试应该快速完成
4. **清晰**: 测试失败时信息明确
5. **覆盖**: 测试正常和异常情况

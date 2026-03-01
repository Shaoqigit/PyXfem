贡献指南
========

感谢你对PyAcoustiX的贡献！

如何贡献
--------

报告问题
~~~~~~~~

如果你发现了bug或有功能建议:

1. 检查是否已有相关issue
2. 创建新的issue，包含:
   - 问题描述
   - 复现步骤
   - 期望行为
   - 实际行为
   - 环境信息 (Python版本, OS等)

提交代码
~~~~~~~~

1. Fork仓库
2. 创建feature分支
3. 编写代码和测试
4. 提交Pull Request

开发环境设置
------------

.. code-block:: bash

   # 克隆仓库
   git clone https://github.com/Shaoqigit/PyXfem.git
   cd PyXfem
   
   # 创建虚拟环境
   python -m venv venv
   source venv/bin/activate
   
   # 安装开发依赖
   pip install -e ".[dev]"
   
   # 运行测试
   python run_test.py

待开发功能列表
--------------

以下是PyAcoustiX计划实现但尚未完成的功能，欢迎社区贡献！

扩展有限元方法 (X-FEM)
~~~~~~~~~~~~~~~~~~~~~~

- **难度**: 高
- **描述**: 实现扩展有限元方法，用于处理裂纹扩展、材料界面等不连续问题
- **技能要求**: FEM理论、Python数值计算
- **相关文件**: `SAcouS/acxfem/Basis.py`, `SAcouS/acxfem/Assembly.py`

Ghost Node 技术
~~~~~~~~~~~~~~~

- **难度**: 中高
- **描述**: 实现Ghost Node技术，用于处理多材料界面和复杂几何
- **技能要求**: 计算几何、FEM基础
- **相关文件**: `SAcouS/Mesh.py`, `SAcouS/acxfem/DofHandler.py`

完整的 FEM-TMM 耦合
~~~~~~~~~~~~~~~~~~~

- **难度**: 中
- **描述**: 深度集成传递矩阵法(TMM)与有限元方法，实现真正的耦合计算
- **技能要求**: 声学理论、矩阵运算
- **相关文件**: `SAcouS/acxtmm/`, `SAcouS/acxfem/Assembly.py`

壳单元 (Shell Element)
~~~~~~~~~~~~~~~~~~~~~~

- **难度**: 中高
- **描述**: 为3D薄板结构实现高效壳单元
- **技能要求**: 结构力学、FEM高级主题
- **相关文件**: `SAcouS/acxfem/Basis.py`

完美匹配层 (PML)
~~~~~~~~~~~~~~~~

- **难度**: 中
- **描述**: 实现PML吸收边界条件，用于无限域声学问题
- **技能要求**: 声学理论、边界条件处理
- **相关文件**: `SAcouS/acxfem/BCsImpose.py`

无限元 (Infinite Element)
~~~~~~~~~~~~~~~~~~~~~~~~~

- **难度**: 中高
- **描述**: 实现无限元方法，用于自由场辐射和散射问题
- **技能要求**: 波动理论、数值方法
- **相关文件**: `SAcouS/acxfem/Basis.py`, `SAcouS/acxfem/Assembly.py`

并行计算支持
~~~~~~~~~~~~

- **难度**: 高
- **描述**: 基于MPI实现频率扫描的并行计算
- **技能要求**: 并行计算、MPI
- **相关文件**: `SAcouS/acxfem/Solver.py`, `run_test.py`

JIT 编译加速
~~~~~~~~~~~~

- **难度**: 中
- **描述**: 使用Numba对单元矩阵计算进行即时编译加速
- **技能要求**: Numba, Python性能优化
- **相关文件**: `SAcouS/acxfem/Assembly.py`, `SAcouS/acxfem/Basis.py`

阻尼系统的 MoR
~~~~~~~~~~~~~~

- **难度**: 高
- **描述**: 针对多孔材料实现RB+EIM降阶算法
- **技能要求**: 降阶方法、线性代数
- **相关文件**: `SAcouS/acxmor/`

有限导纳法 (FAM)
~~~~~~~~~~~~~~~~

- **难度**: 中高
- **描述**: 实现FAM方法，用于高精度多层结构配置
- **技能要求**: 声学理论、多层结构建模
- **相关文件**: `SAcouS/acxtmm/`

代码重构任务
~~~~~~~~~~~~

- **1D Lobatto形函数重构**: 使其API与2D/3D实现保持一致
- **材料类修正**: 完善Limp和BiotMaterial类
- **3D网格可视化**: 完成3D网格的绘制功能
- **3D二次单元**: 完整的3D二次单元支持

贡献流程
--------

1. **选择功能**: 从上面的列表中选择感兴趣的功能
2. **创建Issue**: 在GitHub上创建issue，说明你要开发的功能
3. **讨论方案**: 与维护者讨论实现方案
4. **开发实现**: 编写代码和测试
5. **提交PR**: 提交Pull Request，包含详细说明

代码规范
--------

Python风格
~~~~~~~~~~

- 遵循PEP 8
- 使用4空格缩进
- 最大行长度100字符
- 使用类型注解

文档
~~~~

- 所有公共API必须有docstrings
- 使用Google风格
- 包含参数和返回值说明
- 提供使用示例

测试
~~~~

- 新功能必须包含测试
- 测试文件放在 ``tests/`` 目录
- 命名规范: ``test_*.py``
- 运行 ``python run_test.py`` 确保所有测试通过

提交信息规范
~~~~~~~~~~~~

使用清晰的提交信息:

.. code-block:: text

   type(scope): subject
   
   body (optional)
   
   footer (optional)

类型:

- ``feat``: 新功能
- ``fix``: 修复bug
- ``docs``: 文档更新
- ``test``: 测试相关
- ``refactor``: 重构
- ``style``: 代码格式

示例:

.. code-block:: text

   feat(materials): add new porous material model
   
   Add Johnson-Lafarge model for materials with 
   non-circular pores.
   
   Closes #123

Pull Request流程
----------------

1. **创建PR**: 从feature分支创建PR到main分支
2. **填写模板**: 描述改动内容和动机
3. **代码审查**: 等待维护者审查
4. **修改反馈**: 根据反馈修改代码
5. **合并**: PR被批准后合并

代码审查清单
~~~~~~~~~~~~

- [ ] 代码符合项目风格
- [ ] 包含适当的测试
- [ ] 文档已更新
- [ ] 所有测试通过
- [ ] 没有引入回归问题

行为准则
--------

- 尊重他人
- 接受建设性批评
- 关注项目最佳利益
- 包容不同的观点和经验

许可证
------

贡献的代码将采用与项目相同的MIT许可证。

联系方式
--------

- 邮箱: shaoqiwu@outlook.com
- GitHub Issues: https://github.com/Shaoqigit/PyXfem/issues

我们非常期待你的贡献！🚀

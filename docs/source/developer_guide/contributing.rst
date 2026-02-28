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

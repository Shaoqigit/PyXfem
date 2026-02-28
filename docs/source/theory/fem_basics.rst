有限元方法基础
==============

本章节介绍有限元方法(Finite Element Method, FEM)的基本原理。

什么是有限元方法？
------------------

有限元方法是一种数值技术，用于求解偏微分方程(PDE)。其核心思想是:

1. **离散化**: 将连续域划分为有限个小单元
2. **近似**: 在每个单元上用简单函数近似解
3. **组装**: 将单元方程组合成全局方程组
4. **求解**: 求解线性方程组得到数值解

基本流程
--------

.. code-block:: text

   连续问题 → 弱形式 → 离散化 → 组装 → 求解 → 后处理

弱形式
------

考虑一般的边值问题:

.. math::

   Lu = f \quad \text{in } \Omega

其中 :math:`L` 是微分算子，:math:`u` 是未知函数，:math:`f` 是源项。

**强形式**: 要求解在每一点满足微分方程。

**弱形式**: 通过加权残差法，要求:

.. math::

   \int_\Omega (Lu - f)v \, d\Omega = 0, \quad \forall v \in V

其中 :math:`v` 是测试函数，:math:`V` 是适当的函数空间。

通过分部积分，可以降低对解的光滑性要求，得到弱形式。

离散化
------

**网格划分**

将计算域 :math:`\Omega` 划分为 :math:`N_e` 个单元:

.. math::

   \Omega \approx \Omega_h = \bigcup_{e=1}^{N_e} \Omega_e

常用单元类型:

- 1D: 线段
- 2D: 三角形、四边形
- 3D: 四面体、六面体

**基函数**

在每个单元上，解用基函数的线性组合近似:

.. math::

   u_h(x) = \sum_{j=1}^{N_n} u_j \phi_j(x)

其中 :math:`\phi_j` 是基函数，:math:`u_j` 是节点值。

常用基函数:

- **Lagrange多项式**: 节点在单元内部均匀分布
- **Lobatto多项式**: 节点包含端点，适合谱元法
- **Legendre多项式**: 正交多项式，数值稳定性好

单元方程
--------

将近似解代入弱形式，得到单元方程:

.. math::

   K^{(e)} u^{(e)} = F^{(e)}

其中:

- **刚度矩阵**: :math:`K^{(e)}_{ij} = \int_{\Omega_e} \nabla \phi_i \cdot \nabla \phi_j \, d\Omega`
- **载荷向量**: :math:`F^{(e)}_i = \int_{\Omega_e} f \phi_i \, d\Omega`

数值积分
~~~~~~~~

单元积分通常用高斯积分计算:

.. math::

   \int_{\Omega_e} g(x) \, d\Omega \approx \sum_{q=1}^{N_q} w_q g(x_q)

其中 :math:`w_q` 是权重，:math:`x_q` 是积分点。

全局组装
--------

将所有单元方程组装成全局系统:

.. math::

   K u = F

其中:

.. math::

   K = \bigoplus_{e=1}^{N_e} K^{(e)}, \quad F = \bigoplus_{e=1}^{N_e} F^{(e)}

组装过程通过**自由度映射**将单元自由度映射到全局自由度。

边界条件
--------

**本质边界条件 (Dirichlet)**

直接指定节点值:

.. math::

   u_i = g_i, \quad i \in \Gamma_D

实现方法:

- 直接代入法
- 大数法
- 消去法

**自然边界条件 (Neumann)**

出现在弱形式中:

.. math::

   \int_{\Gamma_N} g \phi_i \, d\Gamma

直接加到载荷向量中。

**混合边界条件 (Robin)**

.. math::

   \frac{\partial u}{\partial n} + \alpha u = \beta

同时影响刚度矩阵和载荷向量。

线性求解
--------

得到线性系统 :math:`Ku = F` 后，需要求解:

**直接求解器**

- LU分解
- Cholesky分解 (对称正定)

适用于中小规模问题 (:math:`N < 10^5`)

**迭代求解器**

- 共轭梯度法 (CG, 对称正定)
- GMRES (一般矩阵)
- BiCGSTAB (非对称)

适用于大规模稀疏系统。

误差分析
--------

**离散误差**

FEM解的误差与网格尺寸 :math:`h` 和多项式阶数 :math:`p` 有关:

.. math::

   \|u - u_h\|_{H^1} \leq C h^p |u|_{H^{p+1}}

这意味着:

- 细化网格 (:math:`h` 减小) 可以提高精度
- 提高阶数 (:math:`p` 增大) 可以提高精度

**收敛性检验**

通过网格细化检验收敛阶:

.. code-block:: python

   errors = []
   for h in [0.1, 0.05, 0.025, 0.0125]:
       uh = solve(h)
       error = compute_error(uh, u_exact)
       errors.append(error)
   
   # 计算收敛阶
   for i in range(1, len(errors)):
       rate = np.log(errors[i-1]/errors[i]) / np.log(2)
       print(f"收敛阶: {rate}")

PyAcoustiX中的实现
------------------

PyAcoustiX的FEM实现包括:

- **网格模块** (:mod:`SAcouS.Mesh`): 1D/2D/3D网格处理
- **基函数** (:mod:`SAcouS.acxfem.Basis`): Lagrange/Lobatto多项式
- **组装器** (:mod:`SAcouS.acxfem.Assembly`): 矩阵组装
- **求解器** (:mod:`SAcouS.acxfem.Solver`): 线性和迭代求解器
- **边界条件** (:mod:`SAcouS.acxfem.BCsImpose`): 边界条件施加

进一步阅读
----------

- Hughes, T.J.R. "The Finite Element Method"
- Zienkiewicz, O.C. and Taylor, R.L. "The Finite Element Method"
- Brenner, S.C. and Scott, L.R. "The Mathematical Theory of Finite Element Methods"

Helmholtz方程
=============

Helmholtz方程是声学仿真的核心控制方程，描述了时谐声波的传播。

波动方程
--------

从基本的流体动力学方程出发，可以得到声波方程:

.. math::

   \nabla^2 p - \frac{1}{c^2} \frac{\partial^2 p}{\partial t^2} = 0

其中:

- :math:`p`: 声压 [Pa]
- :math:`c`: 声速 [m/s]
- :math:`t`: 时间 [s]

时谐假设
--------

对于简谐振动，声压可以表示为:

.. math::

   p(\mathbf{x}, t) = \tilde{p}(\mathbf{x}) e^{j\omega t}

其中:

- :math:`\tilde{p}`: 复声压幅值
- :math:`\omega = 2\pi f`: 角频率 [rad/s]
- :math:`j`: 虚数单位

代入波动方程，得到**Helmholtz方程**:

.. math::

   \nabla^2 \tilde{p} + k^2 \tilde{p} = 0

其中 :math:`k = \omega/c = 2\pi/\lambda` 是**波数** [rad/m]。

弱形式
------

Helmholtz方程的弱形式:

.. math::

   \int_\Omega \nabla \tilde{p} \cdot \nabla v \, d\Omega - k^2 \int_\Omega \tilde{p} v \, d\Omega = \int_{\Gamma} \frac{\partial \tilde{p}}{\partial n} v \, d\Gamma

其中 :math:`v` 是测试函数。

边界条件
--------

**刚性壁 (Neumann)**

法向速度为零:

.. math::

   \frac{\partial p}{\partial n} = 0

对应于固体边界。

**压力边界 (Dirichlet)**

指定声压:

.. math::

   p = p_0

**速度边界 (Neumann)**

指定法向速度:

.. math::

   \frac{\partial p}{\partial n} = -j\omega\rho_0 v_n

其中 :math:`\rho_0` 是介质密度，:math:`v_n` 是法向速度。

**阻抗边界 (Robin)**

声压与速度的关系:

.. math::

   p = Z_s v_n = \frac{Z_s}{j\omega\rho_0} \frac{\partial p}{\partial n}

其中 :math:`Z_s` 是表面阻抗。

**完全匹配层 (PML)**

用于模拟无限域，通过在边界添加吸收层实现。

1D Helmholtz方程
----------------

简化形式:

.. math::

   \frac{d^2 p}{dx^2} + k^2 p = 0

通解:

.. math::

   p(x) = A e^{-jkx} + B e^{+jkx}

其中:

- :math:`A e^{-jkx}`: 右行波
- :math:`B e^{+jkx}`: 左行波

2D/3D Helmholtz方程
-------------------

**平面波解**

.. math::

   p(\mathbf{x}) = A e^{-j\mathbf{k} \cdot \mathbf{x}}

其中 :math:`\mathbf{k}` 是波矢量，:math:`|\mathbf{k}| = k`。

**球面波解 (3D)**

.. math::

   p(r) = \frac{A}{r} e^{-jkr}

表示从点源发出的波。

**柱面波解 (2D)**

.. math::

   p(r) = A H_0^{(2)}(kr)

其中 :math:`H_0^{(2)}` 是第二类Hankel函数。

数值求解
--------

**FEM离散**

使用Galerkin方法，得到离散方程:

.. math::

   (K - k^2 M) u = F

其中:

- **刚度矩阵**: :math:`K_{ij} = \int_\Omega \nabla \phi_i \cdot \nabla \phi_j \, d\Omega`
- **质量矩阵**: :math:`M_{ij} = \int_\Omega \phi_i \phi_j \, d\Omega`
- **载荷向量**: :math:`F_i = \int_\Gamma g \phi_i \, d\Gamma`

**波数分辨率**

为了准确捕捉波动现象，需要满足:

.. math::

   kh \leq C

其中 :math:`h` 是单元尺寸。经验法则:

- 线性单元 (:math:`p=1`): :math:`kh \leq 0.2` (每波长至少10个单元)
- 二次单元 (:math:`p=2`): :math:`kh \leq 0.5`
- 更高阶: 可以使用更大的 :math:`h`

**色散误差**

FEM解存在数值色散，即数值波数 :math:`k_h` 与真实波数 :math:`k` 不同:

.. math::

   k_h \approx k - \frac{k^3 h^2}{24}

高阶单元可以显著减小色散误差。

应用
----

Helmholtz方程应用于:

- **室内声学**: 房间声学响应
- **噪声控制**: 消声器、隔声设计
- **超声**: 医学成像、无损检测
- **水下声学**: 声纳、海洋声学
- **地震学**: 波传播分析

进一步阅读
----------

- Ihlenburg, F. "Finite Element Analysis of Acoustic Scattering"
- Wu, T.W. (Ed.) "Boundary Element Acoustics"
